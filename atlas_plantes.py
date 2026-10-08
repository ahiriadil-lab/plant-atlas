#!/usr/bin/env python3
"""Plant atlas: WFO taxa, LOTUS natural products, GBIF occurrences and environmental descriptors.

Publication edition of the pipeline used for release 20261005T141832_961774.

Pipeline
  1. World Flora Online Plant List (June 2026) -> queue of accepted species.
  2. LOTUS frozen metadata (April 2026) -> species-compound associations (exact WFO names).
  3. GBIF occurrence snapshot (2026-09-01, Parquet) -> reproducible per-species sample.
  4. SoilGrids 1 km (0-30 cm), Open-Meteo / NASA POWER 30-day weather, Copernicus DEM GLO-90.
  5. Validation, quarantine, SQLite release, summaries and map.

Usage
  python atlas_plantes.py --global-root GLOBAL --contact EMAIL              # collect a batch (repeat)
  python atlas_plantes.py --global-root GLOBAL --stage climate --contact EMAIL
  python atlas_plantes.py --global-root GLOBAL --atlas                      # assemble the release
"""

from pathlib import Path
import argparse, contextlib, datetime, hashlib, html, io, json, math, os, platform, re, shutil, sqlite3, time, zipfile
import email.utils, http.client, urllib.error, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
import numpy as np, pandas as pd, scipy, sklearn, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCHEMAS = {
    "taxa": ["taxon_id", "published_name", "accepted_name", "family", "genus", "source"],
    "compounds": ["compound_id", "compound_name", "smiles", "chemical_class", "source"],
    "samples": [
        "sample_id",
        "study_id",
        "site_id",
        "taxon_id",
        "organ",
        "stage",
        "cultivation",
        "latitude",
        "longitude",
        "geo_uncertainty_m",
        "date",
        "source",
    ],
    "chemistry": [
        "record_id",
        "analysis_id",
        "sample_id",
        "compound_id",
        "status",
        "searched",
        "method",
        "extraction",
        "lod",
        "value",
        "unit",
        "basis",
        "quantifier",
        "source",
    ],
    "occurrences": [
        "occurrence_id",
        "study_id",
        "site_id",
        "taxon_id",
        "latitude",
        "longitude",
        "geo_uncertainty_m",
        "date",
        "status",
        "survey_complete",
        "source",
    ],
    "environment": [
        "env_id",
        "point_id",
        "variable",
        "value",
        "unit",
        "origin",
        "source",
        "period_start",
        "period_end",
        "depth_top_cm",
        "depth_bottom_cm",
        "resolution_m",
    ],
    "associations": ["association_id", "taxon_id", "compound_id", "source"],
}

ENV_SPEC = {
    "soil_ph": ("pH", 0, 14),
    "soil_carbon": ("g/kg", 0, 1000),
    "soil_clay": ("%", 0, 100),
    "soil_nitrogen": ("g/kg", 0, 1000),
    "precip_30d": ("mm", 0, 5000),
    "temp_30d": ("degC", -80, 65),
    "altitude": ("m", -500, 9000),
    "temp_30d_cell": ("degC", -80, 65),
}

MAX_COARSE_UNCERTAINTY_M = 30000

POINT_UNCERTAINTY_M = 1000

LAPSE_RATE_C_PER_M = 0.0065

ELEVATION_QC_TOLERANCE_M = 300

NUM_ENV = [v for v in ENV_SPEC if v != "temp_30d_cell"]


# ---------------- Validation and features


def json_write(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def flag_bool(s):
    return s.astype(str).str.lower().isin(["true", "1", "yes", "oui"])


def event(log, stage, status, message, count=None):
    log.append(dict(stage=stage, status=status, message=message, count=count))


def nonempty(s):
    return s.notna() & s.astype(str).str.strip().ne("")


def validate_tables(raw, cfg, out, log):
    """Check schemas and units, drop duplicates and orphan links, and write clean and quarantine tables."""
    out = Path(out)
    clean = {}
    quarantines = {}
    keys = {
        "taxa": "taxon_id",
        "compounds": "compound_id",
        "samples": "sample_id",
        "chemistry": "record_id",
        "occurrences": "occurrence_id",
        "environment": "env_id",
        "associations": "association_id",
    }
    for name, cols in SCHEMAS.items():
        d = raw[name].copy()
        missing = set(cols) - set(d.columns)
        if missing:
            raise ValueError(f"{name}: missing columns {sorted(missing)}")
        for col in d.select_dtypes(include=["object", "string"]).columns:
            d[col] = d[col].map(lambda x: x.strip() if isinstance(x, str) else x)
        before = len(d)
        d = d.drop_duplicates()
        if len(d) < before:
            event(log, name, "CORRECTED", "Exact duplicates removed", before - len(d))
        key = keys[name]
        if not nonempty(d[key]).all():
            raise ValueError(f"{name}: empty identifier")
        if d[key].duplicated().any():
            raise ValueError(f"{name}: conflicting identifier, resolve before continuing")
        for col in [c for c in d if c.endswith("_id")]:
            d[col] = d[col].astype(str)
        if not nonempty(d["source"]).all():
            raise ValueError(f"{name}: missing required source")
        required = {
            "samples": ["study_id", "site_id", "taxon_id"],
            "occurrences": ["study_id", "site_id", "taxon_id"],
            "chemistry": ["analysis_id", "sample_id", "compound_id"],
            "environment": ["point_id"],
            "associations": ["taxon_id", "compound_id"],
        }.get(name, [])
        for col in required:
            if d[col].astype(str).isin(["", "nan", "None", "<NA>"]).any():
                raise ValueError(f"{name}: empty required key {col}")
        clean[name] = d.reset_index(drop=True)
    for name in ["samples", "occurrences"]:
        d = clean[name]
        for col in ["latitude", "longitude", "geo_uncertainty_m"]:
            d[col] = pd.to_numeric(d[col], errors="coerce")
        d["date"] = pd.to_datetime(d["date"], errors="coerce")
        d["geo_valid"] = (
            d.latitude.between(-90, 90)
            & d.longitude.between(-180, 180)
            & d.geo_uncertainty_m.between(0, cfg["max_geo_uncertainty_m"])
        )
        event(
            log,
            name,
            "REVIEW" if (~d.geo_valid).any() else "VALIDATED",
            "Inadmissible coordinates",
            int((~d.geo_valid).sum()),
        )
    c = clean["chemistry"]
    if c.duplicated(["analysis_id", "sample_id", "compound_id"], keep=False).any():
        raise ValueError(
            "Several results for one analysis–sample–compound: resolve or distinguish the analyses"
        )
    if not c.status.isin(["detected", "not_detected", "not_searched", "unknown"]).all():
        raise ValueError("Unknown chemical status")
    c["searched"] = flag_bool(c.searched)
    for col in ["value", "lod"]:
        c[col] = pd.to_numeric(c[col], errors="coerce")
    for unit, factor in [("ug/g", 0.001), ("µg/g", 0.001), ("mg/kg", 0.001)]:
        mask = c.unit.eq(unit)
        c.loc[mask, "value"] *= factor
        c.loc[mask, "lod"] *= factor
        c.loc[mask, "unit"] = "mg/g"
        if mask.any():
            event(log, "units", "CORRECTED", f"{unit} to mg/g", int(mask.sum()))
    valid_quant = (c.value >= 0) & c.quantifier.eq("exact") & c.status.eq("detected")
    c["quant_valid"] = valid_quant
    documented = nonempty(c.method) & c.lod.gt(0)
    c["label"] = np.where(
        c.status.eq("detected") & documented,
        1,
        np.where(c.status.eq("not_detected") & c.searched & documented, 0, np.nan),
    )
    o = clean["occurrences"]
    o["survey_complete"] = flag_bool(o.survey_complete)
    if not o.status.isin(["present", "survey_absent", "unknown", "background"]).all():
        raise ValueError("Unknown botanical status")
    o["label"] = np.where(
        o.status.eq("present"), 1, np.where(o.status.eq("survey_absent") & o.survey_complete, 0, np.nan)
    )
    relationships = [
        ("samples", "taxon_id", "taxa", "taxon_id"),
        ("occurrences", "taxon_id", "taxa", "taxon_id"),
        ("chemistry", "sample_id", "samples", "sample_id"),
        ("chemistry", "compound_id", "compounds", "compound_id"),
        ("associations", "taxon_id", "taxa", "taxon_id"),
        ("associations", "compound_id", "compounds", "compound_id"),
    ]
    for table, col, parent, pk in relationships:
        d = clean[table]
        bad = ~d[col].isin(clean[parent][pk])
        if bad.any():
            q = d.loc[bad].copy()
            q["reason"] = f"Orphan reference {col}"
            quarantines.setdefault(table, []).append(q)
            clean[table] = d.loc[~bad].copy()
            event(log, table, "QUARANTINE", f"Orphan references {col}", int(bad.sum()))
    for name, d in clean.items():
        d.to_csv(out / "clean" / (name + ".csv"), index=False)
    for name, parts in quarantines.items():
        pd.concat(parts).to_csv(out / "quarantine" / (name + ".csv"), index=False)
    json_write(out / "reports" / "quality.json", log)
    return clean


def enrich_environment(clean, cfg, out, log):
    """Select one value per point and variable, preferring measured over gridded values."""
    s = clean["samples"]
    o = clean["occurrences"]
    cols = ["point_id", "date", "latitude", "longitude", "geo_uncertainty_m"]
    points = pd.concat(
        [
            s.rename(columns={"sample_id": "point_id"})[cols],
            o.rename(columns={"occurrence_id": "point_id"})[cols],
        ],
        ignore_index=True,
    )
    if points.point_id.duplicated().any():
        raise ValueError("sample_id and occurrence_id must be globally distinct")
    e = clean["environment"].copy()
    for c in ["value", "resolution_m", "depth_top_cm", "depth_bottom_cm"]:
        e[c] = pd.to_numeric(e[c], errors="coerce")
    for c in ["period_start", "period_end"]:
        e[c] = pd.to_datetime(e[c], errors="coerce")
    e = e.merge(points, on="point_id", how="left", validate="many_to_one")
    coarse_cap = cfg.get("max_coarse_uncertainty_m", cfg["max_geo_uncertainty_m"])
    reasons = []
    for _, r in e.iterrows():
        reason = ""
        if r.variable not in ENV_SPEC:
            reason = "variable not configured"
        elif not (
            pd.notna(r.latitude)
            and -90 <= r.latitude <= 90
            and pd.notna(r.longitude)
            and -180 <= r.longitude <= 180
            and pd.notna(r.geo_uncertainty_m)
            and r.geo_uncertainty_m >= 0
        ):
            reason = "missing point or inadmissible geography"
        else:
            unit, low, high = ENV_SPEC[r.variable]
            if r.unit != unit or not pd.notna(r.value) or not low <= r.value <= high:
                reason = "invalid unit or value"
            elif (
                r.origin not in ["measured", "gridded"]
                or not isinstance(r.source, str)
                or not r.source.strip()
            ):
                reason = "missing provenance"
            elif not pd.notna(r.resolution_m) or not 0 <= r.resolution_m <= cfg["max_env_resolution_m"]:
                reason = "inadmissible resolution"
            elif r.geo_uncertainty_m > min(max(cfg["max_geo_uncertainty_m"], r.resolution_m), coarse_cap):
                reason = "point precision insufficient for this variable"
            elif r.variable.startswith("soil") and not (r.depth_top_cm == 0 and r.depth_bottom_cm == 30):
                reason = "depth other than 0–30 cm"
            elif r.variable in ["precip_30d", "temp_30d", "temp_30d_cell"]:
                if pd.isna(r.date) or pd.isna(r.period_start) or pd.isna(r.period_end):
                    reason = "missing period"
                elif (r.date - r.period_end).days != 0 or (r.period_end - r.period_start).days != 30:
                    reason = "window other than the 30 days before sampling"
        reasons.append(reason)
    e["reason"] = reasons
    rejected = e[e.reason.ne("")]
    valid = e[e.reason.eq("")].copy()
    valid["priority"] = valid.origin.map({"measured": 0, "gridded": 1})
    valid = valid[valid.priority.eq(valid.groupby(["point_id", "variable"]).priority.transform("min"))]
    ambiguous = valid.duplicated(["point_id", "variable"], keep=False)
    ambiguous_rows = valid.loc[ambiguous].copy()
    ambiguous_rows["reason"] = "several values of equal priority"
    valid = valid.loc[~ambiguous]
    pd.concat([rejected, ambiguous_rows]).to_csv(Path(out) / "quarantine" / "environment.csv", index=False)
    valid.to_csv(Path(out) / "clean" / "environment_selected.csv", index=False)
    wide = (
        valid.pivot(index="point_id", columns="variable", values="value")
        .reindex(columns=list(ENV_SPEC))
        .reset_index()
    )
    event(log, "environment", "VALIDATED", "Values retained", len(valid))
    event(log, "environment", "REVIEW", "Rejected or ambiguous values", len(rejected) + len(ambiguous_rows))
    return wide


def build_features(clean, env):
    """Join taxa, compounds and environmental descriptors into chemical and geographic views."""
    tax = clean["taxa"][["taxon_id", "accepted_name", "family", "genus"]]
    mol = clean["compounds"][["compound_id", "chemical_class", "compound_name", "smiles"]]
    s = clean["samples"].merge(tax, on="taxon_id", validate="many_to_one")
    s = s.merge(
        env.rename(columns={"point_id": "sample_id"}), on="sample_id", how="left", validate="one_to_one"
    )
    chem = clean["chemistry"].merge(s, on="sample_id", validate="many_to_one", suffixes=("", "_sample"))
    chem = chem.merge(mol, on="compound_id", validate="many_to_one")
    geo = clean["occurrences"].merge(tax, on="taxon_id", validate="many_to_one")
    geo = geo.merge(
        env.rename(columns={"point_id": "occurrence_id"}),
        on="occurrence_id",
        how="left",
        validate="one_to_one",
    )
    return chem, geo


ATLAS_MAP_PAGE = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Atlas map</title><style>
:root{--bg:#f6f8f5;--fg:#18241d;--mut:#5b6b61;--card:#fff;--line:#dfe6e0;--acc:#2f6f4e;--pt:#d9660a}
@media (prefers-color-scheme:dark){:root{--bg:#111a15;--fg:#e8efe9;--mut:#a3b3a8;--card:#18241d;--line:#2a3a31;--acc:#7cc39b;--pt:#ffb454}}
*{box-sizing:border-box}body{margin:0;padding:16px;background:var(--bg);color:var(--fg);font:15px/1.45 system-ui,sans-serif}main{max-width:1200px;margin:auto}
h1{margin:0 0 4px;font-size:22px}p{color:var(--mut);margin:4px 0 12px}.bar{display:flex;flex-wrap:wrap;gap:8px;align-items:flex-start;margin-bottom:10px}
label{display:block;font-size:12px;color:var(--mut);margin-bottom:2px}input,select,button{font:inherit;padding:8px 10px;border:1px solid var(--line);border-radius:8px;background:var(--card);color:var(--fg)}
button{cursor:pointer}button:hover{border-color:var(--acc)}#q{width:260px}#res{position:absolute;z-index:5;background:var(--card);border:1px solid var(--line);border-radius:8px;max-height:300px;overflow:auto;width:340px}
#res button{display:block;width:100%;text-align:left;border:0;border-radius:0;border-bottom:1px solid var(--line)}#res small{color:var(--mut)}.w{position:relative}
canvas{width:100%;height:520px;display:block;background:var(--card);border:1px solid var(--line);border-radius:10px;touch-action:none;cursor:grab}
#info{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px;margin-top:10px;min-height:64px}#hover{color:var(--mut);font-size:13px;min-height:20px;margin-top:6px}
#legend{display:flex;align-items:center;gap:8px;font-size:13px;color:var(--mut);margin-top:8px}#ramp{flex:0 0 160px;height:10px;border-radius:5px;border:1px solid var(--line)}
</style><main><h1>Plant atlas map</h1>
<p>Overview: number of distinct plant species observed in each cell (logarithmic scale). Choose a species to see all its precise
occurrences (uncertainty ≤ 1 km). Scroll or use the buttons to zoom, drag to pan. No external service: the map works offline.</p>
<div class="bar"><div class="w"><label for="q">Species</label><input id="q" type="search" placeholder="Scientific name…" autocomplete="off"><div id="res" hidden></div></div>
<div><label for="fam">Family</label><select id="fam"></select></div>
<div><label>&nbsp;</label><button id="zi" aria-label="Zoom in">＋</button> <button id="zo" aria-label="Zoom out">－</button> <button id="rs">World view</button></div></div>
<canvas id="map" role="img" aria-label="Density map of plant occurrences"></canvas>
<div id="legend"><span id="lmin">1</span><span id="ramp"></span><span id="lmax"></span><span>species per cell</span></div><div id="hover" aria-live="off"></div>
<div id="info" aria-live="polite">Choose a species in the search box.</div></main>
<script>const D=__DATA__;
const N=D.taxa.length,LAT=new Float32Array(D.lat),LON=new Float32Array(D.lon),OFF=Int32Array.from(D.off),NAMES=D.taxa.map(t=>t[0].toLowerCase());
const cv=document.getElementById('map'),ctx=cv.getContext('2d'),q=document.getElementById('q'),res=document.getElementById('res'),famSel=document.getElementById('fam'),info=document.getElementById('info'),hover=document.getElementById('hover');
const dark=matchMedia('(prefers-color-scheme:dark)').matches,css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
let z=1,tx=0,ty=0,sel=-1,fam=-1,cells=null,dpr=1;
const shade=t=>`hsl(160,55%,${dark?18+47*t:90-62*t}%)`;
document.getElementById('ramp').style.background=`linear-gradient(90deg,${shade(0)},${shade(.5)},${shade(1)})`;
famSel.append(new Option('All families',-1));D.fam.forEach((f,i)=>famSel.append(new Option(f,i)));
function size(){dpr=window.devicePixelRatio||1;cv.width=Math.round(cv.clientWidth*dpr);cv.height=Math.round(520*dpr);reset()}
function reset(redraw=true){z=1;tx=0;ty=-(90-75)/180*(cv.width/2);if(redraw)draw()}
function fit(a,b,c,d){const w=Math.max(b-a,1),h=Math.max(d-c,1);z=Math.min(300,.8*360/w,.8*cv.height*360/(h*cv.width));
 tx=cv.width/2-((a+b)/2+180)/360*cv.width*z;ty=cv.height/2-(90-(c+d)/2)/180*(cv.width/2)*z}
function zoomAt(f,x,y){const nz=Math.min(400,Math.max(1,z*f));f=nz/z;tx=x-(x-tx)*f;ty=y-(y-ty)*f;z=nz;draw()}
function draw(){
 const cs=Math.round(8*dpr),cols=Math.ceil(cv.width/cs),rows=Math.ceil(cv.height/cs),cnt=new Uint32Array(cols*rows),dist=new Uint32Array(cols*rows),last=new Int32Array(cols*rows).fill(-1);
 const sx=cv.width/360*z,sy=cv.width/360*z,ox=180*sx+tx,oy=90*sy+ty;let maxd=1;
 for(let i=0;i<N;i++){if(fam>=0&&D.taxa[i][1]!==fam)continue;
  for(let j=OFF[i];j<OFF[i+1];j++){const x=LON[j]*sx+ox,y=oy-LAT[j]*sy;if(x<0||y<0||x>=cv.width||y>=cv.height)continue;
   const k=((y/cs)|0)*cols+((x/cs)|0);cnt[k]++;if(last[k]!==i){last[k]=i;if(++dist[k]>maxd)maxd=dist[k]}}}
 cells={cnt,dist,cols,cs};ctx.clearRect(0,0,cv.width,cv.height);
 ctx.strokeStyle=css('--line');ctx.lineWidth=1;ctx.beginPath();
 for(let lo=-180;lo<=180;lo+=30){const x=lo*sx+ox;ctx.moveTo(x,0);ctx.lineTo(x,cv.height)}
 for(let la=-60;la<=90;la+=30){const y=oy-la*sy;ctx.moveTo(0,y);ctx.lineTo(cv.width,y)}ctx.stroke();
 ctx.globalAlpha=sel>=0?.3:1;
 for(let k=0;k<dist.length;k++)if(dist[k]){ctx.fillStyle=shade(Math.log(1+dist[k])/Math.log(1+maxd));ctx.fillRect((k%cols)*cs,((k/cols)|0)*cs,cs-1,cs-1)}
 ctx.globalAlpha=1;
 if(sel>=0){ctx.fillStyle=css('--pt');ctx.strokeStyle=css('--bg');ctx.lineWidth=1.5*dpr;
  for(let j=OFF[sel];j<OFF[sel+1];j++){ctx.beginPath();ctx.arc(LON[j]*sx+ox,oy-LAT[j]*sy,4*dpr,0,7);ctx.fill();ctx.stroke()}}
 document.getElementById('lmax').textContent=maxd.toLocaleString('en-US');
 document.getElementById('lmin').textContent='1'}
function choose(i){sel=i;res.hidden=true;const t=D.taxa[i],n=OFF[i+1]-OFF[i];q.value=t[0];
 if(n){let a=180,b=-180,c=90,d=-90;for(let j=OFF[i];j<OFF[i+1];j++){a=Math.min(a,LON[j]);b=Math.max(b,LON[j]);c=Math.min(c,LAT[j]);d=Math.max(d,LAT[j])}fit(a,b,c,d)}
 const f=(v,u)=>v===null?'—':v.toLocaleString('en-US')+u,line=[`Family: ${D.fam[t[1]]}`,`${t[2].toLocaleString('en-US')} known molecules (literature, not local evidence)`,
  `${n} precise occurrences shown out of ${t[3].toLocaleString('en-US')} retained`,`Medians — elevation: ${f(t[4],' m')} · 30-day temperature: ${f(t[5],' °C')} · soil pH: ${f(t[6],'')}`];
 info.replaceChildren();const h=document.createElement('b');h.textContent=t[0];info.append(h);line.forEach(l=>{const d=document.createElement('div');d.textContent=l;info.append(d)});
 if(!n){const d=document.createElement('div');d.textContent="No precise occurrence (≤ 1 km) for this species.";info.append(d)}draw()}
q.addEventListener('input',()=>{const v=q.value.trim().toLowerCase();res.replaceChildren();if(v.length<2){res.hidden=true;return}
 let shown=0;for(let i=0;i<N&&shown<12;i++)if(NAMES[i].includes(v)){const b=document.createElement('button'),s=document.createElement('small');b.type='button';b.textContent=D.taxa[i][0]+' ';
  s.textContent=`${D.fam[D.taxa[i][1]]} · ${D.taxa[i][2]} molecules`;b.append(s);b.onclick=()=>choose(i);res.append(b);shown++}res.hidden=!shown});
q.addEventListener('keydown',e=>{if(e.key==='Enter'&&res.firstChild)res.firstChild.click();if(e.key==='Escape')res.hidden=true});
famSel.onchange=()=>{fam=+famSel.value;sel=-1;info.textContent=fam<0?'Choose a species in the search box.':`Family: ${D.fam[fam]}`;draw()};
document.getElementById('zi').onclick=()=>zoomAt(1.6,cv.width/2,cv.height/2);document.getElementById('zo').onclick=()=>zoomAt(1/1.6,cv.width/2,cv.height/2);
document.getElementById('rs').onclick=()=>{sel=-1;q.value='';info.textContent='Choose a species in the search box.';reset()};
cv.addEventListener('wheel',e=>{e.preventDefault();const r=cv.getBoundingClientRect();zoomAt(e.deltaY<0?1.25:.8,(e.clientX-r.left)*dpr,(e.clientY-r.top)*dpr)},{passive:false});
let drag=null;cv.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY,tx,ty];cv.setPointerCapture(e.pointerId);cv.style.cursor='grabbing'});
cv.addEventListener('pointerup',()=>{drag=null;cv.style.cursor='grab'});
cv.addEventListener('pointermove',e=>{if(drag){tx=drag[2]+(e.clientX-drag[0])*dpr;ty=drag[3]+(e.clientY-drag[1])*dpr;draw();return}
 const r=cv.getBoundingClientRect(),x=(e.clientX-r.left)*dpr,y=(e.clientY-r.top)*dpr;if(!cells)return;const k=((y/cells.cs)|0)*cells.cols+((x/cells.cs)|0);
 hover.textContent=cells.dist[k]?`${cells.dist[k].toLocaleString('en-US')} species · ${cells.cnt[k].toLocaleString('en-US')} occurrences in this cell`:''});
addEventListener('resize',size);size();
const m=decodeURIComponent(location.hash.slice(1)).match(/^(plant|family)=(.*)$/);   // direct links: #plant=Thymus vulgaris, #family=Lamiaceae
if(m){const v=m[2].toLowerCase();if(m[1]==='plant'){const i=NAMES.indexOf(v);if(i>=0)choose(i)}else{const k=D.fam.findIndex(f=>f.toLowerCase()===v);if(k>=0){famSel.value=k;famSel.onchange()}}}
</script></html>"""


def write_atlas_map(sp, geo, path):
    """Self-contained HTML map: species richness per cell, or the precise points of one species."""
    sp = sp.reset_index(drop=True)
    index = {t: i for i, t in enumerate(sp.taxon_id)}
    pts = geo[geo.geo_valid & geo.status.eq("present") & geo.taxon_id.isin(index)][
        ["taxon_id", "latitude", "longitude"]
    ].copy()
    pts["i"] = pts.taxon_id.map(index)
    pts = pts.sort_values("i", kind="stable")
    off = (
        np.concatenate([[0], np.cumsum(np.bincount(pts.i.to_numpy(), minlength=len(sp)))])
        .astype(int)
        .tolist()
    )
    family = sp.family.fillna("unknown")
    fams = family.value_counts().index.tolist()
    fi = {f: k for k, f in enumerate(fams)}
    med = lambda c: [
        None if pd.isna(v) else round(float(v), 1) for v in (sp[c] if c in sp else [np.nan] * len(sp))
    ]
    alt, temp, ph = med("median_altitude"), med("median_temp_30d"), med("median_soil_ph")
    taxa = [
        [str(n), fi[f], int(m), int(o), a, t, h]
        for n, f, m, o, a, t, h in zip(sp.accepted_name, family, sp.molecules, sp.observations, alt, temp, ph)
    ]
    data = dict(
        taxa=taxa, fam=fams, lat=pts.latitude.round(3).tolist(), lon=pts.longitude.round(3).tolist(), off=off
    )
    page = ATLAS_MAP_PAGE.replace(
        "__DATA__", json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
    )
    Path(path).write_text(page, encoding="utf-8")
    return len(pts)


def export_graph(clean, out):
    """Taxon-compound graph as node and edge tables."""
    edges = []
    for _, r in clean["taxa"].iterrows():
        edges.append(("taxon:" + r.taxon_id, "belongs_to", "genus:" + str(r.genus), r.source))
        edges.append(("genus:" + str(r.genus), "belongs_to", "family:" + str(r.family), r.source))
    for _, r in clean["associations"].iterrows():
        edges.extend(
            [
                ("association:" + r.association_id, "taxon", "taxon:" + r.taxon_id, r.source),
                ("association:" + r.association_id, "compound", "compound:" + r.compound_id, r.source),
            ]
        )
    for _, r in clean["samples"].iterrows():
        edges.extend(
            [
                ("sample:" + r.sample_id, "taxon", "taxon:" + r.taxon_id, r.source),
                ("sample:" + r.sample_id, "site", "site:" + r.site_id, r.source),
            ]
        )
    for _, r in clean["chemistry"].iterrows():
        edges.extend(
            [
                ("result:" + r.record_id, "sample", "sample:" + r.sample_id, r.source),
                ("result:" + r.record_id, r.status, "compound:" + r.compound_id, r.source),
            ]
        )
    edge = pd.DataFrame(edges, columns=["subject", "predicate", "object", "source"]).drop_duplicates()
    edge.to_csv(Path(out) / "graph_edges.csv", index=False)
    nodes = pd.unique(pd.concat([edge.subject, edge.object]))
    pd.DataFrame({"node_id": nodes}).to_csv(Path(out) / "graph_nodes.csv", index=False)
    return edge


# ---------------- Public services

CONTACT = os.environ.get("ATLAS_CONTACT", "")

MAX_RETRY_WAIT = 120


def digest(s):
    return hashlib.sha256(str(s).encode()).hexdigest()


def db_connect(path, **kw):
    """SQLite connection with a long busy timeout (collection and climate passes share the database)."""
    return sqlite3.connect(path, **{"timeout": 120, **kw})


class NetworkDown(RuntimeError):
    """No network response; the species stays queued instead of being marked partial."""


def check_payload(data, accept):
    """Reject error pages served with HTTP 200."""
    if "json" in accept:
        json.loads(data)
    elif "tiff" in accept and data[:4] not in (b"II*\x00", b"MM\x00*"):
        raise ValueError("Non-GeoTIFF response: " + data[:120].decode("utf-8", "replace"))


def retry_after_seconds(value, default=10):
    """Parse a Retry-After header (seconds or HTTP date)."""
    if not value:
        return default
    if str(value).strip().isdigit():
        return int(value)
    try:
        return max(
            0,
            (
                email.utils.parsedate_to_datetime(value) - datetime.datetime.now(datetime.timezone.utc)
            ).total_seconds(),
        )
    except (TypeError, ValueError):
        return default


class PublicClient:
    """HTTP client with a raw response cache, per-host pacing and a request log."""

    def __init__(self, root, offline=False, timeout=25, retries=2, contact=None):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        contact = contact or CONTACT
        self.offline = offline
        self.timeout = timeout
        self.retries = retries
        self.last = {}
        self.events = []
        self.interval = {}
        self.disabled = set()
        self.base = {}
        self.backoff = 5.0
        self.ua = "AtlasPlantesResearch/3.1 (python-urllib" + ("; " + contact if contact else "") + ")"

    def log(self, source, status, message, **kwargs):
        row = dict(
            time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            source=source,
            status=status,
            message=message,
            **kwargs,
        )
        self.events.append(row)
        with (self.root / "requests.jsonl").open("a", encoding="utf-8") as f:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")

    def get(
        self, base, params=None, source="http", accept="application/json", retries=None, fail_status="FAILED"
    ):
        retries = self.retries if retries is None else retries
        url = base + ("?" + urllib.parse.urlencode(params, doseq=True) if params else "")
        key = digest(url)
        p = self.root / (key + ".bin")
        meta = self.root / (key + ".json")
        if p.exists() and meta.exists():
            data = p.read_bytes()
            m = json.loads(meta.read_text())
            if hashlib.sha256(data).hexdigest() != m.get("sha256"):
                raise ValueError("Corrupted cache: " + str(p))
            try:
                check_payload(data, accept)
                self.log(source, "CACHE", "Raw response reused", url=url)
                return data
            except ValueError:
                p.unlink()
                meta.unlink()
                self.log(source, "CACHE_INVALID", "Invalid cached response removed", url=url)
        if self.offline:
            raise RuntimeError("Response missing from offline cache: " + source)
        host = urllib.parse.urlparse(base).netloc
        attempts = 1 if retries == 0 else max(retries, 4) + 1
        for attempt in range(attempts):
            delay = max(0, self.interval.get(host, 1.1) - (time.monotonic() - self.last.get(host, 0)))
            if delay:
                time.sleep(delay)
            self.last[host] = time.monotonic()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": self.ua, "Accept": accept})
                with urllib.request.urlopen(req, timeout=self.timeout) as response:
                    data = response.read(25_000_001)
                    if len(data) > 25_000_000:
                        raise ValueError("Response >25 MB: use a bulk export")
                    content_type = response.headers.get("Content-Type", "")
                    m = {
                        "url": url,
                        "retrieved_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        "content_type": content_type,
                        "sha256": hashlib.sha256(data).hexdigest(),
                    }
                check_payload(data, accept)
                p.write_bytes(data)
                meta.write_text(json.dumps(m, indent=2), encoding="utf-8")
                if host in self.interval:
                    self.interval[host] = max(self.base.get(host, 1.1), self.interval[host] * 0.97)
                self.log(source, "OK", "Downloaded", url=url, bytes=len(data))
                return data
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last = exc
                if isinstance(exc, urllib.error.HTTPError) and exc.code not in [429, 500, 502, 503, 504]:
                    break
                wait = min(2**attempt, 4)
                if isinstance(exc, urllib.error.HTTPError) and exc.code in [429, 503]:
                    asked = retry_after_seconds(exc.headers.get("Retry-After"), default=0)
                    wait = min(max(asked, self.backoff * 2**attempt), MAX_RETRY_WAIT)
                    self.interval[host] = min(max(self.interval.get(host, 1.1) * 2, asked), MAX_RETRY_WAIT)
                    self.log(source, "THROTTLED", f"HTTP {exc.code}; waiting {wait:.0f} s", url=url)
                if attempt < attempts - 1:
                    time.sleep(wait)
                last = exc
            except Exception as exc:
                last = exc
                break
        last = locals().get("last")
        if (
            isinstance(last, urllib.error.URLError)
            and not isinstance(last, urllib.error.HTTPError)
            and isinstance(last.reason, OSError)
            and not isinstance(last.reason, TimeoutError)
        ):
            self.log(source, "OFFLINE", str(last), url=url)
            raise NetworkDown(source + ": Internet connection unavailable (" + str(last.reason) + ")")
        msg = str(last or "Request refused by the service")
        self.log(source, fail_status, msg, url=url)
        raise RuntimeError(source + ": " + msg)

    def json(self, base, params=None, source="http", accept="application/json", **kw):
        return json.loads(self.get(base, params, source, accept, **kw))


_CLIENTS = {}


def public_client(config):
    """One shared client per cache directory."""
    key = (str(config["cache_dir"]), config.get("offline", False), config.get("timeout", 25))
    if key not in _CLIENTS:
        _CLIENTS[key] = PublicClient(config["cache_dir"], offline=key[1], timeout=key[2])
    client = _CLIENTS[key]
    client.events = []
    return client


def parse_gbif_occurrence(r, taxon):
    """GBIF record -> occurrence row."""
    lat = r.get("decimalLatitude")
    lon = r.get("decimalLongitude")
    if lat is None or lon is None or not (-90 <= lat <= 90 and -180 <= lon <= 180):
        return None
    date = str(r.get("eventDate", ""))
    date = date[:10] if re.fullmatch(r"\d{4}-\d{2}-\d{2}(T.*)?", date) else ""
    grid = f"{math.floor(lat / 0.25)}_{math.floor(lon / 0.25)}"
    return dict(
        occurrence_id="GBIF:" + str(r["key"]),
        study_id="GBIFDATA:" + str(r.get("datasetKey", "unknown")),
        site_id="BLOCK025:" + grid,
        taxon_id=taxon["taxon_id"],
        latitude=lat,
        longitude=lon,
        geo_uncertainty_m=r.get("coordinateUncertaintyInMeters", np.nan),
        date=date,
        status="present",
        survey_complete=False,
        source="https://www.gbif.org/occurrence/" + str(r["key"]),
        license=r.get("license", ""),
        dataset_key=r.get("datasetKey", ""),
        basis_of_record=r.get("basisOfRecord", ""),
        establishment_means=r.get("establishmentMeans", ""),
        event_date_original=r.get("eventDate", ""),
        occurrence_issues=json.dumps(r.get("issues", [])),
        reported_elevation_m=r.get("elevation", np.nan),
        reported_elevation_accuracy_m=r.get("elevationAccuracy", np.nan),
    )


def environment_row(
    point,
    variable,
    value,
    unit,
    source,
    start="",
    end="",
    depth_top=np.nan,
    depth_bottom=np.nan,
    resolution=250,
    **extra,
):
    """Build one environment row with its source, period, depth and resolution."""
    return dict(
        env_id="ENV:"
        + digest(point["point_id"] + "|" + variable + "|" + source + "|" + str(start) + "|" + str(end))[:24],
        point_id=point["point_id"],
        variable=variable,
        value=value,
        unit=unit,
        origin="gridded",
        source=source,
        period_start=start,
        period_end=end,
        depth_top_cm=depth_top,
        depth_bottom_cm=depth_bottom,
        resolution_m=resolution,
        **extra,
    )


def fetch_weather_nasa(client, point):
    """NASA POWER fallback: 30-day temperature and precipitation for the grid cell."""
    day = pd.to_datetime(point.get("date"), errors="coerce")
    if pd.isna(day) or day < pd.Timestamp("1981-02-01") or day.date() >= datetime.date.today():
        return []
    start = day - pd.Timedelta(days=30)
    last = day - pd.Timedelta(days=1)
    lat = round(float(point["latitude"]) / 0.5) * 0.5
    lon = round(float(point["longitude"]) / 0.625) * 0.625
    lon = max(-180, min(179.375, lon))
    p = {
        "parameters": "T2M,PRECTOTCORR",
        "community": "AG",
        "longitude": lon,
        "latitude": lat,
        "start": start.strftime("%Y%m%d"),
        "end": last.strftime("%Y%m%d"),
        "format": "JSON",
        "time-standard": "UTC",
    }
    d = client.json("https://power.larc.nasa.gov/api/temporal/daily/point", p, "NASA POWER")
    dates = pd.date_range(start, last).strftime("%Y%m%d").tolist()
    parameters = d.get("properties", {}).get("parameter", {})
    fill = d.get("header", {}).get("fill_value", -999)
    coords = (d.get("geometry") or {}).get("coordinates") or []
    grid_elevation = float(coords[2]) if len(coords) > 2 and coords[2] is not None else np.nan
    rows = []
    for key, var, unit, op in [
        ("T2M", "temp_30d_cell", "degC", np.mean),
        ("PRECTOTCORR", "precip_30d", "mm", np.sum),
    ]:
        series = parameters.get(key, {})
        vals = [series.get(date) for date in dates]
        if any(v is None or v == fill or not np.isfinite(v) for v in vals):
            client.log(
                "NASA POWER", "SKIP", "Incomplete 30-day series", point=point["point_id"], variable=var
            )
            continue
        actual_unit = d.get("parameters", {}).get(key, {}).get("units")
        if actual_unit not in ({"C", "degC"} if key == "T2M" else {"mm/day"}):
            raise ValueError("Unexpected NASA unit: " + str(actual_unit))
        rows.append(
            environment_row(
                point,
                var,
                float(op(vals)),
                unit,
                "https://power.larc.nasa.gov/ | MERRA2 daily UTC",
                start.date().isoformat(),
                day.date().isoformat(),
                resolution=max(70000, float(point.get("geo_uncertainty_m") or 0)),
                aggregation="30 days [start,end): " + ("mean" if key == "T2M" else "sum"),
                grid_lat=lat,
                grid_lon=lon,
                grid_elevation_m=grid_elevation,
                source_start=d.get("header", {}).get("start"),
                source_end=d.get("header", {}).get("end"),
            )
        )
    return rows


def fetch_weather(client, point):
    """Open-Meteo archive: 30-day temperature and precipitation before the occurrence date."""
    day = pd.to_datetime(point.get("date"), errors="coerce")
    if pd.isna(day):
        return []
    in_window = pd.Timestamp("1950-02-01") <= day <= pd.Timestamp.now() - pd.Timedelta(days=7)
    if "open-meteo-archive" not in client.disabled and in_window:
        for d in (client.base, client.interval):
            d.setdefault("archive-api.open-meteo.com", 1.6)
        start = day - pd.Timedelta(days=30)
        last = day - pd.Timedelta(days=1)
        try:
            d = client.json(
                "https://archive-api.open-meteo.com/v1/archive",
                {
                    "latitude": f"{float(point['latitude']):.5f}",
                    "longitude": f"{float(point['longitude']):.5f}",
                    "start_date": start.date().isoformat(),
                    "end_date": last.date().isoformat(),
                    "daily": "temperature_2m_mean,precipitation_sum",
                    "timezone": "UTC",
                },
                "Open-Meteo archive",
                retries=0,
                fail_status="FALLBACK",
            )
            daily = d.get("daily", {})
            units = d.get("daily_units", {})
            if units.get("temperature_2m_mean") != "°C" or units.get("precipitation_sum") != "mm":
                raise ValueError("Unexpected Open-Meteo units: " + str(units))
            rows = []
            for key, var, op in [
                ("temperature_2m_mean", "temp_30d", np.mean),
                ("precipitation_sum", "precip_30d", np.sum),
            ]:
                vals = daily.get(key) or []
                if len(vals) != 30 or any(v is None for v in vals):
                    client.log(
                        "Open-Meteo archive",
                        "SKIP",
                        "Incomplete 30-day series",
                        point=point["point_id"],
                        variable=var,
                    )
                    continue
                rows.append(
                    environment_row(
                        point,
                        var,
                        float(op(vals)),
                        "degC" if var == "temp_30d" else "mm",
                        "https://open-meteo.com/ | best_match ERA5-Land/ERA5 daily UTC",
                        start.date().isoformat(),
                        day.date().isoformat(),
                        resolution=max(28000, float(point.get("geo_uncertainty_m") or 0)),
                        aggregation="30 days [start,end): " + ("mean" if var == "temp_30d" else "sum"),
                        grid_elevation_m=d.get("elevation"),
                        license="ERA5 / ERA5-Land (Copernicus) ; Open-Meteo CC-BY-4.0, free non-commercial use",
                    )
                )
            elevation = pd.to_numeric(d.get("elevation"), errors="coerce")
            if pd.notna(elevation) and float(point.get("geo_uncertainty_m") or 0) <= POINT_UNCERTAINTY_M:
                rows.append(
                    environment_row(
                        point,
                        "altitude",
                        float(elevation),
                        "m",
                        "https://open-meteo.com/ | Copernicus DEM GLO-90 (weather response)",
                        resolution=90,
                        license="Copernicus DEM GLO-90 free licence ; Open-Meteo CC-BY-4.0",
                    )
                )
            client.om_failures = 0
            return rows
        except (RuntimeError, ValueError) as exc:
            client.om_failures = getattr(client, "om_failures", 0) + 1
            if "429" in str(exc) or client.om_failures >= 3:
                client.disabled.add("open-meteo-archive")
                print("Open-Meteo archive unavailable for this session:", str(exc)[:120], flush=True)
    if in_window and not getattr(client, "nasa_fallback", True):
        client.om_failed = True
        return []
    return fetch_weather_nasa(client, point)


SOIL_PROPS = {
    "phh2o": ("soil_ph", 10, "pH"),
    "soc": ("soil_carbon", 10, "g/kg"),
    "clay": ("soil_clay", 10, "%"),
    "nitrogen": ("soil_nitrogen", 100, "g/kg"),
}


def decode_soil_tiff(data, lat, lon):
    """Read the pixel value at a coordinate from a single-band GeoTIFF."""
    from PIL import Image

    im = Image.open(io.BytesIO(data))
    tags = im.tag_v2
    tie = tags.get(33922)
    scale = tags.get(33550)
    if not tie or not scale:
        raise ValueError("Missing GeoTIFF georeferencing")
    keys = tags.get(34735, ())
    geographic = [keys[i + 3] for i in range(4, len(keys), 4) if keys[i] == 2048 and keys[i + 1] == 0]
    if geographic != [4326]:
        raise ValueError("The WCS must return an EPSG:4326 GeoTIFF")
    col = int(math.floor((lon - tie[3]) / scale[0] + tie[0]))
    row = int(math.floor((tie[4] - lat) / scale[1] + tie[1]))
    if not (0 <= col < im.width and 0 <= row < im.height):
        raise ValueError("Point outside the returned raster")
    value = float(im.getpixel((col, row)))
    nodata = tags.get(42113)
    if value < 0 or (nodata is not None and value == float(str(nodata).strip("\x00"))):
        return None
    return value


def fetch_soil(client, point, properties=None):
    """SoilGrids 250 m WCS fallback for one point (0-30 cm, thickness-weighted)."""
    lat = float(point["latitude"])
    lon = float(point["longitude"])
    rows = []
    for prop in properties or list(SOIL_PROPS):
        variable, divisor, unit = SOIL_PROPS[prop]
        parts = []
        for depth, weight in [("0-5cm", 5), ("5-15cm", 10), ("15-30cm", 15)]:
            params = [
                ("map", "/map/" + prop + ".map"),
                ("SERVICE", "WCS"),
                ("VERSION", "2.0.1"),
                ("REQUEST", "GetCoverage"),
                ("COVERAGEID", prop + "_" + depth + "_mean"),
                ("FORMAT", "GEOTIFF_INT16"),
                ("SUBSET", f"X({lon - 0.001},{lon + 0.001})"),
                ("SUBSET", f"Y({lat - 0.001},{lat + 0.001})"),
                ("SUBSETTINGCRS", "http://www.opengis.net/def/crs/EPSG/0/4326"),
                ("OUTPUTCRS", "http://www.opengis.net/def/crs/EPSG/0/4326"),
            ]
            data = client.get(
                "https://maps.isric.org/mapserv", params, "SoilGrids WCS", accept="image/tiff, */*"
            )
            v = decode_soil_tiff(data, lat, lon)
            if v is None:
                parts = []
                break
            parts.append((v / divisor, weight))
        if len(parts) == 3:
            val = sum(v * w for v, w in parts) / 30
            rows.append(
                environment_row(
                    point,
                    variable,
                    val,
                    unit,
                    "https://maps.isric.org/ | SoilGrids " + prop + " mean",
                    depth_top=0,
                    depth_bottom=30,
                    resolution=250,
                    aggregation="thickness-weighted profile descriptor; NOT bulk measurement",
                    depth_values=json.dumps(parts),
                    license="CC-BY-4.0",
                )
            )
    return rows


def admissible_points(occurrences, max_uncertainty_m):
    """Occurrences with a known coordinate uncertainty not above max_uncertainty_m."""
    points = []
    for r in occurrences:
        uncertainty = pd.to_numeric(r.get("geo_uncertainty_m"), errors="coerce")
        if pd.notna(uncertainty) and 0 <= uncertainty <= max_uncertainty_m:
            points.append(
                dict(
                    point_id=r["occurrence_id"],
                    latitude=r["latitude"],
                    longitude=r["longitude"],
                    date=r["date"],
                    geo_uncertainty_m=float(uncertainty),
                )
            )
    return points


def climate_points(occurrences, config):
    """Weather point per species: most recent precise point, else a coarse one (<= 30 km)."""
    points = admissible_points(occurrences, config["max_geo_uncertainty_m"]) or admissible_points(
        occurrences, config.get("max_coarse_uncertainty_m", config["max_geo_uncertainty_m"])
    )
    now = pd.Timestamp.now()

    def rank(p):
        d = pd.to_datetime(p["date"], errors="coerce")
        if pd.isna(d) or d < pd.Timestamp("1950-02-01") or d > now - pd.Timedelta(days=8):
            return (1, 0)
        return (0, -d.value)

    points = [points[i] for i in np.random.default_rng(42).permutation(len(points))]
    return sorted(points, key=rank)[: config["max_environment_points"]]


DEM_GLO90 = "https://copernicus-dem-90m.s3.amazonaws.com/"

_RASTERS = {}


def open_raster(path):
    """Cached rasterio dataset."""
    import rasterio

    if str(path) not in _RASTERS:
        _RASTERS[str(path)] = rasterio.open(path)
    return _RASTERS[str(path)]


def elevation_from_tiles(folder, points, download=True):
    """Copernicus DEM GLO-90 elevation from local 1-degree tiles."""
    rows = []
    tiles = {}
    for p in points:
        la, lo = math.floor(float(p["latitude"])), math.floor(float(p["longitude"]))
        tiles.setdefault(
            f"Copernicus_DSM_COG_30_{'N' if la >= 0 else 'S'}{abs(la):02d}_00_{'E' if lo >= 0 else 'W'}{abs(lo):03d}_00_DEM",
            [],
        ).append(p)
    for name, pts in tiles.items():
        path = Path(folder) / (name + ".tif")
        absent = Path(folder) / (name + ".absent")
        if absent.exists() or (not download and not path.exists()):
            continue
        try:
            fetch_file(DEM_GLO90 + name + "/" + name + ".tif", path)
        except FileNotFoundError:
            absent.touch()
            continue
        src = open_raster(path)
        for p, v in zip(pts, src.sample([(float(p["longitude"]), float(p["latitude"])) for p in pts])):
            if np.isfinite(v[0]) and v[0] != src.nodata:
                rows.append(
                    environment_row(
                        p,
                        "altitude",
                        float(v[0]),
                        "m",
                        "https://copernicus-dem-90m.s3.amazonaws.com/ | Copernicus DEM GLO-90",
                        resolution=90,
                        license="Copernicus DEM GLO-90, free licence",
                    )
                )
    return rows


def lapse_corrected_temperature(environment, rate=LAPSE_RATE_C_PER_M):
    """Correct NASA cell temperature to point elevation with a constant lapse rate."""
    if not len(environment) or "grid_elevation_m" not in environment:
        return []
    dem = (
        environment[environment.variable.eq("altitude")]
        .drop_duplicates("point_id")
        .set_index("point_id")
        .value
    )
    rows = []
    for _, r in environment[environment.variable.eq("temp_30d_cell")].iterrows():
        grid = pd.to_numeric(r.grid_elevation_m, errors="coerce")
        if pd.isna(grid) or r.point_id not in dem.index:
            continue
        point_dem = float(dem[r.point_id])
        rows.append(
            environment_row(
                dict(point_id=r.point_id),
                "temp_30d",
                float(r.value) - rate * (point_dem - float(grid)),
                "degC",
                "derived | NASA POWER T2M + lapse " + str(rate * 1000) + " degC/km + Copernicus DEM GLO-90",
                r.period_start,
                r.period_end,
                resolution=r.resolution_m,
                raw_cell_temp_30d=float(r.value),
                dem_elevation_m=point_dem,
                grid_elevation_m=float(grid),
            )
        )
    return rows


def elevation_check(occurrences, environment, tolerance=ELEVATION_QC_TOLERANCE_M):
    """Flag occurrences whose reported elevation differs from the DEM beyond the tolerance."""
    if "reported_elevation_m" not in occurrences or not len(environment):
        return pd.DataFrame()
    q = pd.DataFrame(
        {
            "occurrence_id": occurrences.occurrence_id,
            "reported_elevation_m": pd.to_numeric(occurrences.reported_elevation_m, errors="coerce"),
            "reported_accuracy_m": pd.to_numeric(
                occurrences.get("reported_elevation_accuracy_m", pd.Series(np.nan, index=occurrences.index)),
                errors="coerce",
            ),
        }
    ).dropna(subset=["reported_elevation_m"])
    dem = environment[
        environment.variable.eq("altitude")
        & pd.to_numeric(environment.resolution_m, errors="coerce").le(1000)
    ][["point_id", "value"]]
    q = q.merge(
        dem.rename(columns={"point_id": "occurrence_id", "value": "dem_elevation_m"}), on="occurrence_id"
    )
    q["difference_m"] = q.reported_elevation_m - q.dem_elevation_m
    q["suspect"] = q.difference_m.abs() > tolerance + q.reported_accuracy_m.fillna(0)
    return q


def try_source(client, label, fn, default):
    """Log a failed source (partial status) and continue without inventing values."""
    try:
        return fn()
    except NetworkDown:
        raise
    except Exception as exc:
        client.log(label, "FAILED", str(exc))
        print(label, "— unavailable:", str(exc)[:180], flush=True)
        return default


def persist_database(path, tables, chem, geo, environment_selected=None):
    """Write the release SQLite database with foreign keys."""
    path = Path(path)
    if path.exists():
        raise FileExistsError("Database already exists: choose a new run, " + str(path))
    conn = db_connect(path)
    conn.execute("PRAGMA foreign_keys=ON")
    keys = {
        "taxa": "taxon_id",
        "compounds": "compound_id",
        "samples": "sample_id",
        "occurrences": "occurrence_id",
        "chemistry": "record_id",
        "associations": "association_id",
        "environment": "env_id",
    }
    refs = {
        "samples": [("taxon_id", "taxa", "taxon_id")],
        "occurrences": [("taxon_id", "taxa", "taxon_id")],
        "chemistry": [("sample_id", "samples", "sample_id"), ("compound_id", "compounds", "compound_id")],
        "associations": [("taxon_id", "taxa", "taxon_id"), ("compound_id", "compounds", "compound_id")],
    }

    def sqlname(s):
        return '"' + s.replace('"', '""') + '"'

    def value(v):
        if isinstance(v, (dict, list, tuple)):
            return json.dumps(v, ensure_ascii=False)
        if pd.isna(v):
            return None
        if isinstance(v, (pd.Timestamp, datetime.datetime)):
            return v.isoformat()
        if isinstance(v, np.generic):
            return v.item()
        return v

    try:
        with conn:
            for name in [
                "taxa",
                "compounds",
                "samples",
                "occurrences",
                "chemistry",
                "associations",
                "environment",
            ]:
                df = (
                    environment_selected
                    if name == "environment" and environment_selected is not None
                    else tables[name]
                ).copy()
                cols = []
                for c in df:
                    kind = (
                        "INTEGER"
                        if pd.api.types.is_bool_dtype(df[c])
                        else "REAL"
                        if pd.api.types.is_numeric_dtype(df[c])
                        else "TEXT"
                    )
                    cols.append(sqlname(c) + " " + kind + (" PRIMARY KEY" if c == keys[name] else ""))
                cols += [
                    "FOREIGN KEY ("
                    + sqlname(col)
                    + ") REFERENCES "
                    + sqlname(parent)
                    + "("
                    + sqlname(pk)
                    + ")"
                    for col, parent, pk in refs.get(name, [])
                ]
                conn.execute("CREATE TABLE " + sqlname(name) + " (" + ",".join(cols) + ")")
                placeholders = ",".join("?" for _ in df.columns)
                conn.executemany(
                    "INSERT INTO " + sqlname(name) + " VALUES (" + placeholders + ")",
                    [[value(v) for v in row] for row in df.itertuples(index=False, name=None)],
                )
            chem.to_sql("chemical_features", conn, index=False, if_exists="fail")
            geo.to_sql("geographic_features", conn, index=False, if_exists="fail")
            conn.executescript("""CREATE INDEX idx_occ_taxon ON occurrences(taxon_id);
            CREATE INDEX idx_assoc_taxon ON associations(taxon_id);
            CREATE INDEX idx_assoc_compound ON associations(compound_id);
            CREATE VIEW regional_candidates AS
              SELECT o.occurrence_id,o.taxon_id,o.latitude,o.longitude,a.compound_id,
                     o.source AS botanical_source,a.source AS association_source,
                     'DOCUMENTED_ELSEWHERE_NOT_LOCAL_CHEMISTRY' AS evidence_status
              FROM occurrences o JOIN associations a ON o.taxon_id=a.taxon_id
              WHERE o.status='present';""")
        violations = conn.execute("PRAGMA foreign_key_check").fetchall()
        if violations:
            raise ValueError("SQL relation violation: " + str(violations))
        return {
            "path": str(path),
            "foreign_key_violations": 0,
            "counts": {n: conn.execute("SELECT COUNT(*) FROM " + sqlname(n)).fetchone()[0] for n in tables},
            "regional_candidates": conn.execute("SELECT COUNT(*) FROM regional_candidates").fetchone()[0],
        }
    finally:
        conn.close()


def features_from_database(path):
    """Read the chemical and geographic views back from the release database."""
    with db_connect(path) as conn:
        chem = pd.read_sql_query("SELECT * FROM chemical_features", conn)
        geo = pd.read_sql_query("SELECT * FROM geographic_features", conn)
    for df in [chem, geo]:
        for c in ["geo_valid", "quant_valid", "searched", "survey_complete"]:
            if c in df:
                df[c] = df[c].fillna(0).astype(bool)
        if "date" in df:
            df["date"] = pd.to_datetime(df.date, errors="coerce")
    return chem, geo


# ---------------- Global WFO queue

WFO_RECORD = "20782718"

WFO_FILE = "_DwC_backbone_R.zip"


def file_digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch_file(url, dest, expected=None, attempts=6):
    """Chunked download with HTTP Range resume."""
    dest = Path(dest)
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".part")
    for attempt in range(attempts):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": "AtlasPlantesResearch/3.1"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=60) as r:
                resumed = r.status == 206
                if resumed and not r.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                    raise ValueError("Inconsistent HTTP resume: " + url)
                if expected is None:
                    size = (
                        r.headers.get("Content-Range", "").rsplit("/", 1)[-1]
                        if resumed
                        else r.headers.get("Content-Length", "")
                    )
                    expected = int(size) if size.isdigit() else None
                total = offset if resumed else 0
                reported = total
                with partial.open("ab" if resumed else "wb") as f:
                    while block := r.read(1024 * 1024):
                        total += len(block)
                        if expected and total > expected:
                            raise ValueError("Downloaded size exceeds the manifest: " + url)
                        f.write(block)
                        if total - reported >= 20_000_000:
                            print(
                                f"  {dest.name} : {total / 1e6:.0f}"
                                + (f"/{expected / 1e6:.0f}" if expected else "")
                                + " MB",
                                flush=True,
                            )
                            reported = total
            if expected is None or partial.stat().st_size == expected:
                partial.replace(dest)
                return dest
            print(
                f"  {dest.name} : transfer interrupted at {partial.stat().st_size / 1e6:.0f} MB, resuming",
                flush=True,
            )
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise FileNotFoundError(url) from exc
            if exc.code == 416 and partial.exists():
                partial.replace(dest)
                return dest
            if attempt == attempts - 1:
                raise RuntimeError("Download interrupted, run again to resume: " + url) from exc
            time.sleep(2)
        except (OSError, urllib.error.URLError, TimeoutError, http.client.HTTPException) as exc:
            if attempt == attempts - 1:
                raise RuntimeError("Download interrupted, run again to resume: " + url) from exc
            time.sleep(2)
    raise RuntimeError("Incomplete download; run again to resume: " + url)


def download_catalog(root, offline=False, record=WFO_RECORD, filename=WFO_FILE):
    """Download a file from a frozen Zenodo record and verify its published checksum."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    metadata_path = root / ("zenodo_" + record + ".json")
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
    else:
        if offline:
            raise RuntimeError("Zenodo metadata missing in offline mode: " + record)
        with urllib.request.urlopen("https://zenodo.org/api/records/" + record, timeout=45) as r:
            metadata = json.load(r)
        metadata_path.write_text(json.dumps(metadata, ensure_ascii=False), encoding="utf-8")
    item = next((x for x in metadata["files"] if x["key"] == filename), None)
    if item is None:
        raise ValueError(filename + " not found in Zenodo record " + record)
    expected = item["size"]
    algorithm, checksum = item["checksum"].split(":", 1)
    if expected > 500_000_000:
        raise ValueError("Archive >500 MB: review the storage budget explicitly")
    dest = root / (("wfo_" + record + ".zip") if filename == WFO_FILE else filename)
    if dest.exists():
        if dest.stat().st_size == expected and file_digest(dest, algorithm) == checksum:
            return dest, metadata
        raise ValueError("Corrupted existing file: move " + str(dest) + " before downloading again")
    if offline:
        raise RuntimeError("File missing in offline mode: " + filename)
    print("Download:", filename, f"({expected / 1e6:.0f} MB)", flush=True)
    tmp = fetch_file(item["links"]["self"], dest.with_name(dest.name + ".download"), expected)
    if file_digest(tmp, algorithm) != checksum:
        raise ValueError("Checksum mismatch: delete " + str(tmp))
    tmp.replace(dest)
    return dest, metadata


def dwc_chunks(archive, chunksize=25000):
    """Stream the Darwin Core core file described by meta.xml."""
    z = zipfile.ZipFile(archive)
    try:
        metas = [n for n in z.namelist() if n.endswith("meta.xml")]
        if not metas and "classification.csv" in z.namelist():
            with z.open("classification.csv") as f:
                for chunk in pd.read_csv(
                    f,
                    sep="\t",
                    dtype=str,
                    keep_default_na=False,
                    chunksize=chunksize,
                    low_memory=False,
                    encoding="utf-8",
                    encoding_errors="replace",
                ):
                    yield chunk
            return
        if len(metas) != 1:
            raise ValueError("Expected one meta.xml, or classification.csv with a header")
        root = ET.fromstring(z.read(metas[0]))
        ns = {"d": "http://rs.tdwg.org/dwc/text/"}
        core = root.find("d:core", ns)
        if core is None:
            raise ValueError("Missing DwC core")
        location = core.find("d:files/d:location", ns).text
        member = str(Path(metas[0]).parent / location)
        if member not in z.namelist():
            member = location
        fieldmap = {
            int(f.attrib["index"]): f.attrib["term"].rsplit("/", 1)[-1]
            for f in core.findall("d:field", ns)
            if "index" in f.attrib
        }
        idel = core.find("d:id", ns)
        if idel is not None and int(idel.attrib["index"]) not in fieldmap:
            fieldmap[int(idel.attrib["index"])] = "taxonID"
        sep = core.attrib.get("fieldsTerminatedBy", "\t").replace("\\t", "\t")
        quote = core.attrib.get("fieldsEnclosedBy", '"') or '"'
        skip = int(core.attrib.get("ignoreHeaderLines", 0))
        with z.open(member) as f:
            reader = pd.read_csv(
                f,
                sep=sep,
                header=None,
                skiprows=skip,
                dtype=str,
                keep_default_na=False,
                encoding=core.attrib.get("encoding", "UTF-8"),
                quotechar=quote,
                chunksize=chunksize,
                low_memory=False,
            )
            for chunk in reader:
                yield chunk.rename(columns=fieldmap)
    finally:
        z.close()


def registry_status(db):
    """Inventory and queue counts of the global database."""
    with db_connect(db) as c:
        total = c.execute("SELECT COUNT(*) FROM catalog").fetchone()[0]
        statuses = dict(c.execute("SELECT status,COUNT(*) FROM jobs GROUP BY status"))
        ranks = pd.read_sql_query(
            "SELECT taxonomic_status,rank,COUNT(*) AS n FROM catalog GROUP BY taxonomic_status,rank", c
        )
        stages = pd.read_sql_query("SELECT status,COUNT(*) AS n FROM jobs GROUP BY status", c)
        covered = c.execute("SELECT COUNT(*) FROM jobs WHERE associations_count>0").fetchone()[0]
        geolocated = c.execute("SELECT COUNT(*) FROM jobs WHERE occurrences_count>0").fetchone()[0]
        entities = dict(c.execute("SELECT table_name,COUNT(*) FROM records GROUP BY table_name"))
    return (
        dict(
            catalog_names=total,
            accepted_species=sum(statuses.values()),
            queue=statuses,
            species_with_associations=covered,
            species_with_occurrences=geolocated,
            imported_records=entities,
        ),
        ranks,
        stages,
    )


def build_registry(archive, db, release=WFO_RECORD):
    """Load the WFO backbone and create one queue job per accepted species."""
    db = Path(db)
    db.parent.mkdir(parents=True, exist_ok=True)
    if db.exists() and db.stat().st_size == 0:
        db.unlink()
    fingerprint = file_digest(archive)
    if db.exists():
        with db_connect(db) as c:
            meta = dict(c.execute("SELECT key,value FROM meta"))
        if meta.get("catalog_sha256") != fingerprint:
            raise ValueError("Different WFO version: use another global directory")
        if meta.get("catalog_complete") == "true":
            return registry_status(db)[0]
        raise ValueError(
            "Incomplete inventory: rebuild into a new file without modifying the existing database"
        )
    temp = db.with_suffix(".building")
    temp.unlink(missing_ok=True)
    c = db_connect(temp)
    c.executescript("""CREATE TABLE meta(key TEXT PRIMARY KEY,value TEXT);
    CREATE TABLE catalog(wfo_id TEXT PRIMARY KEY,scientific_name TEXT,canonical_name TEXT,rank TEXT,
      taxonomic_status TEXT,accepted_id TEXT,parent_id TEXT,family TEXT,genus TEXT,authorship TEXT);
    CREATE TABLE jobs(wfo_id TEXT PRIMARY KEY,query_name TEXT,sort_key TEXT,status TEXT DEFAULT 'pending',
      attempts INTEGER DEFAULT 0,updated TEXT,error TEXT,gbif_taxon_id TEXT,occurrences_count INTEGER DEFAULT 0,
      associations_count INTEGER DEFAULT 0,environment_count INTEGER DEFAULT 0,
      FOREIGN KEY(wfo_id) REFERENCES catalog(wfo_id));
    CREATE TABLE records(table_name TEXT,record_id TEXT,payload TEXT,PRIMARY KEY(table_name,record_id));
    CREATE TABLE history(id INTEGER PRIMARY KEY AUTOINCREMENT,wfo_id TEXT,status TEXT,time TEXT,detail TEXT);
    CREATE TABLE crosswalk(wfo_id TEXT,gbif_taxon_id TEXT,PRIMARY KEY(wfo_id,gbif_taxon_id));""")
    try:
        c.execute("INSERT INTO meta VALUES (?,?)", ("catalog_sha256", fingerprint))
        c.execute("INSERT INTO meta VALUES (?,?)", ("release", release))
        total = 0
        encoding_issues = 0
        for chunk in dwc_chunks(archive):
            required = ["taxonID", "scientificName", "taxonRank", "taxonomicStatus"]
            if any(k not in chunk for k in required):
                raise ValueError("Unrecognized DwC schema: " + str(chunk.columns.tolist()))
            rows = []
            jobs = []
            for r in chunk.to_dict("records"):
                wid = str(r["taxonID"])
                name = str(r["scientificName"]).strip()
                rank = str(r["taxonRank"]).lower()
                status = str(r["taxonomicStatus"]).lower()
                if not wid or not name:
                    raise ValueError("Missing WFO name or identifier")
                genus = str(r.get("genus", ""))
                epithet = str(r.get("specificEpithet", ""))
                canonical = (
                    (genus + " " + epithet).strip() if genus and epithet and rank == "species" else name
                )
                authorship = str(r.get("scientificNameAuthorship", ""))
                rows.append(
                    (
                        wid,
                        name,
                        canonical,
                        rank,
                        status,
                        str(r.get("acceptedNameUsageID", "")),
                        str(r.get("parentNameUsageID", "")),
                        str(r.get("family", "")),
                        genus,
                        authorship,
                    )
                )
                if status == "accepted" and rank == "species":
                    query = name if not authorship or name.endswith(authorship) else name + " " + authorship
                    encoding_issue = any(
                        "\ufffd" in str(v)
                        for v in [wid, name, canonical, authorship, r.get("family", ""), genus]
                    )
                    encoding_issues += int(encoding_issue)
                    jobs.append(
                        (
                            wid,
                            query,
                            hashlib.sha256(wid.encode()).hexdigest(),
                            "review" if encoding_issue else "pending",
                            "Invalid encoding in an identity field; raw archive kept"
                            if encoding_issue
                            else "",
                        )
                    )
            with c:
                c.executemany("INSERT INTO catalog VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
                c.executemany(
                    "INSERT INTO jobs(wfo_id,query_name,sort_key,status,error) VALUES (?,?,?,?,?)", jobs
                )
            total += len(rows)
            if total % 100000 == 0:
                print("WFO names indexed:", total, flush=True)
        c.executescript(
            "CREATE INDEX idx_jobs_status ON jobs(status,sort_key); CREATE INDEX idx_catalog_names ON catalog(canonical_name);"
        )
        c.execute("INSERT INTO meta VALUES (?,?)", ("encoding_issues_accepted_species", str(encoding_issues)))
        c.execute("INSERT INTO meta VALUES (?,?)", ("catalog_complete", "true"))
        c.commit()
        c.close()
        temp.replace(db)
    except BaseException:
        c.close()
        raise
    return registry_status(db)[0]


def retry_global_jobs(db, include_partial=True, include_errors=True):
    """Requeue partial and failed species."""
    statuses = (["partial"] if include_partial else []) + (["error"] if include_errors else [])
    if not statuses:
        return 0
    with db_connect(db) as c:
        placeholders = ",".join("?" for _ in statuses)
        n = c.execute(
            "UPDATE jobs SET status=? WHERE status IN (" + placeholders + ")", ["pending", *statuses]
        ).rowcount
    return n


def recover_interrupted_jobs(db):
    """Requeue species left running by an interrupted launch."""
    with db_connect(db) as c:
        return c.execute("UPDATE jobs SET status='pending' WHERE status='running'").rowcount


def canonical_payload(row):
    def clean(v):
        if isinstance(v, dict):
            return {k: clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [clean(x) for x in v]
        if v is pd.NA or v is pd.NaT:
            return None
        if isinstance(v, np.generic):
            v = v.item()
        if isinstance(v, (pd.Timestamp, datetime.datetime)):
            return v.isoformat()
        if isinstance(v, float) and not np.isfinite(v):
            return None
        return v

    return json.dumps(clean(row), ensure_ascii=False, sort_keys=True)


def store_records(connection, frames, wfo_id):
    """Store canonical JSON records for one species."""
    keys = {
        "taxa": "taxon_id",
        "compounds": "compound_id",
        "samples": "sample_id",
        "chemistry": "record_id",
        "occurrences": "occurrence_id",
        "associations": "association_id",
        "environment": "env_id",
    }
    for table, df in frames.items():
        for row in df.to_dict("records"):
            payload = canonical_payload(row)
            rid = str(row[keys[table]])
            existing = connection.execute(
                "SELECT payload FROM records WHERE table_name=? AND record_id=?", (table, rid)
            ).fetchone()
            if existing and existing[0] != payload:
                connection.execute(
                    "INSERT INTO history(wfo_id,status,time,detail) VALUES (?,?,?,?)",
                    (
                        wfo_id,
                        "record_variant",
                        datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        canonical_payload({"table": table, "record_id": rid, "variant": row}),
                    ),
                )
            else:
                connection.execute("INSERT OR IGNORE INTO records VALUES (?,?,?)", (table, rid, payload))


def merge_batch_frames(batches, schemas):
    result = {}
    keys = {
        "taxa": "taxon_id",
        "compounds": "compound_id",
        "samples": "sample_id",
        "chemistry": "record_id",
        "occurrences": "occurrence_id",
        "associations": "association_id",
        "environment": "env_id",
    }
    for name, cols in schemas.items():
        parts = [b[name] for b in batches if len(b[name])]
        d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)
        if len(d):
            d = d.drop_duplicates(keys[name], keep="first")
        result[name] = d
    return result


def ensure_job_columns(c):
    """Add the priority and climate columns to an older queue."""
    cols = [r[1] for r in c.execute("PRAGMA table_info(jobs)")]
    if "priority" not in cols:
        c.execute("ALTER TABLE jobs ADD COLUMN priority INTEGER DEFAULT 1")
    if "climate" not in cols:
        c.execute("ALTER TABLE jobs ADD COLUMN climate TEXT")
        c.execute("UPDATE jobs SET climate='done' WHERE status!='pending'")
    c.execute("CREATE INDEX IF NOT EXISTS idx_jobs_queue ON jobs(status,priority,sort_key)")


def run_global_batch(db, online, global_config, collector, schemas):
    """Process one batch of queued species with a single worker."""
    cfg = dict(online)
    cfg.pop("taxon_names", None)
    profile = dict(cfg)
    for key in ["offline", "cache_dir", "bulk_db", "soil_dir", "dem_dir", "defer_climate"]:
        profile.pop(key, None)
    profile_hash = hashlib.sha256(json.dumps(profile, sort_keys=True, default=str).encode()).hexdigest()
    with db_connect(db) as c:
        ensure_job_columns(c)
        old = c.execute("SELECT value FROM meta WHERE key='collection_profile'").fetchone()
        if old and old[0] != profile_hash:
            raise ValueError("Collection profile changed: use another global directory for this profile")
        c.execute("INSERT OR IGNORE INTO meta VALUES (?,?)", ("collection_profile", profile_hash))
    start = time.monotonic()
    batches = []
    outcomes = []
    for _ in range(int(global_config["batch_size"])):
        if time.monotonic() - start > global_config["time_budget_seconds"]:
            break
        with db_connect(db) as c:
            c.execute("BEGIN IMMEDIATE")
            job = c.execute(
                "SELECT j.wfo_id,j.query_name,c.canonical_name FROM jobs j JOIN catalog c USING(wfo_id) WHERE j.status='pending' ORDER BY j.priority,j.sort_key LIMIT 1"
            ).fetchone()
            if not job:
                break
            wid, name, canonical = job
            c.execute(
                "UPDATE jobs SET status='running',attempts=attempts+1,updated=? WHERE wfo_id=?",
                (datetime.datetime.now(datetime.timezone.utc).isoformat(), wid),
            )
        print("Species:", name, "[", wid, "]", flush=True)
        try:
            one = dict(cfg, taxon_names=[name], wfo_id=wid)
            frames, report = collector(one, schemas)
            taxa = frames["taxa"]
            if taxa.empty:
                raise RuntimeError("No taxon resolved")
            mismatch = bool(canonical and not taxa.accepted_name.eq(canonical).all())
            if mismatch:
                state = "review"
                detail = (
                    "Accepted name differs between WFO and the source: "
                    + "; ".join(taxa.accepted_name)
                    + "; no automatic merge"
                )
                counts = {}
            else:
                errors = [e for e in report.get("events", []) if e.get("status") == "FAILED"]
                counts = {k: len(v) for k, v in frames.items()}
                state = (
                    "partial"
                    if errors
                    else ("no_data" if not counts["occurrences"] and not counts["associations"] else "done")
                )
                detail = json.dumps({"counts": counts, "errors": errors, "capped": True}, ensure_ascii=False)
            with db_connect(db) as c:
                if not mismatch:
                    store_records(c, frames, wid)
                    for tid in taxa.taxon_id:
                        c.execute("INSERT OR IGNORE INTO crosswalk VALUES (?,?)", (wid, tid))
                c.execute(
                    "UPDATE jobs SET status=?,updated=?,error=?,gbif_taxon_id=?,occurrences_count=?,associations_count=?,environment_count=? WHERE wfo_id=?",
                    (
                        state,
                        datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        detail,
                        ";".join(taxa.taxon_id),
                        counts.get("occurrences", 0),
                        counts.get("associations", 0),
                        counts.get("environment", 0),
                        wid,
                    ),
                )
                c.execute(
                    "INSERT INTO history(wfo_id,status,time,detail) VALUES (?,?,?,?)",
                    (wid, state, datetime.datetime.now(datetime.timezone.utc).isoformat(), detail),
                )
            if not mismatch:
                batches.append(frames)
            outcomes.append({"wfo_id": wid, "name": name, "status": state})
        except NetworkDown as exc:
            with db_connect(db) as c:
                c.execute("UPDATE jobs SET status='pending' WHERE wfo_id=?", (wid,))
            print(
                "Internet connection unavailable: batch stopped, species requeued:",
                str(exc)[:120],
                flush=True,
            )
            break
        except Exception as exc:
            state = "error"
            detail = str(exc)
            with db_connect(db) as c:
                c.execute(
                    "UPDATE jobs SET status=?,error=?,updated=? WHERE wfo_id=?",
                    (state, detail, datetime.datetime.now(datetime.timezone.utc).isoformat(), wid),
                )
                c.execute(
                    "INSERT INTO history(wfo_id,status,time,detail) VALUES (?,?,?,?)",
                    (wid, state, datetime.datetime.now(datetime.timezone.utc).isoformat(), detail),
                )
            outcomes.append({"wfo_id": wid, "name": name, "status": state, "error": detail})
            print("To retry:", detail[:180], flush=True)
    frames = merge_batch_frames(batches, schemas)
    status, ranks, stages = registry_status(db)
    return frames, dict(
        scope="all_accepted_species_in_WFO",
        registry=str(db),
        inventory=status,
        batch=outcomes,
        limitations=[
            "done means a capped collection pass, not complete scientific study",
            "WFO covers vascular plants and bryophytes; other botanical scopes need another inventory",
            "No_data means no result from configured queries, not biological absence",
            "Database grows across sessions; current notebook analyses the selected batch by default",
        ],
    )


# ---------------- Bulk sources: LOTUS, GBIF snapshot, SoilGrids, weather

LOTUS_RECORD = "19360665"

LOTUS_FILE = "260413_frozen_metadata.csv.gz"

LOTUS_LICENSE = "CC-BY-4.0 ; LOTUS frozen 2026-04-13, doi:10.5281/zenodo.19360665"

SOILGRIDS_1KM = "https://files.isric.org/soilgrids/latest/data_aggregated/1000m/"

GBIF_BUCKET = "https://gbif-open-data-us-east-1.s3.amazonaws.com/"

GEOSPATIAL_ISSUES = [
    "ZERO_COORDINATE",
    "COORDINATE_OUT_OF_RANGE",
    "COORDINATE_INVALID",
    "COUNTRY_COORDINATE_MISMATCH",
]


def wfo_name_index(c):
    """Canonical name -> accepted WFO species; homonyms are excluded."""
    cat = pd.read_sql_query(
        """SELECT k.canonical_name,k.taxonomic_status='accepted' AS accepted,
                             CASE WHEN k.taxonomic_status='accepted' THEN k.wfo_id ELSE k.accepted_id END AS target
                             FROM catalog k WHERE k.rank='species'""",
        c,
    )
    cat = cat[cat.target.isin(set(r[0] for r in c.execute("SELECT wfo_id FROM jobs")))]
    accepted = cat[cat.accepted.eq(1)]
    pool = pd.concat(
        [accepted, cat[cat.accepted.eq(0) & ~cat.canonical_name.isin(accepted.canonical_name)]]
    ).drop_duplicates(["canonical_name", "target"])
    n = pool.groupby("canonical_name").target.transform("nunique")
    return pool[n.eq(1)].set_index("canonical_name").target, int(pool[n.gt(1)].canonical_name.nunique())


def wfo_taxa(c, ids):
    """Accepted WFO taxa in the release schema."""
    rows = []
    for i in range(0, len(ids), 900):
        part = list(ids[i : i + 900])
        rows += c.execute(
            f"SELECT wfo_id,scientific_name,authorship,canonical_name,family,genus FROM catalog WHERE wfo_id IN ({','.join('?' * len(part))})",
            part,
        ).fetchall()
    return pd.DataFrame(
        [
            dict(
                taxon_id="WFO:" + w,
                published_name=(n + " " + a).strip() if a and not n.endswith(a) else n,
                accepted_name=cn,
                family=f or "unknown",
                genus=g or "unknown",
                source="https://www.worldfloraonline.org/taxon/" + w,
                wfo_id=w,
            )
            for w, n, a, cn, f, g in rows
        ],
        columns=SCHEMAS["taxa"] + ["wfo_id"],
    )


def build_lotus(db, archive):
    """Map LOTUS organism names to WFO species and store associations and compounds."""
    with db_connect(db) as c:
        done = c.execute("SELECT value FROM meta WHERE key='lotus_summary'").fetchone()
        if done:
            return json.loads(done[0])
        index, ambiguous = wfo_name_index(c)
    cols = [
        "structure_wikidata",
        "structure_inchikey",
        "structure_smiles",
        "structure_cid",
        "structure_nameTraditional",
        "structure_nameIupac",
        "structure_molecular_formula",
        "structure_taxonomy_npclassifier_01pathway",
        "structure_taxonomy_npclassifier_03class",
        "organism_name",
        "organism_taxonomy_09species",
        "reference_doi",
        "reference_wikidata",
    ]
    lotus = pd.read_csv(archive, usecols=cols, dtype=str, keep_default_na=False)
    lotus["wfo_id"] = lotus.organism_name.map(index).fillna(lotus.organism_taxonomy_09species.map(index))
    m = lotus[lotus.wfo_id.notna() & lotus.structure_inchikey.ne("")].copy()
    qid = m.structure_wikidata.str.rsplit("/", n=1).str[-1]
    m["compound_id"] = np.where(qid.str.startswith("Q"), "WD:" + qid, "IK:" + m.structure_inchikey)
    m["proof"] = np.where(m.reference_doi.ne(""), "https://doi.org/" + m.reference_doi, m.reference_wikidata)
    m = m[m.proof.ne("")]
    compounds = [
        dict(
            compound_id=r.compound_id,
            compound_name=r.structure_nameTraditional or r.structure_nameIupac or r.structure_inchikey,
            smiles=r.structure_smiles,
            chemical_class=r.structure_taxonomy_npclassifier_03class or "unknown",
            source=r.structure_wikidata or "https://pubchem.ncbi.nlm.nih.gov/compound/" + r.structure_cid,
            inchikey=r.structure_inchikey,
            pubchem_cid=r.structure_cid,
            molecular_formula=r.structure_molecular_formula,
            pathway=r.structure_taxonomy_npclassifier_01pathway,
            license=LOTUS_LICENSE,
        )
        for r in m.drop_duplicates("compound_id").itertuples()
    ]
    associations = {}
    for r in m.itertuples():
        aid = "LOTUS:" + digest(r.wfo_id + "|" + r.compound_id + "|" + r.proof)[:24]
        associations[aid] = dict(
            association_id=aid,
            taxon_id="WFO:" + r.wfo_id,
            compound_id=r.compound_id,
            source=r.proof,
            organism_name_source=r.organism_name,
            evidence="referenced_association",
            local_evidence=False,
            license=LOTUS_LICENSE,
        )
    counts = m.groupby("wfo_id").compound_id.nunique()
    with db_connect(db) as c:
        taxa = wfo_taxa(c, counts.index.tolist())
        rows = [("taxa", t["taxon_id"], canonical_payload(t)) for t in taxa.to_dict("records")]
        rows += [("compounds", x["compound_id"], canonical_payload(x)) for x in compounds]
        rows += [("associations", k, canonical_payload(v)) for k, v in associations.items()]
        c.executemany("INSERT OR IGNORE INTO records VALUES (?,?,?)", rows)
        c.execute(
            "CREATE INDEX IF NOT EXISTS idx_records_taxon ON records(table_name,json_extract(payload,'$.taxon_id'))"
        )
        if "priority" not in [r[1] for r in c.execute("PRAGMA table_info(jobs)")]:
            c.execute("ALTER TABLE jobs ADD COLUMN priority INTEGER DEFAULT 1")
        c.executemany(
            "UPDATE jobs SET priority=0,associations_count=? WHERE wfo_id=?",
            [(int(n), w) for w, n in counts.items()],
        )
        summary = dict(
            lotus_pairs=len(lotus),
            pairs_matched_to_wfo=int(len(m)),
            species_with_chemistry=int(len(counts)),
            compounds=len(compounds),
            associations=len(associations),
            ambiguous_names_excluded=ambiguous,
            source=LOTUS_LICENSE,
        )
        c.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", ("lotus_summary", json.dumps(summary)))
    return summary


def soil_rasters(folder, offline=False):
    """Download the SoilGrids 1 km maps (4 properties x 3 depths) once."""
    for prop in SOIL_PROPS:
        for depth in ["0-5cm", "5-15cm", "15-30cm"]:
            name = f"{prop}_{depth}_mean_1000.tif"
            if not (Path(folder) / name).exists():
                if offline:
                    raise RuntimeError("SoilGrids map missing in offline mode: " + name)
                print("SoilGrids 1 km:", name, flush=True)
                fetch_file(SOILGRIDS_1KM + prop + "/" + name, Path(folder) / name)
    return Path(folder)


def soil_from_rasters(folder, points, props=None):
    """Thickness-weighted 0-30 cm soil values from the local 1 km maps."""
    from rasterio.warp import transform

    rows = []
    if not points:
        return rows
    lons = [float(p["longitude"]) for p in points]
    lats = [float(p["latitude"]) for p in points]
    for prop in props or list(SOIL_PROPS):
        variable, divisor, unit = SOIL_PROPS[prop]
        layers = []
        for depth, weight in [("0-5cm", 5), ("5-15cm", 10), ("15-30cm", 15)]:
            src = open_raster(Path(folder) / f"{prop}_{depth}_mean_1000.tif")
            xs, ys = transform("EPSG:4326", src.crs, lons, lats)
            layers.append(([float(v[0]) for v in src.sample(zip(xs, ys))], weight))
        for i, p in enumerate(points):
            raw = [vals[i] for vals, _ in layers]
            if any(v < 0 or not np.isfinite(v) for v in raw):
                continue
            parts = [(v / divisor, w) for v, (_, w) in zip(raw, layers)]
            rows.append(
                environment_row(
                    p,
                    variable,
                    sum(v * w for v, w in parts) / 30,
                    unit,
                    "https://files.isric.org/soilgrids/ | SoilGrids " + prop + " mean 1 km",
                    depth_top=0,
                    depth_bottom=30,
                    resolution=1000,
                    aggregation="thickness-weighted profile descriptor; NOT bulk measurement",
                    depth_values=json.dumps(parts),
                    license="CC-BY-4.0",
                )
            )
    return rows


def gbif_tables(c):
    c.executescript("""CREATE TABLE IF NOT EXISTS gbif_files(key TEXT PRIMARY KEY,done INTEGER DEFAULT 0,rows INTEGER);
    CREATE TABLE IF NOT EXISTS gbif_sample(wfo_id TEXT,gbifid TEXT,imprecise INTEGER,h TEXT,payload TEXT,PRIMARY KEY(wfo_id,gbifid));""")


def gbif_snapshot(c, offline=False):
    """Freeze the public GBIF snapshot (Parquet on AWS) and its file list."""
    snap = c.execute("SELECT value FROM meta WHERE key='gbif_snapshot'").fetchone()
    if snap:
        return snap[0]
    if offline:
        raise RuntimeError("GBIF snapshot not frozen: a first online launch is required")
    listing = lambda q: urllib.request.urlopen(GBIF_BUCKET + "?list-type=2&" + q, timeout=60).read().decode()
    snapshot = sorted(
        re.findall(
            r"<Prefix>(occurrence/\d{4}-\d{2}-\d{2}/)</Prefix>", listing("prefix=occurrence/&delimiter=/")
        )
    )[-1]
    keys = []
    token = None
    while True:
        x = listing(
            "prefix="
            + urllib.parse.quote(snapshot + "occurrence.parquet/")
            + "&max-keys=1000"
            + ("&continuation-token=" + urllib.parse.quote(token) if token else "")
        )
        keys += [k for k in re.findall(r"<Key>([^<]+)</Key>", x) if not k.endswith("/")]
        token = (re.findall(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", x) or [None])[0]
        if not token:
            break
    c.executemany("INSERT OR IGNORE INTO gbif_files(key) VALUES (?)", [(k,) for k in keys])
    c.execute("INSERT INTO meta VALUES (?,?)", ("gbif_snapshot", snapshot))
    return snapshot


def gbif_status(db):
    with db_connect(db) as c:
        gbif_tables(c)
        snap = c.execute("SELECT value FROM meta WHERE key='gbif_snapshot'").fetchone()
        total, done = c.execute("SELECT COUNT(*),COALESCE(SUM(done),0) FROM gbif_files").fetchone()
        species, rows = c.execute("SELECT COUNT(DISTINCT wfo_id),COUNT(*) FROM gbif_sample").fetchone()
    return dict(
        snapshot=snap[0] if snap else None,
        files_done=int(done),
        files_total=int(total),
        complete=bool(total) and done == total,
        species_sampled=int(species),
        occurrences_sampled=int(rows),
    )


def scan_gbif(
    db,
    max_files,
    k=30,
    max_uncertainty_m=MAX_COARSE_UNCERTAINTY_M,
    offline=False,
    chunk=16,
    uri=lambda key: "s3://gbif-open-data-us-east-1/" + key,
    recent=5,
):
    """Resumable DuckDB scan of the GBIF snapshot into a hashed per-species sample."""
    import duckdb

    with db_connect(db) as c:
        gbif_tables(c)
        snapshot = gbif_snapshot(c, offline)
        todo = [
            r[0]
            for r in c.execute(
                "SELECT key FROM gbif_files WHERE done=0 ORDER BY key LIMIT ?", (int(max_files),)
            )
        ]
        index, _ = wfo_name_index(c) if todo else (None, 0)
    if not todo:
        return gbif_status(db)
    found = re.search(r"\d{4}-\d{2}-\d{2}", snapshot)
    cutoff = found[0] if found else datetime.date.today().isoformat()
    con = duckdb.connect()
    con.execute("SET memory_limit='2GB'")
    with contextlib.suppress(duckdb.Error):
        con.execute("SET enable_external_file_cache=false")
    if any(uri(x).startswith("s3://") for x in todo[:1]):
        con.execute(
            "INSTALL httpfs; LOAD httpfs; SET s3_region='us-east-1'; SET http_timeout=45; SET http_retries=3; SET http_retry_wait_ms=2000; SET http_retry_backoff=2; SET http_keep_alive=false;"
        )
    con.register("names", pd.DataFrame({"species": index.index, "wfo_id": index.values}))
    issues = "[" + ",".join("'" + i + "'" for i in GEOSPATIAL_ISSUES) + "]"
    lo = pd.Timestamp("1950-02-01")

    def day_number(dates):
        d = pd.to_datetime(dates, errors="coerce")
        good = (d >= lo) & (d <= pd.Timestamp(cutoff))
        return ((d - pd.Timestamp("1970-01-01")).dt.days.where(good, -1)).fillna(-1).astype(int)

    for i in range(0, len(todo), chunk):
        part = todo[i : i + chunk]
        files = "[" + ",".join("'" + uri(x).replace("'", "''") + "'" for x in part) + "]"
        query = f"""SELECT n.wfo_id,o.gbifid,(o.coordinateuncertaintyinmeters>1000)::INTEGER AS imprecise,md5(o.gbifid) AS h,
               o.datasetkey,o.decimallatitude,o.decimallongitude,o.coordinateuncertaintyinmeters,
               CASE WHEN o.day IS NOT NULL AND o.month IS NOT NULL THEN strftime(o.eventdate,'%Y-%m-%d') END AS event_date,
               o.basisofrecord,o.license,o.establishmentmeans,o.elevation,o.elevationaccuracy,o.issue
             FROM read_parquet({files}) o JOIN names n ON n.species=o.species
             WHERE o.kingdom='Plantae' AND o.occurrencestatus='PRESENT' AND o.decimallatitude IS NOT NULL AND o.decimallongitude IS NOT NULL
               AND o.coordinateuncertaintyinmeters BETWEEN 0 AND {float(max_uncertainty_m)}
               AND NOT list_has_any(coalesce(o.issue::VARCHAR[],[]::VARCHAR[]),{issues}::VARCHAR[])
             QUALIFY row_number() OVER (PARTITION BY n.wfo_id ORDER BY imprecise,h)<={int(k)}
                  OR row_number() OVER (PARTITION BY n.wfo_id ORDER BY imprecise,
                       (event_date IS NULL OR event_date<'1950-02-01' OR event_date>'{cutoff}'),
                       CASE WHEN event_date>='1950-02-01' AND event_date<='{cutoff}' THEN event_date END DESC NULLS LAST,h)<={int(recent)}"""
        for attempt in range(3):
            try:
                df = con.execute(query).df()
                break
            except (duckdb.IOException, duckdb.HTTPException) as exc:
                print(f"GBIF: network interruption ({str(exc)[:90]}…), retrying in 60 s", flush=True)
                time.sleep(60)
        else:
            print("GBIF: network unavailable, scan paused; run again to resume", flush=True)
            return gbif_status(db)
        df = df.astype(object).where(df.notna(), None)
        with db_connect(db) as c:
            c.execute("PRAGMA journal_mode=WAL")
            c.execute("PRAGMA synchronous=NORMAL")
            c.execute("CREATE TEMP TABLE touched(wfo_id TEXT PRIMARY KEY)")
            c.executemany("INSERT OR IGNORE INTO touched VALUES (?)", [(w,) for w in set(df.wfo_id)])
            old = pd.DataFrame(
                c.execute(
                    "SELECT s.wfo_id,s.gbifid,s.imprecise,s.h,json_extract(s.payload,'$.eventDate') FROM gbif_sample s JOIN touched USING(wfo_id)"
                ).fetchall(),
                columns=["wfo_id", "gbifid", "imprecise", "h", "event_date"],
            ).assign(new=False)
            cand = df[["wfo_id", "gbifid", "imprecise", "h", "event_date"]].assign(new=True)
            merged = pd.concat([old, cand]).astype({"imprecise": int}).drop_duplicates(["wfo_id", "gbifid"])
            merged["day"] = day_number(merged.event_date)
            by_hash = merged.sort_values(["wfo_id", "imprecise", "h"]).groupby("wfo_id").head(int(k))
            by_date = (
                merged.sort_values(["wfo_id", "imprecise", "day", "h"], ascending=[True, True, False, True])
                .groupby("wfo_id")
                .head(int(recent))
            )
            keep = pd.concat([by_hash, by_date]).drop_duplicates(["wfo_id", "gbifid"])
            kept = set(zip(keep.wfo_id, keep.gbifid))
            winners = set(zip(keep[keep.new].wfo_id, keep[keep.new].gbifid))
            c.executemany(
                "DELETE FROM gbif_sample WHERE wfo_id=? AND gbifid=?",
                [x for x in zip(old.wfo_id, old.gbifid) if x not in kept],
            )
            c.executemany(
                "INSERT OR IGNORE INTO gbif_sample VALUES (?,?,?,?,?)",
                [
                    (
                        r.wfo_id,
                        r.gbifid,
                        int(r.imprecise),
                        r.h,
                        canonical_payload(
                            dict(
                                key=r.gbifid,
                                datasetKey=r.datasetkey,
                                decimalLatitude=r.decimallatitude,
                                decimalLongitude=r.decimallongitude,
                                coordinateUncertaintyInMeters=r.coordinateuncertaintyinmeters,
                                eventDate=r.event_date or "",
                                basisOfRecord=r.basisofrecord,
                                license=r.license,
                                establishmentMeans=r.establishmentmeans,
                                elevation=r.elevation,
                                elevationAccuracy=r.elevationaccuracy,
                                issues=list(r.issue) if r.issue is not None else [],
                            )
                        ),
                    )
                    for r in df.itertuples()
                    if (r.wfo_id, r.gbifid) in winners
                ],
            )
            c.executemany(
                "UPDATE gbif_files SET done=1,rows=? WHERE key=?",
                [(len(df) if j == 0 else 0, x) for j, x in enumerate(part)],
            )
        st = gbif_status(db)
        print(
            f"GBIF: {st['files_done']}/{st['files_total']} files read; {st['species_sampled']} species, {st['occurrences_sampled']} occurrences retained",
            flush=True,
        )
    return gbif_status(db)


def climate_rows(client, point, config):
    """30-day weather (Open-Meteo, else NASA POWER with lapse-rate correction)."""
    climate = try_source(client, "Weather", lambda: fetch_weather(client, point), [])
    if (
        config.get("dem_dir")
        and point["geo_uncertainty_m"] <= config["max_geo_uncertainty_m"]
        and any(r["variable"] == "temp_30d_cell" for r in climate)
    ):
        altitude = try_source(
            client, "Elevation", lambda: elevation_from_tiles(config["dem_dir"], [point], download=False), []
        )
        climate += altitude + lapse_corrected_temperature(pd.DataFrame(climate + altitude))
    return climate


def climate_pass(db, config, budget_s, wait_s=900):
    """Second pass: weather for species already collected, within service quotas."""
    client = public_client(config)
    client.nasa_fallback = False
    start = time.monotonic()
    done = 0
    failures = 0
    with db_connect(db) as c:
        ensure_job_columns(c)
    while time.monotonic() - start < budget_s:
        with db_connect(db) as c:
            todo = [
                r[0]
                for r in c.execute(
                    "SELECT wfo_id FROM jobs WHERE status IN ('done','partial','no_data') AND climate IS NULL ORDER BY priority,sort_key LIMIT 200"
                )
            ]
        if not todo:
            break
        for wid in todo:
            if time.monotonic() - start >= budget_s:
                break
            with db_connect(db) as c:
                occ = [
                    o
                    for o in (
                        parse_gbif_occurrence(json.loads(p), dict(taxon_id="WFO:" + wid))
                        for (p,) in c.execute("SELECT payload FROM gbif_sample WHERE wfo_id=?", (wid,))
                    )
                    if o
                ]
            client.events = []
            rows = []
            client.om_failed = False
            try:
                for point in climate_points(occ, config):
                    rows += climate_rows(client, point, config)
            except NetworkDown as exc:
                print("Internet connection unavailable: weather pass paused:", str(exc)[:100], flush=True)
                return done
            if client.om_failed:
                failures += 1
                if "open-meteo-archive" in client.disabled or failures >= 3:
                    print(
                        f"Open-Meteo unavailable (quota or outage): waiting {wait_s / 60:.0f} min",
                        flush=True,
                    )
                    time.sleep(min(wait_s, max(0, budget_s - (time.monotonic() - start))))
                    client.disabled.discard("open-meteo-archive")
                    client.om_failures = 0
                    failures = 0
                break
            failures = 0
            state = (
                "failed"
                if any(e["status"] == "FAILED" for e in client.events) and not rows
                else ("done" if rows else "none")
            )
            with db_connect(db) as c:
                if rows:
                    store_records(c, {"environment": pd.DataFrame(rows)}, wid)
                c.execute(
                    "UPDATE jobs SET climate=?,environment_count=environment_count+? WHERE wfo_id=?",
                    (state, len(rows), wid),
                )
            done += 1
            if done % 200 == 0:
                print(f"Weather: {done} species processed in this pass", flush=True)
    return done


def collect_bulk_sources(config, schemas=SCHEMAS):
    """Collect one species from the local copies; only weather queries a service."""
    client = public_client(config)
    wid = config["wfo_id"]
    tid = "WFO:" + wid
    with db_connect(config["bulk_db"]) as c:
        taxa = wfo_taxa(c, [wid])
        assoc = [
            json.loads(p)
            for (p,) in c.execute(
                "SELECT payload FROM records WHERE table_name='associations' AND json_extract(payload,'$.taxon_id')=?",
                (tid,),
            )
        ]
        cids = sorted({a["compound_id"] for a in assoc})
        mols = (
            [
                json.loads(p)
                for (p,) in c.execute(
                    f"SELECT payload FROM records WHERE table_name='compounds' AND record_id IN ({','.join('?' * len(cids))})",
                    cids,
                )
            ]
            if cids
            else []
        )
        occ = [
            o
            for o in (
                parse_gbif_occurrence(json.loads(p), dict(taxon_id=tid))
                for (p,) in c.execute(
                    "SELECT payload FROM gbif_sample WHERE wfo_id=? ORDER BY imprecise,h", (wid,)
                )
            )
            if o
        ]
    tables = {k: [] for k in schemas}
    tables["taxa"] = taxa.to_dict("records")
    tables["associations"] = assoc
    tables["compounds"] = mols
    tables["occurrences"] = occ
    precise = admissible_points(occ, config["max_geo_uncertainty_m"])
    folder = Path(config["soil_dir"]) if config.get("soil_dir") else None
    local = [
        p
        for p in SOIL_PROPS
        if folder
        and all((folder / f"{p}_{d}_mean_1000.tif").exists() for d in ["0-5cm", "5-15cm", "15-30cm"])
    ]
    if precise and local:
        tables["environment"] += try_source(
            client, "SoilGrids 1 km", lambda: soil_from_rasters(folder, precise, local), []
        )
    for point in [] if config.get("defer_climate") else climate_points(occ, config):
        tables["environment"] += climate_rows(client, point, config)
        if point["geo_uncertainty_m"] <= config["max_geo_uncertainty_m"]:
            for prop in [p for p in SOIL_PROPS if p not in local]:
                tables["environment"] += try_source(
                    client, "SoilGrids WCS " + prop, lambda pr=prop: fetch_soil(client, point, [pr]), []
                )
    frames = {
        name: pd.DataFrame(tables[name]).reindex(
            columns=list(dict.fromkeys(cols + [k for r in tables[name] for k in r]))
        )
        for name, cols in schemas.items()
    }
    return frames, {
        "config": {k: v for k, v in config.items() if k != "taxon_names"},
        "counts": {k: len(v) for k, v in frames.items()},
        "events": client.events,
    }


def load_selection(db, priority_only=True, schemas=SCHEMAS):
    """Load the collected records of the selected species for the release."""
    with db_connect(db) as c:
        ids = [
            r[0]
            for r in c.execute(
                "SELECT wfo_id FROM jobs WHERE status IN ('done','partial')"
                + (" AND priority=0" if priority_only else "")
            )
        ]
        c.execute("CREATE TEMP TABLE sel(taxon_id TEXT PRIMARY KEY)")
        c.executemany("INSERT INTO sel VALUES (?)", [("WFO:" + w,) for w in ids])
        part = lambda table: [
            json.loads(p)
            for (p,) in c.execute(
                "SELECT r.payload FROM sel JOIN records r ON r.table_name=? AND json_extract(r.payload,'$.taxon_id')=sel.taxon_id",
                (table,),
            )
        ]
        tables = {name: part(name) for name in ["taxa", "associations", "occurrences"]}
        c.execute("CREATE TEMP TABLE ids(id TEXT PRIMARY KEY)")
        c.executemany(
            "INSERT OR IGNORE INTO ids VALUES (?)", [(a["compound_id"],) for a in tables["associations"]]
        )
        tables["compounds"] = [
            json.loads(p)
            for (p,) in c.execute(
                "SELECT r.payload FROM ids JOIN records r ON r.table_name='compounds' AND r.record_id=ids.id"
            )
        ]
        c.execute("DELETE FROM ids")
        c.executemany(
            "INSERT OR IGNORE INTO ids VALUES (?)", [(o["occurrence_id"],) for o in tables["occurrences"]]
        )
        tables["environment"] = [
            json.loads(p)
            for (p,) in c.execute(
                "SELECT payload FROM records WHERE table_name='environment' AND json_extract(payload,'$.point_id') IN (SELECT id FROM ids)"
            )
        ]
    print("Atlas selection:", {k: len(v) for k, v in tables.items()}, flush=True)
    return {
        name: pd.DataFrame(tables.get(name, [])).reindex(
            columns=list(dict.fromkeys(cols + [k for r in tables.get(name, [])[:1] for k in r]))
        )
        for name, cols in schemas.items()
    }


def atlas_summaries(clean, geo, out):
    """Species, family, compound and class summaries, and the summary page."""
    out = Path(out)
    taxa = clean["taxa"]
    assoc = clean["associations"]
    mol = clean["compounds"]
    occ = clean["occurrences"]
    env = [v for v in NUM_ENV if v in geo]
    precise = geo[geo.geo_valid]
    sp = taxa[["taxon_id", "accepted_name", "family", "genus"]].copy()
    sp["molecules"] = sp.taxon_id.map(assoc.groupby("taxon_id").compound_id.nunique()).fillna(0).astype(int)
    sp["references"] = sp.taxon_id.map(assoc.groupby("taxon_id").source.nunique()).fillna(0).astype(int)
    sp["observations"] = sp.taxon_id.map(occ.groupby("taxon_id").size()).fillna(0).astype(int)
    sp["precise_observations"] = sp.taxon_id.map(precise.groupby("taxon_id").size()).fillna(0).astype(int)
    sp = sp.merge(
        precise.groupby("taxon_id")[env].median().round(2).add_prefix("median_"),
        left_on="taxon_id",
        right_index=True,
        how="left",
    )
    sp = sp.sort_values(["molecules", "observations"], ascending=False)
    sp.to_csv(out / "atlas_species.csv", index=False)
    map_points = write_atlas_map(sp, geo, out / "atlas.html")
    fam_mol = assoc.merge(taxa[["taxon_id", "family"]], on="taxon_id").groupby("family").compound_id.nunique()
    fam = (
        sp.groupby("family")
        .agg(
            species=("taxon_id", "size"),
            observations=("observations", "sum"),
            **{"median_" + v: ("median_" + v, "median") for v in env},
        )
        .round(2)
    )
    fam.insert(1, "distinct_molecules", fam.index.map(fam_mol).fillna(0).astype(int))
    fam = fam.sort_values("distinct_molecules", ascending=False)
    fam.to_csv(out / "atlas_families.csv")
    links = assoc.merge(taxa[["taxon_id", "family"]], on="taxon_id")
    mo = mol[
        [
            c
            for c in [
                "compound_id",
                "compound_name",
                "chemical_class",
                "pathway",
                "inchikey",
                "pubchem_cid",
                "smiles",
            ]
            if c in mol
        ]
    ].copy()
    mo["species"] = mo.compound_id.map(links.groupby("compound_id").taxon_id.nunique()).fillna(0).astype(int)
    mo["families"] = mo.compound_id.map(links.groupby("compound_id").family.nunique()).fillna(0).astype(int)
    mo = mo.sort_values("species", ascending=False)
    mo.to_csv(out / "atlas_molecules.csv", index=False)
    cl = (
        mo.groupby("chemical_class")
        .agg(molecules=("compound_id", "size"), max_species=("species", "max"))
        .sort_values("molecules", ascending=False)
    )
    cl.to_csv(out / "atlas_chemical_classes.csv")
    stats = {
        "plants": len(sp),
        "families": int(sp.family.nunique()),
        "molecules": len(mo),
        "plant_molecule_links": len(assoc),
        "observations": len(occ),
        "precise_observations": len(precise),
        "plants_with_weather": int(sp["median_temp_30d"].notna().sum()) if "median_temp_30d" in sp else 0,
        "plants_with_soil": int(sp["median_soil_ph"].notna().sum()) if "median_soil_ph" in sp else 0,
    }
    stats["map_points"] = map_points
    json_write(out / "atlas_coverage.json", stats)
    table = lambda df, n=20: df.head(n).to_html(index=False, border=0, classes="t", na_rep="—")
    labels = {
        "plants": "Species",
        "families": "Families",
        "molecules": "Molecules",
        "plant_molecule_links": "Species–molecule links",
        "observations": "Observations",
        "precise_observations": "Precise occurrences (≤ 1 km)",
        "plants_with_weather": "Species with weather",
        "plants_with_soil": "Species with soil",
        "map_points": "Points on the map",
    }
    cards = "".join(
        '<div class="c"><b>' + f"{v:,}" + f"</b><span>{html.escape(labels[k])}</span></div>"
        for k, v in stats.items()
    )
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Plant atlas</title><style>:root{{--bg:#f6f8f5;--fg:#18241d;--mut:#5b6b61;--card:#fff;--line:#dfe6e0;--acc:#2f6f4e}}
@media (prefers-color-scheme:dark){{:root{{--bg:#111a15;--fg:#e8efe9;--mut:#a3b3a8;--card:#18241d;--line:#2a3a31;--acc:#7cc39b}}}}
body{{margin:0;padding:24px 16px;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}}main{{max-width:1100px;margin:auto}}
h1{{margin:0 0 4px}}p{{color:var(--mut)}}.g{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin:18px 0}}
.c{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px}}.c b{{display:block;font-size:22px;color:var(--acc)}}
.c span{{color:var(--mut);font-size:13px}}h2{{margin-top:28px}}.w{{overflow-x:auto}}table.t{{border-collapse:collapse;width:100%;background:var(--card)}}
.t th,.t td{{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}.t th{{color:var(--mut);font-weight:600}}a{{color:var(--acc)}}</style>
<main><h1>Plant atlas — soil, weather, molecules</h1>
<p>Public data (WFO, LOTUS, GBIF, SoilGrids, Open-Meteo / NASA POWER, Copernicus). Literature associations are not evidence
of the local presence of a molecule. Map: <a href="atlas.html">atlas.html</a> · full tables: atlas_*.csv files.</p>
<div class="g">{cards}</div>
<h2>Species with the most known molecules</h2><div class="w">{table(sp[["accepted_name", "family", "molecules", "observations"] + [c for c in ["median_altitude", "median_temp_30d", "median_soil_ph"] if c in sp]])}</div>
<h2>Families</h2><div class="w">{table(fam.reset_index()[["family", "species", "distinct_molecules", "observations"]])}</div>
<h2>Most represented chemical classes</h2><div class="w">{table(cl.reset_index())}</div>
<h2>Molecules found in the most species</h2><div class="w">{table(mo[["compound_name", "chemical_class", "species", "families"]])}</div></main></html>"""
    (out / "atlas_summary.html").write_text(page, encoding="utf-8")
    return stats


def show(title, obj):
    print("\n== " + title)
    print(obj.to_string() if isinstance(obj, pd.DataFrame) and len(obj) else obj)


# ---------------------------------------------------------------- Pipeline


def run(args):
    """Collect one batch into the global queue or, with --atlas, assemble the release."""
    global CONTACT
    if args.contact:
        CONTACT = args.contact
    root = Path(args.root)
    cfg = dict(
        seed=42,
        mode="real",
        phase="diagnostic",
        final_approved=False,
        cv_folds=3,
        min_groups=12,
        min_class=15,
        max_dense_cells=20000000,
        test_size=0.20,
        calibration_size=0.20,
        max_geo_uncertainty_m=POINT_UNCERTAINTY_M,
        max_coarse_uncertainty_m=MAX_COARSE_UNCERTAINTY_M,
        max_env_resolution_m=80000,
        decision_threshold=0.5,
        max_missing_fraction=0.5,
        alpha=0.10,
        extra_group_columns=[],
        concentration_unit="mg/g",
        concentration_basis="dry_mass",
        acceptance_ap_gain=None,
        acceptance_min_precision=None,
    )
    if not CONTACT:
        print(
            "Warning: --contact (or ATLAS_CONTACT) is recommended for the User-Agent of public services.",
            flush=True,
        )
    root.mkdir(parents=True, exist_ok=True)

    run_id = (
        datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S")
        + "_"
        + str(time.time_ns())[-6:]
    )
    run_dir = root / "runs" / run_id
    for folder in ["raw", "clean", "quarantine", "reports", "models", "partitions"]:
        (run_dir / folder).mkdir(parents=True, exist_ok=True)
    log = []
    json_write(run_dir / "config.json", cfg)
    import PIL

    versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit-learn": sklearn.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "Pillow": PIL.__version__,
    }
    json_write(run_dir / "versions.json", versions)
    (root / "templates").mkdir(exist_ok=True)
    for table, cols in SCHEMAS.items():
        pd.DataFrame(columns=cols).to_csv(root / "templates" / (table + ".csv"), index=False)
    json_write(root / "templates" / "schema.json", SCHEMAS)
    json_write(root / "templates" / "environment_units.json", ENV_SPEC)
    print("Run:", run_dir)

    # Collection profile: changing these quotas requires a new global queue.
    online = dict(
        taxon_names=[],
        country=None,
        occurrences_per_taxon=30,
        association_rows_per_taxon=20,
        pubchem_max_compounds=5,
        max_environment_points=1,
        max_geo_uncertainty_m=cfg["max_geo_uncertainty_m"],
        use_nasa=True,
        use_soil=True,
        soil_properties=["phh2o", "soc", "clay", "nitrogen"],
        cache_dir=str(root / "public_sources" / "cache"),
        offline=args.offline,
        timeout=25,
    )
    shared = root / "public_sources"
    groot = Path(args.global_root or root / "global_WFO_2026_06")
    global_db = groot / "global_catalog_bulk.sqlite"
    archive_wfo, _ = download_catalog(groot / "catalogue", offline=args.offline)
    print("WFO inventory:", build_registry(archive_wfo, global_db))
    if args.stage != "climate":
        print("Interrupted jobs requeued:", recover_interrupted_jobs(global_db))
    if args.retry_failed:
        print("Failed jobs requeued:", retry_global_jobs(global_db))
    gcfg = dict(batch_size=args.batch_size, time_budget_seconds=args.time_budget)
    lotus, _ = download_catalog(
        groot / "lotus", offline=args.offline, record=LOTUS_RECORD, filename=LOTUS_FILE
    )
    show("LOTUS chemistry", build_lotus(global_db, lotus))
    scan = scan_gbif(
        global_db,
        args.gbif_files,
        k=online["occurrences_per_taxon"],
        offline=args.offline,
        recent=args.gbif_recent,
    )
    show("GBIF snapshot", scan)
    if args.atlas:
        raw = load_selection(global_db, priority_only=args.atlas_scope == "molecules")
        connection_report = dict(
            stage="atlas", scope=args.atlas_scope, inventory=registry_status(global_db)[0]
        )
    elif scan["complete"]:
        soil_dir = shared / "soilgrids_1km"
        try:
            soil_rasters(soil_dir, args.offline)
        except (RuntimeError, FileNotFoundError) as exc:
            print(
                "Incomplete soil maps, SoilGrids 250 m service used for missing properties:",
                str(exc)[:120],
                flush=True,
            )
        bulk = dict(
            collector="bulk",
            lotus=LOTUS_RECORD + "/" + LOTUS_FILE,
            occurrences="gbif snapshot",
            gbif_snapshot=scan["snapshot"],
            occurrences_per_taxon=online["occurrences_per_taxon"],
            max_environment_points=online["max_environment_points"],
            max_geo_uncertainty_m=cfg["max_geo_uncertainty_m"],
            max_coarse_uncertainty_m=cfg["max_coarse_uncertainty_m"],
            climate="open-meteo best_match, nasa power fallback",
            soil="soilgrids 1 km, else wcs 250 m",
            timeout=online["timeout"],
            cache_dir=online["cache_dir"],
            offline=args.offline,
            bulk_db=str(global_db),
            soil_dir=str(soil_dir),
            dem_dir=str(shared / "dem_glo90"),
        )
        if args.stage == "climate":
            print(
                "Climate pass:",
                climate_pass(global_db, bulk, args.time_budget),
                "species processed",
                flush=True,
            )
            connection_report = dict(stage="climate")
        else:
            _, connection_report = run_global_batch(
                global_db,
                dict(bulk, defer_climate=args.stage == "local"),
                gcfg,
                collect_bulk_sources,
                SCHEMAS,
            )
    else:
        # A per-species sample is only final once every snapshot file has been read.
        print(
            f"GBIF snapshot partly read ({scan['files_done']}/{scan['files_total']} files): run again to continue."
        )
        connection_report = dict(
            stage="GBIF scan in progress", gbif=scan, inventory=registry_status(global_db)[0]
        )
    json_write(run_dir / "reports" / "global_coverage.json", connection_report)
    _, ranks, stages = registry_status(global_db)
    ranks.to_csv(run_dir / "reports" / "global_taxonomy_counts.csv", index=False)
    stages.to_csv(run_dir / "reports" / "global_queue_counts.csv", index=False)
    show("Global queue", stages)
    if not args.atlas:
        return run_dir

    if not args.no_elevation:
        env = raw["environment"]
        known = set(env.env_id) if "env_id" in env else set()
        lapse = [r for r in lapse_corrected_temperature(env) if r["env_id"] not in known]
        if lapse:
            env = pd.concat([env, pd.DataFrame(lapse)], ignore_index=True)
        raw["environment"] = env
        check = elevation_check(raw["occurrences"], env)
        if len(check):
            check.to_csv(run_dir / "reports" / "elevation_check.csv", index=False)
            event(
                log,
                "elevation",
                "REVIEW" if check.suspect.any() else "VALIDATED",
                f"Reported elevation ≠ terrain model (> {ELEVATION_QC_TOLERANCE_M} m + precision)",
                int(check.suspect.sum()),
            )
    json_write(run_dir / "reports" / "connection_report.json", connection_report)
    json_write(run_dir / "online_config.json", online)
    manifest = {}
    for name, df in raw.items():
        path = run_dir / "raw" / (name + ".csv")
        df.to_csv(path, index=False)
        manifest[name] = {"rows": len(df), "sha256": sha(path)}
    json_write(run_dir / "input_manifest.json", manifest)
    show("Inputs", pd.DataFrame(manifest).T)

    clean = validate_tables(raw, cfg, run_dir, log)
    show("Quality log", pd.DataFrame(log))
    env = enrich_environment(clean, cfg, run_dir, log)
    chem, geo = build_features(clean, env)
    chem.to_csv(run_dir / "clean" / "chemical_features.csv", index=False)
    geo.to_csv(run_dir / "clean" / "geographic_features.csv", index=False)
    coverage = pd.DataFrame(
        {"chem_missing": chem[NUM_ENV].isna().mean(), "geo_missing": geo[NUM_ENV].isna().mean()}
    )
    coverage.to_csv(run_dir / "reports" / "coverage.csv")
    show("Missing fraction", coverage)
    if len(chem):
        fig, ax = plt.subplots(1, 2, figsize=(11, 3))
        chem.status.value_counts().plot.bar(ax=ax[0], title="Chemical status")
        coverage.plot.bar(ax=ax[1], title="Missing fraction")
        fig.tight_layout()
        fig.savefig(run_dir / "reports" / "coverage.png", dpi=130)
        plt.close(fig)

    db_path = run_dir / "atlas.sqlite"
    database_report = persist_database(
        db_path, clean, chem, geo, pd.read_csv(run_dir / "clean" / "environment_selected.csv")
    )
    json_write(run_dir / "reports" / "database.json", database_report)
    chem, geo = features_from_database(db_path)
    show("Database", pd.DataFrame([database_report]))
    show("Atlas", atlas_summaries(clean, geo, run_dir / "reports"))
    graph = export_graph(clean, run_dir / "reports")
    print(len(graph), "graph edges. Map:", run_dir / "reports" / "atlas.html")
    assert (
        clean["chemistry"]
        .loc[clean["chemistry"].status.isin(["unknown", "not_searched"]), "label"]
        .isna()
        .all()
    )

    shutil.copy(__file__, run_dir / "atlas_plantes.py")
    (run_dir / "GLOBAL_STATE_LOCATION.txt").write_text(
        str(global_db.resolve()) + "\nCache : " + online["cache_dir"], encoding="utf-8"
    )
    json_write(run_dir / "reports" / "quality.json", log)
    report = {
        "run_id": run_id,
        "mode": "real",
        "phase": "diagnostic",
        "scientific_validation": False,
        "warning": "External scientific validation required",
        "inputs": manifest,
        "versions": versions,
        "final_results": {},
        "models": {},
        "concentration": None,
        "limitations": [
            "Capped, non-exhaustive public import",
            "No validated PU or GNN model",
            "Heuristic abstention",
            "Soil 0–30 cm and 30-day weather in this schema",
            "No causal effect established",
            "Regression conditional on exact detections",
            "Online mode: no absences and no local chemical measurements, so no trainable model",
            "DEM 90 m elevation at the reported point: in rugged terrain, 1 km GPS uncertainty can shift it by hundreds of metres",
            "temp_30d_lapse: mean lapse rate 6.5 °C/km, wrong during temperature inversions",
        ],
    }
    json_write(run_dir / "reports" / "run_summary.json", report)
    (run_dir / "requirements-recorded.txt").write_text(
        "\n".join(f"{k}=={v}" for k, v in versions.items() if k != "python"), encoding="utf-8"
    )
    (run_dir / "README.txt").write_text(
        "PLANT ATLAS — real\nSee reports/run_summary.json and config.json.\n"
        "Map: reports/atlas.html. Summary: reports/atlas_summary.html.\n",
        encoding="utf-8",
    )
    shutil.copytree(root / "templates", run_dir / "templates", dirs_exist_ok=True)
    print("\nRelease:", run_dir)
    return run_dir


def parser():
    p = argparse.ArgumentParser(
        description="Plant atlas: WFO taxa, LOTUS natural products, GBIF occurrences and environmental descriptors."
    )
    p.add_argument(
        "--root", default=os.environ.get("ATLAS_ROOT", "atlas_plantes"), help="working directory (local disk)"
    )
    p.add_argument(
        "--global-root", default=None, help="global WFO queue directory (default: <root>/global_WFO_2026_06)"
    )
    p.add_argument(
        "--stage",
        choices=["all", "local", "climate"],
        default="all",
        help="local: everything except weather; climate: weather pass only; all: both",
    )
    p.add_argument(
        "--atlas", action="store_true", help="assemble the release from the collected data (no requests)"
    )
    p.add_argument(
        "--atlas-scope",
        choices=["molecules", "all"],
        default="molecules",
        help="species with natural products only, or all species",
    )
    p.add_argument(
        "--gbif-recent",
        type=int,
        default=5,
        help="most recent occurrences kept per species, in addition to the hashed sample",
    )
    p.add_argument(
        "--gbif-files",
        type=int,
        default=800,
        help="GBIF Parquet files read per launch (about 9,900 in total)",
    )
    p.add_argument("--batch-size", type=int, default=500, help="species processed per launch")
    p.add_argument("--time-budget", type=float, default=3600, help="seconds, checked between species")
    p.add_argument("--retry-failed", action="store_true", help="requeue partial and failed species")
    p.add_argument("--offline", action="store_true", help="use cached responses only")
    p.add_argument(
        "--no-elevation", action="store_true", help="skip elevation-based temperature correction and checks"
    )
    p.add_argument("--contact", default="", help="e-mail or URL added to the User-Agent")
    return p


if __name__ == "__main__":
    run(parser().parse_args())
