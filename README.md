# Plant atlas pipeline

Code for the dataset described in *From molecules to maps: a reproducible global atlas linking 23,820 plant species, their natural products, and their environments* (Ahiri, submitted to *Biodiversity Informatics*).

The pipeline links accepted species of the World Flora Online Plant List (June 2026) to LOTUS natural-product reports (April 2026), a reproducible sample of the GBIF occurrence snapshot (1 September 2026), SoilGrids 1 km soil descriptors, 30-day weather (Open-Meteo, NASA POWER fallback) and Copernicus DEM GLO-90 elevation.

## Requirements

Python 3.13 and the packages in `requirements.txt`:

```
pip install -r requirements.txt
```

## Usage

Collection is resumable and runs in batches against a global SQLite queue (`GLOBAL`).

```
# 1. Read the GBIF snapshot and collect species from the local copies (repeat until the queue is empty)
python atlas_plantes.py --global-root GLOBAL --stage local --contact you@example.org

# 2. Weather pass (rate-limited public services; repeat until complete)
python atlas_plantes.py --global-root GLOBAL --stage climate --contact you@example.org

# 3. Assemble the release (no network requests)
python atlas_plantes.py --global-root GLOBAL --atlas
```

The release is written to `<root>/runs/<run_id>/`: `clean/` (normalized tables), `quarantine/`, `reports/`, `atlas.sqlite`, the configuration and an input manifest with SHA-256 checksums.

## Relation to the published release

This is a cleaned English edition of the script that produced release `20261005T141832_961774`. Model training, demonstration data and the per-species API collector were removed. Messages, log labels and report names were translated; the data logic is unchanged. Before translation, this edition was run in `--atlas` mode against the original script on the same global database and produced identical tables, database and reports. The translation changes only text literals (and the thousands separator of the summary page).

Report files renamed in this edition:

| Release (original script) | This edition |
|---|---|
| `atlas_especes.csv` | `atlas_species.csv` |
| `atlas_familles.csv` | `atlas_families.csv` |
| `atlas_classes_chimiques.csv` | `atlas_chemical_classes.csv` |
| `atlas_couverture.json` | `atlas_coverage.json` |
| `atlas_resume.html` | `atlas_summary.html` |

Column prefixes changed accordingly (`mediane_` → `median_`, `especes` → `species`, `familles` → `families`).

The verbatim run script is included in the data deposit (`atlas_plantes.py`, SHA-256 `6f40c27200bb9b5ca949c7225bd42f030e3a8714117563739773a5e4bcbb8777`).

## Data sources and licences

WFO Plant List (doi:10.5281/zenodo.20782718), LOTUS (doi:10.5281/zenodo.19360665), GBIF.org, SoilGrids 2.0, Open-Meteo (CC BY 4.0), NASA POWER, Copernicus DEM GLO-90. Occurrence records keep their individual CC0, CC BY or CC BY-NC licences.

## Licence

Code: MIT (see `LICENSE`). Data obtained from the sources above remain under their own licences.
