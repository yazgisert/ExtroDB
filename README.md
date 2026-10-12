# Extremophile prediction from protein language model embeddings

Code, curated data and results of my Master's thesis in Bioinformatics (TUM / LMU), and
**ExtroDB**, a web resource to browse the curated extremophile proteomes and the predictor models.

- 192 prokaryotic extremophile proteomes (5 classes: thermophile, halophile, acidophile,
  psychrophile, alkaliphile) and 13 eukaryotic algae
- Embeddings: ProtT5, ESM2 650M, ESM2 8M and one-hot (amino acid composition), each with the same small FNN
- Models: 5-class and 12-class phenotype prediction, extremophile vs. non-extremophile (binary),
  6-class, transfer from prokaryotes to eukaryotic algae, homology-based inference as baseline

## Repository structure

| Folder | Content |
|---|---|
| `searchable_resource_app/` | ExtroDB website: Proteome Index (`index.html`) and Predictor Models (`predictors.html`) |
| `searchable_resource/` | Master index of all curated proteomes (data behind the website) |
| `backend/`, `frontend/`, `docker-compose.yml` | Docker setup of the website |
| `curation/` | How the dataset was curated: criteria, source tables, literature screening logs, core taxon tables, sequence retrieval trackers, download manifests |
| `scripts/` | All code, grouped by step (see below) |
| `dataset_splits/` | Split summaries and the protein IDs of the clustered binary / 6-class splits |
| `results/training_results/` | Final result files (metrics with bootstrap CIs) and trained models (`.pt`) |
| `results/plots/` | Final figures |
| `results/` (other folders) | Leakage checks, homology search summaries, InterProScan and TM analyses, exploration |

## Website (ExtroDB)

```bash
docker compose up --build      # first start downloads the proteome FASTAs from Zenodo (~200 MB)
open http://localhost:8080
docker compose down            # stop
```

`--build` is only needed after changes to the website or backend code; otherwise `docker compose up` is enough.
The proteome FASTAs come from Zenodo ([record 21864329](https://zenodo.org/records/21864329)).
The predictor input FASTAs are on Zenodo ([doi:10.5281/zenodo.23298950](https://doi.org/10.5281/zenodo.23298950))
and linked on the Predictor Models page (record ID in `searchable_resource_app/predictors.js`, `ZENODO_RECORD_ID`).

## Pipeline

All scripts are run from the repository root.

| Step | Folder | Main scripts |
|---|---|---|
| 1. Curation and download | `scripts/data_collection/` | `curate_*_bacdive.py`, `initialize_sequence_tracker.py`, `update_sequence_tracker_from_gtdb.py`, `download_nonextremophile_uniprot.py` |
| 2. Master index | `scripts/indexing/` | `build_extremophile_manifests.py`, `build_searchable_extremophile_index.py` |
| 3. Train / val / test splits | `scripts/dataset_splits/`, `scripts/data_collection/` | `prepare_combined_fasta.py` → MMseqs2 `easy-cluster --min-seq-id 0.3` → `assign_splits_mmseqs_representatives.py`; non-extremophiles: `cluster_nonextremophiles.py` → `make_clustered_binary_split.py` |
| 4. Embeddings | `scripts/embeddings/` | `compute_prott5_embeddings.py`, `compute_baseline_embeddings.py` (ESM2 650M / 8M), via the biocentral API |
| 5. Training | `scripts/training/` | `train_local.py` (ProtT5 5/12-class), `train_esm2.py`, `train_onehot.py`, `train_cv.py`, `train_binary_embedder.py`, `train_sixclass_extremophile.py` |
| 6. Evaluation | `scripts/evaluation/` | leakage checks (`check_eukaryote_leakage.py`, `search_test_vs_train.py`, `analyze_test_leakage.py`), algae transfer (`eval_eukaryote_embedders.py`, `eval_eukaryote_fiveclass_filtered.py`), homology-based inference (`eval_homology_*.py`), archaea checks (`eval_binary_by_domain.py`, `eval_fiveclass_by_domain.py`, `acidic_residue_check.py`), THERMO comparison (`run_thermo_*.py`) |
| 7. Figures | `scripts/plotting/` | one script per figure group |
| 8. Website data | `scripts/indexing/` | `export_predictor_resources.py` (predictor FASTAs for Zenodo + results summary for the website) |

## Data that is not in this repository

| Data | Where |
|---|---|
| Proteome FASTAs | Zenodo [record 21864329](https://zenodo.org/records/21864329) |
| Predictor input FASTAs (train / val / test per model) | Zenodo [doi:10.5281/zenodo.23298950](https://doi.org/10.5281/zenodo.23298950) |
| Embeddings (`.h5`) | regenerate with `scripts/embeddings/` |
| Non-extremophile sequences | UniProt Swiss-Prot; the exact sample is listed in `curation/source_tables/nonextremophile_uniprot_sample.tsv` |
| `Pfam-A.hmm` | [Pfam / InterPro](https://www.ebi.ac.uk/interpro/download/pfam/) (for the THERMO comparison) |
| `thermo_profiles.hmm`, BacDive-AI models | THERMO (Koblitz et al. 2025) and [BacDive-AI](https://doi.org/10.1101/2024.08.12.607695) |

### Setup for the THERMO comparison

```bash
pip install pyhmmer
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz
git clone https://github.com/JKoblitz/bacdive-AI.git
cd bacdive-AI && unzip models/thermophile_data.p.zip -d models/ && cd ..

python scripts/evaluation/run_thermo_model.py \
    --pfam Pfam-A.hmm \
    --model bacdive-AI/models/thermophile_data.p \
    --index searchable_resource/extremophile_master_index.tsv \
    --fasta_dir . \
    --out results/training_results/thermo_bacdive_result.json
```

## Environments

- `biocentral` (Python 3.11): training, evaluation, plotting
- `biocentral2` (Python 3.12, `biocentral-api` 2.0): embedding computation via the biocentral server
- MMseqs2 for clustering and homology searches
