# BRAF V600E Inhibitor Discovery Pipeline
**Rohs Lab / QBIO 460 | USC, 2024**

## Overview
This repository contains a computational drug discovery pipeline for identifying
potential inhibitors of the BRAF V600E oncogenic mutation — a key driver in
melanoma, thyroid, colorectal, and ovarian cancers.

The project has two layers. The docking pipeline was built as independent research
in the Rohs Lab, running on USC's CARC HPC cluster. The machine learning
classification layer was developed as a final project for QBIO 460 (Intro to
Machine Learning in Biology), directly using the docking results as input data.

---

## Background
The BRAF V600E mutation replaces valine with glutamate at position 600, causing
BRAF to become constitutively active as a monomer and driving unregulated cell
growth. Current FDA-approved inhibitors (Vemurafenib, Dabrafenib, Encorafenib)
are limited by resistance mechanisms including RAF dimerization. This pipeline
screens all three BRAF V600E conformations to more comprehensively identify
candidate inhibitors:

| PDB ID | Conformation |
|--------|-------------|
| 6p7g   | alpha-C helix in / DFG out |
| 4xv2   | alpha-C helix out / DFG in |
| 4e26   | alpha-C helix in / DFG in  |

PDB structures are available at [rcsb.org](https://www.rcsb.org).

---

## Part 1: Molecular Docking Pipeline (Rohs Lab)

### Goal
Screen ChEMBL and ZINC molecules for binding affinity against all three BRAF
V600E conformations using a modular, parallelized docking pipeline on HPC.

### Pipeline Steps

**Step 1 — `prepare_molecules.py`**
Reads ChEMBL SMILES filtered for BRAF V600E activity, enumerates stereoisomers
and conformers per molecule, and saves all conformer molecules to an SDF file.

**Step 2 — `optimize_molecules.py`**
Reads the conformer SDF, applies MMFF force field optimization to each molecule,
clusters conformers by RMSD using Butina clustering, and saves the
cluster-representative molecules for docking.

**Step 3 — `run_docking.py`**
Docks optimized molecules against a receptor pocket using AutoDock QVina.
Runs in parallel using Python multiprocessing. Run once per conformation
by passing the appropriate pocket and crystal ligand as arguments.

```bash
# Example: dock against 6p7g conformation
python run_docking.py \
  --clp crystal_ligand_6p7g.sdf \
  --p pocket_6p7g.pdbqt \
  --conf 20 \
  --op docked_molecules/6p7g \
  --dsf docking_scores.csv
```

**Step 4 — `process_scores.py`**
Takes raw docking scores, retains only the best score per molecule, and
plots binding affinity vs pChEMBL value for validation.

**`validation_test.py`**
Plots the top 100 docking scores comparing known ChEMBL binders vs random
ZINC molecules to validate that the pipeline recovers known good binders.

### Supporting Modules

| File | Purpose |
|------|---------|
| `smiles.py` | ChEMBL SMILES parsing, stereoisomer/conformer generation, SDF writing |
| `molecules.py` | MMFF optimization, crystal ligand centroid calculation, parallel docking |
| `clustering.py` | Butina clustering, RMSD matrix construction, MMFF wrapper |
| `subprocess_calls.py` | OpenBabel and AutoDock QVina subprocess wrappers |

### HPC Job Scripts
All docking jobs were run on USC's CARC HPC cluster using SLURM.
Job scripts are in `job_scripts/` and map directly to each pipeline step.

| Script | Runs |
|--------|------|
| `prepare_molecules.job` | `prepare_molecules.py` |
| `optimize_molecules.job` | `optimize_molecules.py` |
| `run_docking.job` | `run_docking.py` (pass pocket path as argument per conformation) |
| `process_scores.job` | `process_scores.py` |
| `validation.job` | `validation_test.py` |

Intermediate files (conformer SDFs, docked PDBQT files) were stored on the
cluster during execution and are not included in this repository.

---

## Part 2: ML Classification (QBIO 460)

### Goal
Apply machine learning on top of the docking results to classify molecules by
binding conformation and screen the ZINC database for novel candidates.

### Approach
- Morgan fingerprints (2048-bit) generated for all molecules using RDKit
- UMAP dimensionality reduction with hyperparameter search
- KMeans clustering to group molecules by structural similarity, evaluated
  by silhouette score
- 1D CNN trained to classify molecules into four classes: non-binders and
  binders for each of the three BRAF V600E conformations
- Systematic hyperparameter tuning across filter size, filter count, batch
  size, and dense layer dimensions
- Final model: testing AUROC > 0.9, testing accuracy > 0.9

### Results

![Figure 1](ml-classification/figures/figure1_docking_scores.png)
*Figure 1. Top 100 docking scores comparing known ChEMBL binders (blue) vs.
random ZINC molecules (yellow) for 6p7g. ChEMBL molecules dominate the top ranks despite
being outnumbered 20-to-1, validating the docking pipeline.*

![Figure 2](ml-classification/figures/figure2_cnn_training.png)
*Figure 2. CNN training curves showing testing accuracy and AUROC both
stabilizing above 0.9 after approximately 20 epochs.*

The trained CNN was applied to ZINC database molecules with known docking scores,
identifying **27 candidate compounds** predicted to bind the alpha-C helix in /
DFG out conformation (6p7g). Top hit: ZINC000000001547, binding affinity
-9.4 kcal/mol, Class 1 probability 0.86. Full results in `output.csv`.

### Experimental Extension
`vae.ipynb` contains an exploratory implementation of a Variational Autoencoder
trained on molecular fingerprints. As noted in the paper's conclusion, a VAE
could augment the training dataset with structurally diverse samples and improve
CNN performance. This component was not fully completed.

---

## Limitations
- CNN loss plateaued around 0.2; further regularization tuning needed to
  improve loss without overfitting
- ZINC screening was focused on the 6p7g conformation; cross-conformation
  comparison of CNN predictions remains as future work
- Molecular fingerprints may not be the optimal CNN input; alternative
  vectorized molecular representations could improve classification performance

---

## Tools & Libraries
Python, RDKit, TensorFlow/Keras, scikit-learn, UMAP, pandas, numpy, matplotlib,
AutoDock QVina, OpenBabel, MODELLER

**HPC:** USC Center for Advanced Research Computing (CARC)

## Data Sources
- [ChEMBL](https://www.ebi.ac.uk/chembl/) BRAF V600E inhibitor dataset
- [ZINC database](https://zinc.docking.org/) (DBAP subset)
- PDB structures: [6p7g](https://www.rcsb.org/structure/6p7g),
  [4xv2](https://www.rcsb.org/structure/4xv2),
  [4e26](https://www.rcsb.org/structure/4e26)

---

## Repository Structure
```
├── docking-pipeline/
│   ├── prepare_molecules.py      # Step 1: SMILES → stereoisomers → conformers
│   ├── optimize_molecules.py     # Step 2: MMFF optimization + RMSD clustering
│   ├── run_docking.py            # Step 3: parallel QVina docking
│   ├── process_scores.py         # Step 4: score processing + plotting
│   ├── validation_test.py        # ChEMBL vs ZINC validation plot
│   ├── smiles.py                 # SMILES utilities
│   ├── molecules.py              # Docking and optimization functions
│   ├── clustering.py             # Clustering and RMSD functions
│   ├── subprocess_calls.py       # QVina and OpenBabel wrappers
│   └── job_scripts/
│       ├── prepare_molecules.job
│       ├── optimize_molecules.job
│       ├── run_docking.job       # Used for all three conformations
│       ├── process_scores.job
│       └── validation.job
├── ml-classification/
│   ├── final_460.ipynb           # Full ML pipeline: UMAP → CNN → ZINC screening
│   ├── vae.ipynb                 # Experimental VAE extension
│   └── figures/
│       ├── figure1_docking_scores.png
│       └── figure2_cnn_training.png
├── output.csv                    # 27 candidate ZINC compounds with CNN predictions
├── docking_scores.csv            # Raw ZINC docking results
└── paper.pdf                     # Full project writeup (QBIO 460)
```

---

*Docking pipeline developed in the Rohs Lab at USC.
ML classification developed for QBIO 460.*