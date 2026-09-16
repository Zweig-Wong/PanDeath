# <h1 align="center">PanDeath</h1>

<h4 align="center">A streamlined framework for classifying cancer samples into programmed cell death (PCD) subtypes</h4>

<p align="center">
  <img src="Graphical Abstract.png" alt="Graphical Abstract" width="650">
</p>

## Table of Contents

* [Introduction](#introduction)
* [Quick Start](#quick-start)
  * [Setup](#setup)
  * [Input Requirements](#input-requirements)
* [PCD Subtype Prediction Workflow](#pcd-subtype-prediction-workflow)
  * [Step 1: Predict PCD Subtypes](#step-1-predict-pcd-subtypes)
  * [Step 2: Calculate PCD Mode Scores](#step-2-calculate-pcd-mode-scores)
  * [Step 3: Evaluate Concordance with NMF-derived Subtypes](#step-3-evaluate-concordance-with-nmf-derived-subtypes)
* [Publication](#publication)
* [Maintainers](#maintainers)

---

## Introduction

This repository provides an easy-to-run implementation of the pan-cancer programmed cell death (PCD) subtyping framework developed in our study.

Through integrative pan-cancer transcriptomic analysis using non-negative matrix factorization (NMF) encompassing four canonical PCD modes—pyroptosis, necroptosis, ferroptosis, and apoptosis—we identified **five conserved PCD subtypes**:

* **PCD-P** — Pyroptosis-dominant
* **PCD-PN** — Pyroptosis- and necroptosis-dominant
* **PCD-NF** — Necroptosis- and ferroptosis-dominant
* **PCD-FA** — Ferroptosis- and apoptosis-dominant
* **PCD-A** — Apoptosis-dominant

These subtypes are characterized by distinct genomic alterations, cellular stresses, malignant features, and immune microenvironment landscapes. Notably, they are associated with differential prognoses and immunotherapy outcomes.

To facilitate application to independent datasets, we developed a supervised XGBoost-based classifier for robust assignment of samples to the five PCD subtypes. Here we provide a streamlined pipeline for PCD subtype prediction and downstream PCD mode scoring.

---

## Quick Start

### Setup
The pipeline requires:
* **Jupyter Notebook** for running the provided `.ipynb` scripts
* **Python ≥ 3.10**
* **R ≥ 4.0**

Required Python and R packages are specified in the corresponding scripts.


### Input Requirements

The input RNA-seq expression matrix should meet the following requirements:

1. **Tab-separated values (`.tsv`)** format
2. **Samples as rows** and **genes (HUGO symbols) as columns**
3. Expression values should be **TPM-quantified, log2-transformed, and robust-scaled**
4. For datasets containing multiple cancer types, expression values should be **robust-scaled within each cancer type**

> Example PCD expression matrices from the **TCGA**, **CPTAC**, and **pan-immunotherapy (pan-IM)** cohorts used in our study are provided in the `Data/` folder.

---

## PCD Subtype Prediction Workflow

### Step 1: Predict PCD Subtypes

Run: `1-Predict_PCD_subtype.ipynb`

**Outputs:**

* `Prediction.tsv` — Predicted PCD subtype and subtype-specific prediction probabilities for each sample.
* `Prediction_probability.png` — Bar plots showing subtype prediction probabilities for each sample.

### Step 2: Calculate PCD Mode Scores

Run: `2-PCD_score.R`

**Outputs:**

* `PCD_Score.tsv` — ssGSEA scores for pyroptosis, necroptosis, ferroptosis, and apoptosis for each sample.
* `PCD-X_Score.png` — Boxplots comparing PCD mode scores across predicted PCD subtypes.


### **Step 3: Evaluate Concordance with NMF-derived Subtypes**
> **Optional:** This step is intended for internal benchmarking when NMF-derived subtype labels are available.

Run: `3-Concordance.ipynb`

**Outputs:**

* `Performance_overall.tsv` — Overall accuracy, AUC, precision, recall, and F1 score
* `Performance_subtype.tsv` — Performance metrics for individual PCD subtypes
* `Confusion_matrix.png` — Confusion matrix
* `ROC.png` — Receiver operating curves (ROC) for each subtype

---

## Publication

**Manuscript currently under submission.**

> *⚠️ Example datasets and the trained predictive model will be released upon publication.*

Scripts for reproducing the figures presented in the manuscript are provided in the `Figure Scripts/` folder.

## Maintainers

* **Zweig Wong** (GitHub: `Zweig-Wong`)
