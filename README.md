# integromics_2 — Exploratory Multi-Omics Integration (PCOS Project)

## Overview

This repository contains the **continuation of my MSc Bioinformatics & Biostatistics final project (TFM)** and represents the **first exploratory phase** of my PhD research.
It focuses on the **integrative analysis of Polycystic Ovary Syndrome (PCOS)** using three complementary omics layers:

* **Microbiome** (amplicon/relative abundance profiles)
* **Blood proteomics**
* **Untargeted/targeted metabolomics**

The goal of this project is **inferential and exploratory multi-omic characterisation**, rather than biomarker discovery.
It documents the early development of modelling strategies, data harmonization, and cross-omics inference that later evolved into more structured frameworks (e.g., *metatest_final*).

Multi-omics factor analysis v2 (MOFA+) reveals specific co-variation patterns in women with PMOS that are strongly influenced by obesity Edmond Géraud-Aguilar , M Ángeles Martínez-García , María Insenser , Susana Barceló-Cerdá , Manuel Luque-Ramírez , Francisco García-García , Héctor F Escobar-Morreal
https://academic.oup.com/hropen/article/2026/4/hoag080/8788611


---

## Main Features

### 1. **Multi-omics data integration with MOFA2**

* Synchronization of sample IDs across microbiome, proteomics and metabolomics
* Basic pre-processing and normalization strategies per omic
* Construction of combined matrices and matched subject-level datasets
* Preliminary batch exploration and variance decomposition

---

### 2. **Exploratory inferential analyses**

This repository includes a broad set of early analytical attempts:

* Univariate statistical comparisons
* Multivariate exploratory analyses (PCA, PLS-based methods)
* Correlation networks across omics layers
* Initial latent factor modelling concepts (pre-MOFA2)
* Cross-omics association scans

These analyses formed the conceptual foundation for later MOFA2-based integration in *metatest_final*.

---

### 3. **Focus and scope**

This project is **not** intended for:

* biomarker discovery,
* supervised classification,
* or predictive modelling.

Instead, its purpose is methodological exploration, hypothesis generation and identification of multi-layer structure in PCOS.

---

## Rationale

This repository represents the **initial integrative research phase** of a larger multi-omics project on PCOS.
It bridges the gap between:

* the MSc TFM exploratory integration, and
* the structured, MOFA2-based multi-omics modelling implemented later in *metatest_final*.

Its value lies in documenting the **conceptual development** of the later, more formal pipelines.

---

## Status

The project is **exploratory** and reflects early PhD research.
It remains a useful reference for:

* inspecting raw cross-omics relationships,
* trialling integration strategies, and
* deriving hypotheses for more advanced models.

---

## Final scripts for reproducible research

* load_and_preprocess.R
* integration.R
* new_Downstream.R
