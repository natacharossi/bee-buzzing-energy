# First direct quantification of floral handling costs in bees

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16877020.svg)](https://doi.org/10.5281/zenodo.16877020)

Repository containing the processed data, analysis workflows, and figure code for the manuscript:

> **First direct quantification of floral handling costs in bees**  
> Authors: 

**Corresponding author:**  
**ORCID:** 

---

## Repository structure

- **data/**
  - `for_notebooks/` → processed inputs for Jupyter notebooks (Python).
  - `for_R/` → processed inputs for R scripts (figures/tables).
- **notebooks/** → Python notebooks for data processing.
- **R/** → R scripts for statistical analyses and figure generation.
- **figures/** → outputs from the R scripts.

> Figures come from **R scripts**. Notebooks are for **data processing**.

---

## How to reproduce the results

1. **Run the Python notebooks** (in `notebooks/`, Python 3.11.13).  
   Suggested order:  
   1. `synchronising signal detection.ipynb`  
   2. `floral buzz extraction.ipynb`  
   3. `plot CO2 buzzing takeoff.ipynb`  
   4. `metabolic rates computation.ipynb`  
   5. `respiratory quotient analysis.ipynb`  
   6. `stpd corrected metabolic rates.ipynb`  
   7. `compute energy from VCO2 and RQ.ipynb`
   8. `mean energy expenditure and equivalent nectar volume per bee, per behaviour, and per flower species.ipynb`
   9. `volume of nectar required per behaviour (Table 1).ipynb`
   10. `nectar volume computation.ipynb`
   11. `general nectar reference table (Table 2).ipynb`

2. **Run the R scripts** (in `R//`, e.g. in RStudio).  
   These read from `data/for_R/` and produce outputs in `figures/`.

---

## Data availability

Processed data used by the notebooks and R scripts are provided within this repository under `data/for_notebooks/` and `data/for_R/`.

**Data citation**  
 (2025). *Data and code for: First direct quantification of floral handling costs in bees* [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.16877020

---

## Citation

Please cite:
- The manuscript (full reference once available).
- This repository archive: https://doi.org/10.5281/zenodo.16877020

---

## Contact

Questions or issues? Email .
