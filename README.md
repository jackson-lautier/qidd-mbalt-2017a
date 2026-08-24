<h1 align="center"><project-name></h1>

<p align="center"><project-description></p>

## Introduction

This repository is intended as an online supplement to the working paper:

Lautier, J. P., Chiou, S.H. (2026)
"Testing quasi-independence for discrete, left-truncated data."
(see [https://jacksonlautier.com/publications](https://jacksonlautier.com/publications)
for current working papers)

Please attribute any citations of this repository to the original
manuscript.

This repository includes:

- **raw-data** Scraped loan demographic and performance data from the ABS bond
MBALT 2017-A.  Scraped 2024 Fetal Death Data Files (National Center for Health Statistics).
<i> Please contact the corresponding author (Lautier, J.P.) if interested in the Oxford
Conception Study data. </i>
 

- **data-clean** Cleaned raw MBALT 2017a and fetal death data into files used within the
manuscript. These files are identical to the files created by `data-processing.R'
in the **code** folder.

- **code** First run `data-processing.R` to create the
clean data files in a new folder, **processed-data** (alternatively, rename the
**data-clean** folder as **processed-data**).  Second, all results in the
application section of the manuscript can be replicated with `data-analysis.R`.
All results will either print in the R console or be
stored in a new folder, **results**.

| Coding File | Description |
| --- | --- |
| `data-analysis.R` | Replication code to reproduce the MBALT 2017A application results |
| `data-processing.R` | Replication code to process **raw-data** for `data-analysis.R` |
| `lrt-formulas-LT-only.R` | Functions used for QIDD testing in `data-analysis.R`, left-truncation (LT) only |
| `lrt-formulas-LT-RC.R` | Functions used for QIDD testing in `data-analysis.R`, LT + right-censoring|
| `mbalt-term-time-function.R` | Function used to process the MBALT 2017A raw lease data into survival times |
| `2024-date-ref.csv` | Date file to estimate date of death for `data-processing.R` |
| `fet-death-cols.csv` | Columns references to clean fetal data data with `data-processing.R` |
 



## Lead, Corresponding Author

**Jackson P. Lautier**

- [Website](https://jacksonlautier.com/)

## Complete Authors

**Sy Han Chiou**

- [Website](https://www.sychiou.com/)