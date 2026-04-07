<h1 align="center"><project-name></h1>

<p align="center"><project-description></p>

## Introduction

This repository is intended as an online supplement to the working paper:

Lautier, J. P., Chiou, S.H. (2026) "Testing quasi-independence for discrete data subject to left-truncation."
(see [https://jacksonlautier.com/publications](https://jacksonlautier.com/publications)
for current working papers)

Please attribute any citations of this repository to the original
manuscript.

This repository includes:

- **raw-data** Scraped loan demographic and performance data from the ABS bond
MBALT 2017-A.

- **data-clean** Cleaned raw MBALT 2017a data into files used within the
manuscript. These files are identical to the files created by `data-processing.R'
in the **code** folder.

- **code** First run `data-processing.R` to create the
clean data files in a new folder, **processed-data** (alternatively, rename the
**data-clean** folder as **processed-data**).  Second, all results in the
application section of the manuscript can be replicated with `data-analysis.R`.
All results will either print in the R console or be
stored in a new folder, **results**.


%## Lead, Corresponding Author

%**Jackson P. Lautier**

%- [Website](https://jacksonlautier.com/)

%## Complete Authors

%**Sy Han Chiou**

%- [Website](https://www.sychiou.com/)