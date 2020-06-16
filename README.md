
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Code implementation of the subsampling algorithms.

Acronyms: healthcare associated Infections (HAI), antimicrobial (AM),
microorganism (MO), alcoholic hand rub (AHR), point prevalence study
(PPS), likelihood ratio test (LRT).

## Introduction

We implemented 3 algorithms to produce a sample with specific
distributional properties starting from a large convenience sample. We
defined these methods as: “Distance procedure”, “Probability procedure”,
and the “Uniformity procedure”. In the following document we present the
R code to implement the algorithms. The specific implementation
described here is focused on the Italian PPS study on HAI/AMU prevalence
in acute care hospitals, utilizing hospital size and location (Italian
region) as stratifying variables, but can be easily extended to more
general cases. Each function take as input the sample data and the
reference data whose distribution should be reproduced, plus some
procedure specific arguments. The full code is available in the R folder
of this repository together with simulated data based on the Italian PPS
convenience sample and the official list of Italian Hospitals. The data
show the number of beds in the hospital as a number and as size class,
the region of belonging, and an ID. The sample data also contains
simulated Quality Scores (QS).

## Probability Procedure
