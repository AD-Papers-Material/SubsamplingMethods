
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Code implementation of the subsampling algorithms.

## Introduction

We implemented 3 algorithms, defined “Distance procedure”, the
“Probability procedure”, and the “Uniformity procedure” with the aim
of producing sub-samples with specific distributional properties
starting from an initial, collection of data which is possibly biased in
terms of generalizability to a target population. In the following
document we present the R code to implement the algorithms. The specific
implementation described here is focused on the Italian PPS study on
HAI/AMU prevalence in acute care hospitals but can easily translated to
different settings.

The methods try to change the distributional characteristics of a sample
by producing a sub-sample whose units are chosen according to some
characteristics of interest. Two of these methods, the Distance
procedure and the Probability procedure, are aimed at generating a
representative sample of a target population, by using reference data at
the population level with information on the distribution of relevant
characteristics of the observational units. The Uniformity procedure
instead generate a sample which is uniform in relation to such
characteristics, therefore it does not require population reference
data.

This document and all the relative material is available at
<https://github.com/AD-Papers-Material/SubsamplingMethods>.

## Input data

The algorithms take as input a database with a statistical unit for each
row and the characteristics of interest as columns. Furthermore, a
column with an unique ID for each unit is useful for post-hoc checks and
mandatory for the Distance procedure. Our algorithms also accept the
Quality Score (![QS](https://latex.codecogs.com/png.latex?QS "QS"), cfr.
Methods) as additional characteristics, but they can be easily modified
to remove such feature.

The Probability and the Distance procedures also requires population
level reference data, with the same structure. If individual level data
is not available at the population, simulated data can be generated
given access to the joint distribution of the considered
characteristics, ensuring that the simulated dataset is large enough to
limit random variation.

Our test case data use acute care hospitals as observational units and
hospital size (number of acute care beds) and region of location as
characteristics of interest. As a reference, the hospital are also
grouped in three hospital size categories; the same categorization is
used in the manuscript. We provide simulated sample data for the testing
of the procedures at
<https://github.com/AD-Papers-Material/SubsamplingMethods>.

``` r

str(SampleData) # Simulated sample data retaining the real sample characteristics
#> 'data.frame':    143 obs. of  5 variables:
#>  $ Region: chr  "regione piemonte" "regione piemonte" "regione piemonte" "regione piemonte" ...
#>  $ Beds  : int  183 172 157 179 85 133 189 182 185 71 ...
#>  $ QS    : num  100.8 628.2 147 108.7 23.3 ...
#>  $ Class : chr  "< 200" "< 200" "< 200" "< 200" ...
#>  $ Code  : int  1 2 3 4 5 6 7 8 9 10 ...

str(ReferenceData) # Official list of Italian acute hospitals, updated to 2016
#> 'data.frame':    963 obs. of  4 variables:
#>  $ Code  : int  10007 10010 10012 10612 10653 10655 10003 10011 10013 10611 ...
#>  $ Beds  : int  258 73 22 105 9 96 337 368 95 115 ...
#>  $ Class : chr  "200 - 500" "< 200" "< 200" "< 200" ...
#>  $ Region: chr  "regione piemonte" "regione piemonte" "regione piemonte" "regione piemonte" ...
```

## General aspects and notation

The procedures at the moment can utilize only categorical
characteristics/features of a sample. Therefore continuous
characteristics (specifically, hospital size in this case) are
discretized in quantiles by the algorithms before use. The number of
quantiles to split continuous features into is an input to the
algorithms. The hospitals are then grouped into *blocks* according to
the characteristics of interest: in this study, we used
*location/hospital size blocks*, defined by the Italian region
(![Region](https://latex.codecogs.com/png.latex?Region "Region")) and
the quantile of number of acute beds
(![HSize](https://latex.codecogs.com/png.latex?HSize "HSize")) the
hospitals fall into. The region and the hospital size quantile allow the
definition of a joint discrete probability distribution which is used by
the algorithms. For each block
![block\_i](https://latex.codecogs.com/png.latex?block_i "block_i") a
probability ![p\_{i} = Pr(hospital|Region,
HSize)](https://latex.codecogs.com/png.latex?p_%7Bi%7D%20%3D%20Pr%28hospital%7CRegion%2C%20HSize%29
"p_{i} = Pr(hospital|Region, HSize)") is defined, either given a sample
(![p\_{i,sample}](https://latex.codecogs.com/png.latex?p_%7Bi%2Csample%7D
"p_{i,sample}")) or the whole country
(![p\_{i,country}](https://latex.codecogs.com/png.latex?p_%7Bi%2Ccountry%7D
"p_{i,country}")) which indicate the fraction of hospital in a block
over the total. All the algorithms are costrained by a parameter
![N\_{required}](https://latex.codecogs.com/png.latex?N_%7Brequired%7D
"N_{required}") which defines the size of the final sub-sample. As
mentioned above, the algorithms may use the Quality Score
![QS](https://latex.codecogs.com/png.latex?QS "QS") as further
discriminant in the selection; to avoid using the
![QS](https://latex.codecogs.com/png.latex?QS "QS") without code
modification is sufficient to assign the same value (a positive number)
to all hospitals. Note that lower
![QS](https://latex.codecogs.com/png.latex?QS "QS") implies better data
quality.

## Unformity procedure

This procedure sub-samples hospitals trying to obtain an equal
proportion of hospitals in every block, by iteratively choosing one
hospital from each block. This is the general implementation:

  - hospitals in the sample are permuted or ordered by
    ![QS](https://latex.codecogs.com/png.latex?QS "QS");
  - the first hospital (random or with the lower, better,
    ![QS](https://latex.codecogs.com/png.latex?QS "QS")) is selected;
  - all other hospitals belonging to the same block are removed (the
    block is not available anymore);
  - a new hospital is chosen again randomly or by
    ![QS](https://latex.codecogs.com/png.latex?QS "QS") for the
    remaining blocks;
  - once one hospital from each blocks has been chosen all the blocks
    are made available again;
  - continue until
    ![N\_{required}](https://latex.codecogs.com/png.latex?N_%7Brequired%7D
    "N_{required}") is reached.

For the uniform procedure, we discretized the number of beds only into 4
quantiles since it was less relevant to build a precise discrete
probability distribution.

``` r
subsample.uniform <- function(InputSample, n.required, n.quantiles = 4, use.QS = T){
    library(dplyr)
    library(Hmisc)
 
    Hospitals <- InputSample %>%
        transmute(Code, Region, Beds = Hmisc::cut2(Beds, g = n.quantiles), QS = if (use.QS) QS else 1) # Prepare data by discretizing continuous variables like the number of acute beds and by changing QS to a fixed value if not to be used

    Selected.hospitals <- c()
    Hospitals.temp <- data.frame()

    for (i in 1:n.required) { # Until n.required is reached..
        if (nrow(Hospitals.temp) == 0) { # If no hospitals are already selected store all sample hospital here
            Hospitals.temp <- Hospitals %>%
                filter(!(Code %in% Selected.hospitals)) %>% # Remove already selected hospitals
                sample_frac() %>% # Permute order, useful only if QS is not used
                arrange(QS)
        }

        # Extract the first hospital of the temporary list and add it to the list of selected hospitals
        Extracted.hospital <- Hospitals.temp[1,]
        Selected.hospitals <- c(Selected.hospitals, Extracted.hospital$Code)

        # Remove from the temporary list all hospitals in the same location/size block of the extracted hospital
        Hospitals.temp <- Hospitals.temp %>%
            filter(!(Region %in% Extracted.hospital$Region & Beds %in% Extracted.hospital$Beds))
    }

    InputSample %>% filter(Code %in% Selected.hospitals) # Filter the initial data by the selected hospital codes
}

## Esamples
# Create a subsample of 56 hospitals using the QS
# subsample.uniform(SampleData, n.required = 56)

# The same but this time hospitals are chosen randomly
# subsample.uniform(SampleData, n.required = 56)
```
