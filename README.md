
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

str(Sample.Data) # Simulated sample data retaining the real sample characteristics
#> 'data.frame':    143 obs. of  5 variables:
#>  $ Region: chr  "regione piemonte" "regione piemonte" "regione piemonte" "regione piemonte" ...
#>  $ Beds  : int  183 172 157 179 85 133 189 182 185 71 ...
#>  $ QS    : num  100.8 628.2 147 108.7 23.3 ...
#>  $ Class : chr  "< 200" "< 200" "< 200" "< 200" ...
#>  $ Code  : int  1 2 3 4 5 6 7 8 9 10 ...

str(Reference.Data) # Official list of Italian acute hospitals, updated to 2016
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

  - a candidate list is created from the original sample;
  - the hospitals in the list are permuted randomly or ordered in
    ascending order by ![QS](https://latex.codecogs.com/png.latex?QS
    "QS") (the lower the score, the higher the data quality);
  - the first hospital is selected;
  - all other hospitals belonging to the same block of the selected
    hospital are removed from the candidate list;
  - the process is repeated until there are still available blocks in
    the candidate list;
  - once one hospital from each blocks has been chosen, the hospitals
    from all the blocks are made available again in the list, apart from
    those already selected in the sample;
  - continue until
    ![N\_{required}](https://latex.codecogs.com/png.latex?N_%7Brequired%7D
    "N_{required}") is reached.

For the uniform procedure, we discretized the number of beds only into 4
quantiles since it was less relevant to build a precise discrete
probability distribution.

``` r
subsample.uniform <- function(Input.Sample, n.required, n.quantiles = 4, use.QS = T){
    library(dplyr)
    library(Hmisc)
    
    Hospitals <- Input.Sample %>%
        transmute(Code, Region, Beds = Hmisc::cut2(Beds, g = n.quantiles), QS = if (use.QS) QS else 1) # Prepare data by discretizing continuous variables like the number of acute beds and by changing QS to a fixed value if not to be used
    
    Selected.hospitals <- c()
    Candidates <- data.frame()
    
    for (i in 1:n.required) { # Until n.required is reached..
        if (nrow(Candidates) == 0) { # If the candidate list is empty, rebuilt it from the non-selected hospitals
            Candidates <- Hospitals %>%
                filter(!(Code %in% Selected.hospitals)) %>% # Remove already selected hospitals
                sample_frac() %>% # Permute order, useful only if QS is not used
                arrange(QS)
        }
        
        # Extract the first hospital of the temporary list and add it to the list of selected hospitals
        Extracted.hospital <- Candidates[1,]
        Selected.hospitals <- c(Selected.hospitals, Extracted.hospital$Code)
        
        # Remove from the temporary list all hospitals in the same location/size block of the extracted hospital
        Candidates <- Candidates %>%
            filter(!(Region %in% Extracted.hospital$Region & Beds %in% Extracted.hospital$Beds))
    }
    
    Input.Sample %>% filter(Code %in% Selected.hospitals) # Filter the initial data by the selected hospital codes
}

## Esamples
# Create a subsample of 56 hospitals using the QS
# subsample.uniform(Sample.Data, n.required = 56)

# The same but this time hospitals are chosen randomly
# subsample.uniform(Sample.Data, n.required = 56)
```

## Probability procedure

This algorithm uses information from a population level list (Reference
Data) to built a discrete probability distribution representative of the
target population and then uses it to create a representative
sub-sample. The hospitals are selected according to how representative
is the block they belong to at the country level. The
![QS](https://latex.codecogs.com/png.latex?QS "QS") is used to weight
such representativeness. The weight of the
![QS](https://latex.codecogs.com/png.latex?QS "QS") itself can be
weighted.

  - the blocks are identified in the Reference Data and for each block
    ![p\_{i, country} = Pr(hospital|Region,
    HSize)](https://latex.codecogs.com/png.latex?p_%7Bi%2C%20country%7D%20%3D%20Pr%28hospital%7CRegion%2C%20HSize%29
    "p_{i, country} = Pr(hospital|Region, HSize)") is computed as the
    proportion of hospitals in the block over the total;
  - these probabilities are assigned to the relative blocks in the
    sample;
  - a score is computer for each hospital
    ![j](https://latex.codecogs.com/png.latex?j "j") as ![score\_j =
    p\_{j, country}
    (1-scaled.QS\_j)^w](https://latex.codecogs.com/png.latex?score_j%20%3D%20p_%7Bj%2C%20country%7D%20%281-scaled.QS_j%29%5Ew
    "score_j = p_{j, country} (1-scaled.QS_j)^w") where:
  - ![p\_{j,
    country}](https://latex.codecogs.com/png.latex?p_%7Bj%2C%20country%7D
    "p_{j, country}") is the probability of the block of the hospital
    ![j](https://latex.codecogs.com/png.latex?j "j") at the country
    level;
  - ![scaled.QS\_j](https://latex.codecogs.com/png.latex?scaled.QS_j
    "scaled.QS_j") is the ![QS](https://latex.codecogs.com/png.latex?QS
    "QS") of the hospital ![j](https://latex.codecogs.com/png.latex?j
    "j") after that all ![QS](https://latex.codecogs.com/png.latex?QS
    "QS") have been rescaled to the range \[0,1\], with 1 representing
    the worst quality score and 0 the best.
    ![(1-scaled.QS\_j)](https://latex.codecogs.com/png.latex?%281-scaled.QS_j%29
    "(1-scaled.QS_j)") reweighs the probability of being included of a
    hospital using the quality of the data;
  - ![w](https://latex.codecogs.com/png.latex?w "w") allows scaling the
    importance of the ![QS](https://latex.codecogs.com/png.latex?QS
    "QS") in the selection, with ![w
    = 0](https://latex.codecogs.com/png.latex?w%20%3D%200 "w = 0")
    removing its influence;
  - finally, ![score\_j](https://latex.codecogs.com/png.latex?score_j
    "score_j") is used to order the hospitals and the first
    ![N\_{required}](https://latex.codecogs.com/png.latex?N_%7Brequired%7D
    "N_{required}") get selected. In alternative, the score can be used
    as a weight for selecting the hospitals by random sampling.

<!-- end list -->

``` r

subsample.probability <- function(Input.Sample, Reference.Data, n.required, n.quantiles = 10, QS.weight = 1, method = c('arrange', 'random')){
    
    library(dplyr)
    library(Hmisc)
    library(magrittr)
    library(scales)
    
    method <- match.arg(method)
    
    # Definition of quantiles in the distribution of number of beds according to reference data
    quantiles <- quantile(Reference.Data$Beds, seq(0, 1, length.out = n.quantiles)) %>% round
    
    P_country <- Reference.Data %>% 
        count(Block = Hmisc::cut2(Beds, quantiles) %>% paste('-', Region)) %>% 
        mutate(Prob = n / sum(n)) %>% 
        with(magrittr::set_names(Prob, Block))
    
    Selection <- Input.Sample %>%
        mutate(
            Block = Hmisc::cut2(Beds, quantiles) %>% paste('-', Region), # Identification of the country level blocks in the sample
            Prob = P_country[Block], # Assoction of P_country to the hospital in the sample
            QS.rescale = 1 - scales::rescale(QS), # Creation of the quality weight after rescaling of the QS
            Score = Prob * QS.rescale^QS.weight # Definition of the final score
        ) %>% 
        filter(!is.na(Score)) # Remove blocks that do not appear in the national list
    
    if (method == 'arrange') {
        Selection %>%
            arrange(desc(Score)) %>% 
            head(n.required)
    } else {
        slice_sample(Selection, n = n.required, weight_by = Score)
    }
}
```
