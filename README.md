This is the code repository for the paper **Is phoneme inventory size influenced by speaker population size and isolation? (DOI: XX)**

## Data

Datasets used by in the analysis are included in `data` folder within the repository. Phylogenetic and spatial distance matrices, found in `data/distance_matrices` are subsets of matrices from (3). The PHOIBLE phoneme inventory from (1) is found in `data/phoible`. Data for the second phoneme inventory `data/NEE22.csv` is taken from (5), and data for the predictors used in this study `data/NEE24.csv` are from (2) and (4).

## Code

All code and analysis is found and can be run within the `Master.R` file. Due to the time it takes to run the method, weightings of phylogenetic and spatial distance matrices have already been calculated (included in the folder `data`), but can be recalculated by uncommenting the relevant sections (found with the use of the function `best_p` ).

### Script Layout

There are nine sections within the file:

1.  Loading required datasets and packages

2.  Function definitions

3.  Dataset cleaning and formatting

4.  Exporting datasets to disk

5.  Predictor and response variable transformations

6.  Dataset visualisation

7.  OLS and GLS Model Fitting

8.  Result Plotting

9.  References

### Analyses 

There are 11 rounds of analysis (see paper for overview of predictor and response variables):

#### Section 1

-   A1a - PISa \~ L1_pop, N = 766 Languages

-   A1b - PISc \~ L1_pop, N = 1716

#### Section 2

-   A2 - PISc \~ Bordering + Altitude + L1_pop + Island + Range, N = 1710

#### Section 3

-   A4 - PISc \~ Documentation, N = 1716

#### Supplementary / Not Included

-   A3a - PISc \~ Bordering + L1_pop + Island + Mainland + Continent + Range, N = 1716

-   A3b - PISc_i \~ Bordering + L1_pop + Island + Mainland + Continent + Range, N = 153

-   A5a - PISa_c \~ L1_pop, N = 306

-   A5b - PISc_a \~ Bordering + Altitude + L1_pop + Island + Range, N = 306

-   A5c - PISc_a \~ Bordering + L1_pop + Island + Mainland + Continent + Range, N = 306

-   A5d - PISc_a \~ Documentation, N = 306

-   A6 - L1_pop \~ Documentation, N = 1716

Where PISa are phoneme inventories from (1), PISc phoneme inventories from (5), PISc_i are only island endemic languages within PISc, PISa_c are only languages in PISa that are also within PISc, and PISc_a are only languages in PISc that are also within PISa.

### References

(1) Anderson, C., Tresoldi, T., Greenhill, S. J., Forkel, R., Gray, R., & List, J.-M. (2023). Variation in phoneme inventories: Quantifying the problem and improving comparability. Journal of Language Evolution, 8(2), 149–168. https://doi.org/10.1093/jole/lzad011

(2) Bromham, L., Yaxley, K. J., & Cardillo, M. (2024). Islands are engines of language diversity. Nature Ecology & Evolution, 8(10), 1991–2002. https://doi.org/10.1038/s41559-024-02488-4

(3) Hua, X., Greenhill, S. J., Cardillo, M., Schneemann, H., & Bromham, L. (2019). The ecological drivers of variation in global language diversity. Nature Communications, 10(1), 2047. https://doi.org/10.1038/s41467-019-09842-2

(4) Bromham, L., Dinnage, R., Skirgård, H., Ritchie, A., Cardillo, M., Meakins, F., Greenhill, S., & Hua, X. (2022). Global predictors of language endangerment and the future of linguistic diversity. Nature Ecology & Evolution, 6(2), 163–173. https://doi.org/10.1038/s41559-021-01604-y

(5) Creanza, N., Ruhlen, M., Pemberton, T. J., Rosenberg, N. A., Feldman, M. W., & Ramachandran, S. (2015). A comparison of worldwide phonemic and genetic variation in human populations. Proceedings of the National Academy of Sciences, 112(5), 1265–1272. https://doi.org/10.1073/pnas.1424033112

