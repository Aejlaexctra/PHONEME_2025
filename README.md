This is the code repository for the paper **Is phoneme inventory size influenced by speaker population size and isolation? (DOI: XX)**

## Data

Datasets used by in the analysis can be downloaded from **XX (Zenodo Link)**

## Code

All code and analysis is found and can be run within the `Master.R` file.

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

-   A1a - PISa \~ L1_pop, N = 766 Languages

-   A1b - PISc \~ L1_pop, N = 1716

-   A2 - PISc \~ Bordering + Altitude + L1_pop + Island + Range, N = 1710

-   A3a - PISc \~ Bordering + L1_pop + Island + Mainland + Continent + Range, N = 1716

-   A3b - PISc_i \~ Bordering + L1_pop + Island + Mainland + Continent + Range, N = 153

-   A4 - PISc \~ Documentation, N = 1716

-   A5a - PISa_c \~ L1_pop, N = 306

-   A5b - PISc_a \~ L1_pop, N = 306

-   A5c - PISc_a \~ Bordering + L1_pop + Island + Mainland + Continent + Range, N = 306

-   A5d - PISc_a \~ Documentation, N = 306

-   A6 - L1_pop \~ Documentation, N = 1716

Where PISa are phoneme inventories from (1), PISc phoneme inventories from (2), PISc_i are only island endemic languages within PISc, PISa_c are only languages in PISa that are also within PISc, and PISc_a are only languages in PISc that are also within PISa.

### References

(1) Cormac Anderson, Tiago Tresoldi, Simon J Greenhill, Robert Forkel, Russell Gray, Johann-Mattis List, Variation in phoneme inventories: quantifying the problem and improving comparability, Journal of Language Evolution, Volume 8, Issue 2, July 2023, Pages 149–168, <https://doi.org/10.1093/jole/lzad011>

(2) Bromham, L., Yaxley, K.J. & Cardillo, M. Islands are engines of language diversity. Nat Ecol Evol 8, 1991–2002 (2024). <https://doi.org/10.1038/s41559-024-02488-4>
