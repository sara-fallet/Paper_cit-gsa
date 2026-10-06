

# cit_gsa: gene set analysis for scRNA-seq with covariate adjustment

With this repository you can reproduce the simulations and the figures of the article: "*cit_gsa: gene set analysis for scRNA-seq with covariate adjustment*"

## Organisation 

```graphql
.
└── paper_cit-gsa/
    ├── R    
    ├── data
    ├── results
    ├── README.md
    └── paper_cit-gsa.Rproj
```


## Description

### R 

This file contain the R script to 

- Reproduce the simulations with the slurm files to run them in a computational cluster: file `Simulations`
- Reproduce the application on Villani data : `Application_Villani.R`
- Reproduce the figures that are in the paper: `Generate_Figures.R`


### data

This file contain the data used for the application on [Villani et al.](https://www.science.org/doi/10.1126/science.aah4573) data. The raw data can be found on GEO [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820). 

These data are used in the R script `R/Application_Villani.R`
 

### results

This file contain the results for the simulations and the application analysis. 

For the application: 

- The results of the DEA analysis perform in `R/Application_Villani.R`: `Villani_DEA_BTM_adj_cd83.RData` and `Villani_DEA_BTM_no_adj_cd83.RData`
- The results obtained after data preprocessing with `R/Prepare_Data.R`: `Villani_DEA_BTM_pvalues_combine_cd83.RData`

For the simulations:

- In the file `Raw_Pvalues_Simulations`: all the raw outputs for each scenario of simulation obtain with the computational cluster [CURTA](https://www.mcia.fr)
- The computation of the True Discovery Rate for each methods obtain with the script `R/Prepare_Data.R`: `Simulation_TDR_DE.RData` and `Simulation_TDR_DM.RData`

The file `Figures` contain the plots that are in the manuscript generated with `R/Generate_Figures.R` 
