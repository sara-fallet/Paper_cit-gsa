

# cit_gsa: gene set analysis for scRNA-seq with covariate adjustment

With this repository you can reproduce the simulations and the figures of the article: "*cit_gsa: gene set analysis for scRNA-seq with covariate adjustment*"

## Organisation 

```graphql
.
└── paper_cit-gsa/
    ├── Main_Generate_Figures.R
    ├── R
    ├── figures
    ├── data
    ├── results
    ├── README.md
    └── paper_cit-gsa.Rproj
```


## Description

### Script `Main_Generate_Figures.R`

To reproduce all the figures of the paper.

### R 

This file contains the R script to reproduce the results of:

- The simulations with the slurm files to run them in a computational cluster: file `Simulations`
- The application on Villani data : `Application_Villani.R`

And it contains the scripts to reproduce the plot of the main manuscript and the supplementary materials.


### figures

This file contains the figures that are in the paper. 

### data

This file contains the data used for the application on [Villani et al.](https://www.science.org/doi/10.1126/science.aah4573) data. The raw data can be found on GEO [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820). 

These data are used in the R script `R/Application_Villani.R`.

It also contains the gene sets database: the BTM of Stanford. 
 

### results

This file contain the results for the simulations and the application analysis. 

For the application on CD83: 

- The results of the DEA analysis perform in `R/Application_Villani.R`: `Villani_DEA_BTM_adj_cd83.RData` and `Villani_DEA_BTM_no_adj_cd83.RData`
- The results obtained after data preprocessing with `R/Prepare_Data.R`: `Villani_DEA_BTM_pvalues_combine_cd83.RData`

There is the same for IDO analysis.

For the simulations:

- In the file `Raw_Pvalues_Simulations`: all the raw outputs for each scenario of simulation obtain with the computational cluster [CURTA](https://www.mcia.fr)
- The computation of the True Discovery Rate for each methods obtain with the script `R/Prepare_Data.R`: `Simulation_TDR_DE.RData` and `Simulation_TDR_DM.RData`

