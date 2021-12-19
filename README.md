# survivalHelp
Helper functions and workflow to run survival analysis and produce plots and statistics

### Installation:
Install necessary dependencies:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSVA") # v 1.38.3
install.packages("survMisc") # v 0.5.5
install.packages("dplyr") # 1.0.7
install.packages("survival") # 3.2
```
Install survivalHelp package and data:
```
devtools::install_github("hmumme/survivalHelp")
```
