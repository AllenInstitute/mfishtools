# mfishtools

R functions for gene selection and analysis of mFISH data

mfishtools includes many functions that are used for analysis of data for the CZI SpaceTx project, and mostly relies on correlation-based analysis with filtering.  *This library is still in beta testing and may be buggy.  Community involvement is encouraged through both issues and pull requests.*  

## Installation

Install prerequisites:
```
install.packages("devtools")
devtools::install_github("AllenInstitute/scrattch.vis")
devtools::install_github("AllenInstitute/tasic2016data")
```

Note that some people may need to manually install the `GO.db` and `WGCNA` libraries as well:
```
install.packages("BiocManager")
BiocManager::install("GO.db")
BiocManager::install("WGCNA")
```

Install `mfishtools` using:
```
# Quickly, but without the vignettes:
devtools::install_github("AllenInstitute/mfishtools")

# More slowly, but with the vignettes:
install.packages("remotes", repos='http://cran.us.r-project.org')
remotes::install_github("AllenInstitute/mfishtools", build_vignettes = TRUE)
```


## Library use cases

There are two primary use cases for this libary:

1. **Building a combinatorial marker gene panel for spatial transcriptomics.** [LINK TO VIGNETTE](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/inhibitory_marker_selection.html)  This allows the generation of computationally "optimal" marker gene panels based on single cell/nucleus RNA-Seq reference data.  A starting set of manually-selected marker genes is first selected, and then the remaining genes are chosen using a greedy algorithm.  Relevant statistics and plots are generated that show the predicted success for the panel.  
2. **Mapping cells from spatial transcriptomics data sets to reference cell types.** [LINK TO VIGNETTE](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/inhibitory_marker_mapping.html)  This allows for cell type calling of cells in a spatial transcriptomics study, and also predicts the accuracy of the calls based on reference data.  *Note: it is currently unclear how reliable this method is at correctly predicting cell type calls.  Please review results carefully!*  Plots can also be generated to show the results.  

Many functions are currently not included in these vignettes; please use the R help ("?") if additional information is needed, or e-mail me at jeremym@alleninstitute.org.

## License

The license for this package is available on Github at: https://github.com/AllenInstitute/mfishtools/blob/master/LICENSE

## Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.

## Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/mfishtools/blob/master/CONTRIBUTION
