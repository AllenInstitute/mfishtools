# mfishtools

R functions for gene selection and analysis of mFISH data

mfishtools includes many functions that are used for analysis of data for the CZI SpaceTx 
project, and mostly relies on correlation-based analysis with filtering.

Install using:
```
devtools::install_github("AllenInstitute/mfishtools",auth_token="802976690281f1483c40de46d0a07e9d01a3de08")
```

### Library use cases

There are two primary use cases for this libary, both of which will soon have associated vignettes.

1. **Building a combinatorial marker gene panel for spatial transcriptomics.**  This allows the generation of computationally "optimal" marker gene panels based on single cell/nucleus RNA-Seq reference data.  A starting set of manually-selected marker genes is first selected, and then the remaining genes are chosen using a greedy algorithm.  Relevant statistics and plots are generated that show the predicted success for the panel.  
2. **Mapping cells from spatial transcriptomics data sets to reference cell types.** This allows for cell type calling of cells in a spatial transcriptomics study, and also predicts the accuracy of the calls based on reference data.  *Note: it is currently unclear how reliable this method is at correctly predicting cell type calls.  Please review results carefully!*  Plots can also be generated to show the results.  

This README will be updated with some more information soon.  For now, please use the R help ("?") and read the vignettes if additional information is needed, or e-mail me at jeremym@alleninstitute.org.
