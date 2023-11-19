# IPS_code
> **The bulk datasets analysis flow is included in Bulk_analysis.R.**
>
> **The single cell datasets analysis flow is included in GSE149614_annotation.R, GSE151530_annotation.R, singlecell_analysis.R and cellchat.R.**
>
> **The public methods comparison is included in compare_five_models.R**.

- **Bulk_analysis.R:** This contains all the analysis of bulk datasets, details are provided inside the code. The main chapters are as follows: 
  1. Predict immuno-prognostic subtypes in bulk datasets
  2. Survival analysis & the sample distribution of TNM stage & Univariate and multivariate Cox & ROC analysis
  3. The immune landscape
  4. Kegg and Go analysis
  5. Therapy (CTRP2 & Immunotherapy & Sorafenib)
  6. Mutation & CNV analysis
  7. SPP1-related genes bulk analysis (survival analysis & differential analysis & correlation)
- **GSE149614_annotation.R:** This code includes preprocessing, cell annotation and OR ratios for celltypes in single cell data GSE149614. The main chapters are as follows: 
  1. Initial processing HCC_GSE149614 Lu_2022
  2. DoubletFinder (remove double cell)
  3. QC
  4. First clustering
  5. Sample immune subtype
  6. Lymphocyte clustering
  7. Myeloid clustering
- **GSE151530_annotation.R:** This code includes preprocessing, cell annotation and OR ratios for celltypes in single cell data GSE151530. The main chapters are as follows: 
  1. initial processing
  2. The cell categories in the articles were extracted
  3. sample subtype
  4. T+B clustering
  5. Myeloid clustering
- **singlecell_analysis.R:** This code includes construct single cell models and single cell analysis. The main chapters are as follows: 
  1. Construct a single cell model
  2. Single cell subtype classification and functional mapping
  3. Lymphocyte analysis (two single cell datasets)
  4. Myeloid  analysis (two single cell datasets)
  5. Correlations between different single-cell data
  6. infercnv analysis
- **cellchat.R:** This code includes cell communication analysis. The main chapters are as follows: 
  1. GSE149614 cell communication analysis
  2. GSE151530 cell communication analysis
  3. SPP1 expression in single cell datasets
- **compare_five_models.R:** This code compares 29 gene pairs with five published prognostic signatures. The main chapters are as follows: 
  1. method compare (the risk thresholds from the respective dataset)
  2. method compare (the risk thresholds from TCGA as a benchmark)
  3. consistency ratio (from the respective dataset vs from all samples in five datasets merged together)
  4. method compare (the risk thresholds from all samples in five datasets merged together)
- **Supplementary Tables.xlsx is a supplementary tables in article.**

