### Metabolomic analysis of nonalcoholic steatohepatitis (NASH)

- Workflow for univariate analysis
    - Perform volcano plot analysis on metabolomics data
        Need:
        1) Metabolic feature table (e.g. LC-MS/MS features)
        2) Experimental metadata

Project Notes:
Metabolomics β-diversity and α-diversity were calculated using QIIME2. β-diversity measures were compared using pairwise PERMANOVA. For volcano plots, statistical significance was assessed using a t-test (a = 5%), then corrected for multiple hypothesis testing with the FDR method. Significance thresholds for significantly increased and decreased metabolites were set at +/-1 Log2(FoldChange) and p < 0.05.
