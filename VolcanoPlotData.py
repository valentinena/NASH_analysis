# Modified from kevincoakley metabolomicsworkbench/jupyter-notebooks/MWPerformVolcanoPlotAnalysis.ipynb

import pandas as pd
import numpy as np

import scipy
import statsmodels.stats.multitest

def GenerateVolcanoPlotData(feature_table, metadata, condition_col, condition1, condition2):
        
        """Generate data for volcano plots."""

        # Extract data for two specified classes...
        DataA_subject_list = metadata.loc[metadata[condition_col] == condition1].index.to_list()
        DataB_subject_list = metadata.loc[metadata[condition_col] == condition2].index.to_list()
        DataA = feature_table[DataA_subject_list].T
        DataB = feature_table[DataB_subject_list].T


        # Tranform data...
        DataA = np.log2(DataA)
        DataB = np.log2(DataB)

        # Calculate statistics...
        TStatistics, PValues = scipy.stats.ttest_ind(DataA, DataB, equal_var = False, nan_policy = 'omit')
        TStatistics = TStatistics.tolist()
        PValues = PValues.tolist()

        # Adjust P-values...
        Rejects, AdjustedPValues = statsmodels.stats.multitest.fdrcorrection(PValues, alpha=0.05, method='indep', is_sorted=False)
        AdjustedPValues = AdjustedPValues.tolist()

        # Calculate fold change...
        MeanA = DataA.mean(axis = 0)
        MeanB = DataB.mean(axis = 0)
        Log2FoldChange = MeanB.subtract(MeanA)
        Log2FoldChange = Log2FoldChange.tolist()

        # Cast np array to a list...
        Log10PValues = -np.log10(PValues)
        Log10PValues = Log10PValues.tolist()

        VolcanoPlotDataFrame = pd.DataFrame(
            [Log2FoldChange, PValues, Log10PValues, AdjustedPValues, TStatistics],
            columns = DataA.columns,
            index=['log2(FoldChange)', 'P-value', '-log10(P-value)', 'AdjustedP-value', 't-Statistic']).transpose()
    
        return VolcanoPlotDataFrame
