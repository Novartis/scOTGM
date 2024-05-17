# -------------------------------------------------------
# Copyright 2023 Novartis Institutes for BioMedical Research, INC
# Author: Andac Demir
# Licensed under the GNUv3 License (the "License");
# you may not use this file except in compliance with the License.
# -------------------------------------------------------
import numpy as np
import pandas as pd
from scipy.stats import ranksums


class RankSumTestAnalyzer:
    """
    A class used to analyze gene expression data using the Wilcoxon rank-sum test.

    Methods:
        perform_rank_sum_test(): Performs the Wilcoxon rank-sum test on the datasets.
        rank_genes(): Ranks genes based on the likelihood of being perturbed or knocked down.
    """

    def __init__(self, control_data, target_data):
        """
        Args:
            control_data (pd.DataFrame): A dataframe containing gene expression data for control cells.
            target_data (pd.DataFrame): A dataframe containing gene expression data for target cells.
        """
        self.control_data = control_data
        self.target_data = target_data

    def perform_rank_sum_test(self) -> pd.DataFrame:
        """
        Performs the Wilcoxon rank-sum test on the control and target data.

        Returns:
            pd.DataFrame: A dataframe with genes as rows, and p-values & statistics from the test.
        """
        results = []

        for gene in self.control_data.columns:
            control = self.control_data[gene]
            target = self.target_data[gene]
            stat, p_value = ranksums(control, target)
            results.append({"Gene": gene, "Statistic": stat, "p-value": p_value})

        return pd.DataFrame(results)

    def rank_genes(self) -> pd.DataFrame:
        """
        Ranks genes based on the likelihood of being perturbed or knocked down,
        as indicated by the p-values from the Wilcoxon rank-sum test.

        Returns:
            pd.DataFrame: A dataframe with genes and their ranks based on p-values.
        """
        test_results = self.perform_rank_sum_test()
        ranked_genes = test_results.sort_values("p-value").reset_index(drop=True)
        return ranked_genes
