# -------------------------------------------------------
# Copyright 2023 Novartis Institutes for BioMedical Research, INC
# Author: Andac Demir
# Licensed under the GNUv3 License (the "License");
# you may not use this file except in compliance with the License.
# -------------------------------------------------------
import ast
import logging
import warnings
from collections import Counter
from typing import Dict, List, Optional

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats as stats

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# This will ignore all FutureWarnings
warnings.simplefilter(action="ignore", category=FutureWarning)


class CropSeqLoader:
    def __init__(
        self,
        path: str,
        target_cell_barcodes: List[str],
        control_cell_barcodes: List[str],
        verbose: bool = True,
        preprocessing_params: Optional[Dict[str, float]] = None,
    ):
        """
        Initializes the CropSeqLoader with given parameters.

        Args:
            path (str): Path to the dataset in h5ad format.
            target_cell_barcodes (List[str]): List of barcodes corresponding to target cells.
            control_cell_barcodes (List[str]): List of barcodes corresponding to control cells.
            verbose (bool, optional): If True, will print verbose logs. Defaults to True.
            preprocessing_params (Optional[Dict[str, float]], optional): Dictionary containing preprocessing parameters.
                Expected keys and their default values if not provided:
                    - 'min_genes': 200
                    - 'min_cells': 3
                    - 'counts_per_cell_after': 1e4
                    - 'min_mean': 0.0125
                    - 'max_mean': 3.5
                    - 'min_disp': 0.5
                    - 'max_value': 10
                Defaults to None, in which case the default values mentioned above are used.
        """
        self.verbose = verbose
        self.adata = self.read_dataset(path)
        self.adata_raw = self.adata.copy()
        self.target_cell_barcodes = target_cell_barcodes
        self.control_cell_barcodes = control_cell_barcodes
        self.filter_cells()
        self.preprocessing_params = preprocessing_params or {
            "min_genes": 200,
            "min_cells": 20,
            "counts_per_cell_after": 1e4,
            "min_mean": 0.0125,
            "max_mean": 5,
            "min_disp": -2,
            "max_value": 10,
        }
        self.preprocess_data()
        self.label_cells()

    def read_dataset(self, path: str) -> anndata.AnnData:
        if self.verbose:
            logging.info(f"Reading dataset from {path}...")
        try:
            adata = anndata.read_h5ad(path)
            if self.verbose:
                logging.info("All columns: %s", adata.obs.columns)
            return adata
        except Exception as e:
            logging.error("Error reading the dataset from %s: %s", path, e)
            raise ValueError(f"Error reading the dataset from {path}: {e}")

    def filter_cells(self) -> None:
        """
        Filters the cells based on their barcodes. Removes other cells from the data.
        """
        mask = self.adata.obs.index.isin(self.target_cell_barcodes + self.control_cell_barcodes)
        self.adata = self.adata[mask].copy()

    def label_cells(self) -> None:
        """
        Labels the cells based on their barcodes.
        Target cells are labeled as 1 and control cells as 0.
        """
        self.adata.obs["label"] = np.where(self.adata.obs.index.isin(self.target_cell_barcodes), 1, 0)

    def preprocess_data(self) -> None:
        """
        Preprocesses the data using Scanpy preprocessing functions.

        The steps include:
        - Removing cells with less than a threshold number of genes.
        - Removing genes present in less than a threshold number of cells.
        - Normalizing data to a specific count per cell.
        - Applying logarithmic scaling.
        - Identifying highly variable genes based on mean and dispersion thresholds.
        - Scaling data to have zero mean and unit variance.
        """
        steps = [
            (
                "Remove cells with less than {min_genes} genes.",
                sc.pp.filter_cells,
                {"min_genes": self.preprocessing_params["min_genes"]},
            ),
            (
                "Remove genes present in less than {min_cells} cells.",
                sc.pp.filter_genes,
                {"min_cells": self.preprocessing_params["min_cells"]},
            ),
            (
                "Normalize data to {counts_per_cell_after} counts per cell.",
                sc.pp.normalize_per_cell,
                {"counts_per_cell_after": self.preprocessing_params["counts_per_cell_after"]},
            ),
            # ("Apply logarithmic scaling.", sc.pp.log1p, {}),  # -> data already log normalized
            (
                "Identify highly variable genes with {min_mean} < mean < {max_mean} and dispersion > {min_disp}.",
                sc.pp.highly_variable_genes,
                {
                    "min_mean": self.preprocessing_params["min_mean"],
                    "max_mean": self.preprocessing_params["max_mean"],
                    "min_disp": self.preprocessing_params["min_disp"],
                },
            ),
            ("Regress out unwanted variables.", sc.pp.regress_out, {"keys": ["n_counts", "n_genes", "percent_mito"]}),
        ]

        # Plot histogram of gene expressions
        # Check if it's a sparse matrix and convert to dense if so
        if isinstance(self.adata.X, np.ndarray):
            data_for_histogram = np.ravel(self.adata.X)
        else:  # Assuming it's some sparse matrix format
            data_for_histogram = self.adata.X.toarray().ravel()
        plt.hist(data_for_histogram, bins=1000, color="blue", edgecolor="black")
        plt.title("Histogram of Gene Expressions")
        plt.xlabel("Expression Value")
        plt.ylabel("Number of Occurrences")
        plt.yscale("log")
        plt.show()

        logging.info("\nInitial Data Shape: %s", self.adata.shape)
        logging.info(50 * "#")
        for i, (description, func, kwargs) in enumerate(steps):
            func(self.adata, **kwargs)
            if self.verbose:
                logging.info("Preprocessing Step %d: %s", i + 1, description.format(**kwargs))
                logging.info("Data Shape: %s", self.adata.shape)
                logging.info(50 * "#")

        # Retain only the highly variable genes
        self.adata = self.adata[:, self.adata.var["highly_variable"]]
        if self.verbose:
            print(f"Preprocessing Step {i + 2}: Retain only the highly variable genes.")
            print("Data Shape: ", self.adata.shape)
            print(50 * "#")

        # Scale the data
        sc.pp.scale(self.adata, max_value=self.preprocessing_params["max_value"])
        if self.verbose:
            print(f"Preprocessing Step {i + 3}: Scale the data to have zero mean and unit variance.")
            print(50 * "#")

    def compute_t_statistic(self, gene_name: str) -> float:
        """
        Computes the t-statistic for a specific gene between control and targeted cells.

        Args:
            gene_name (str): The name of the gene to compute the t-statistic for.

        Returns:
            t_statistic (float): The t-statistic value for the difference in mean gene expression between control and
                targeted cells.
        """
        # Find the index of the gene in the data
        gene_index = np.where(self.adata.var_names == gene_name)[0][0]

        # Separate the control and targeted cells
        control_data = self.adata.X[self.adata.obs["label"] == 0, gene_index]
        targeted_data = self.adata.X[self.adata.obs["label"] == 1, gene_index]

        # Compute the mean and standard deviation of each group
        control_mean = np.mean(control_data)
        targeted_mean = np.mean(targeted_data)
        control_std = np.std(control_data, ddof=1)
        targeted_std = np.std(targeted_data, ddof=1)

        # Compute the t-statistic
        t_statistic = (targeted_mean - control_mean) / np.sqrt(
            (control_std ** 2 / len(control_data)) + (targeted_std ** 2 / len(targeted_data))
        )

        return t_statistic

    @staticmethod
    def two_sided_ttest(t_statistic: float, df: int, alpha: float = 0.05) -> bool:
        """
        Uses the t-statistic to determine whether to support or reject the null hypothesis.

        Args:
            t_statistic (float): The t-statistic value for the difference in mean gene expression between control and
                targeted cells.
            df (int): Number of samples which are free to vary: n1 + n2 - 2.
            alpha (float): Significance level.

        Returns:
            t_statistic (float): The t-statistic value for the difference in mean gene expression between control and
                targeted cells.
        """
        # Calculate the p-value using the t_statistic and degrees of freedom (df)
        p_val = stats.t.sf(np.abs(t_statistic), df) * 2

        # Compare the p-value to the significance level (alpha)
        if p_val < alpha:
            logging.info("Null hypothesis is rejected")
            return True
        else:
            logging.info("Null hypothesis is not rejected")
            return False


def main():
    # Randomly pick a perturbed gene (either knocked down or overexpressed using CRISPR)
    de_summary = pd.read_csv("../../data/cropSeq/cropSeq_differential_expression_summary.csv")
    gene = de_summary["gene_perturbed"].sample().iloc[0]
    logging.info("Perturbed gene: %s", gene)
    logging.info(50 * "#")

    # Function to transform string representations of lists into actual lists
    def transform_barcodes(series):
        return series.apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Apply the transformation
    de_summary["target_cell_barcodes"] = transform_barcodes(de_summary["target_cell_barcodes"])
    de_summary["control_cell_barcodes"] = transform_barcodes(de_summary["control_cell_barcodes"])
    de_summary["genes_de"] = transform_barcodes(de_summary["genes_de"])
    de_summary["genes_not_de"] = transform_barcodes(de_summary["genes_not_de"])

    # Identify targeted and non-targeted (control) cell IDs for the experiment this randomly selected gene is modified
    target_cell_barcodes = de_summary[de_summary["gene_perturbed"] == gene]["target_cell_barcodes"].explode().tolist()
    control_cell_barcodes = de_summary[de_summary["gene_perturbed"] == gene]["control_cell_barcodes"].explode().tolist()
    logging.info("Target cell barcodes: %d", len(target_cell_barcodes))
    logging.info("Control cell barcodes: %d", len(control_cell_barcodes))
    logging.info(50 * "#")

    # Identify differentially expressed genes when this randomly selected gene is modified
    genes_de = de_summary[de_summary["gene_perturbed"] == gene]["genes_de"].explode().tolist()
    genes_not_de = de_summary[de_summary["gene_perturbed"] == gene]["genes_not_de"].explode().tolist()
    logging.info("Genes differentially expressed: %d", len(genes_de))
    logging.info("Genes not differentially expressed: %d", len(genes_not_de))
    logging.info(50 * "#")

    loader = CropSeqLoader(
        path="../../data/cropSeq/CRISPRi_filtered_filtered_combined_normalized_selected.h5ad",
        target_cell_barcodes=target_cell_barcodes,
        control_cell_barcodes=control_cell_barcodes,
    )

    logging.info("Shape of preprocessed data: %s", loader.adata.X.shape)
    logging.info("Number of labels: %s", loader.adata.obs["label"].shape)
    num_samples, _ = loader.adata.X.shape
    num_labels = loader.adata.obs["label"].size

    if num_samples != num_labels:
        raise ValueError(f"Mismatch in number of samples and labels: {num_samples} samples but {num_labels} labels.")

    print(dict(Counter(loader.adata.obs["label"])))

    # Compute t-statistic value for hypothesis testing
    t_statistic = loader.compute_t_statistic(gene_name=gene)
    logging.info("t-value: %f", t_statistic)
    loader.two_sided_ttest(t_statistic, df=num_samples - 2, alpha=0.05)


if __name__ == "__main__":
    main()
