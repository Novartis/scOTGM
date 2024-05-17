# -------------------------------------------------------
# Copyright 2023 Novartis Institutes for BioMedical Research, INC
# Author: Andac Demir
# Licensed under the GNUv3 License (the "License");
# you may not use this file except in compliance with the License.
# -------------------------------------------------------
from typing import Tuple
import numpy as np


class InSilicoPerturbationPredictor:
    """
    Predicts the effects of in-silico perturbations on gene expression using the covariance matrix
    derived from a perturbation model.
    """
    def __init__(self, perturbation_model, method: str = "hit and run mcmc") -> None:
        """
        Initializes the predictor with a perturbation model and the method to calculate additive parameters.

        Parameters:
            perturbation_model: An object providing the `additive_parameters` method for computing
                the perturbation's mean and covariance.
            method (str): The method used by the perturbation model to estimate additive parameters.
                Available options: "monte carlo" and "hit and run mcmc". Defaults to "hit and run mcmc".
        """
        self.pm = perturbation_model
        self.method = method

    def calculate_delta_xj(self, i: int, j: int, delta_xi: float) -> float:
        """
        Calculates the change in expression of gene j (delta xj) due to a perturbation in gene i (delta xi),
        based on the covariance matrix.

        Parameters:
            i (int): Index of gene i whose expression is perturbed.
            j (int): Index of gene j for which the change in expression is predicted.
            delta_xi (float): Magnitude of perturbation in gene i's expression.

        Returns:
            float: Predicted change in expression of gene j (delta xj).
        """
        perturbation_cov = self.pm.additive_parameters(method=self.method)[1]
        cov = self.pm.pca.components_.T @ (perturbation_cov * self.pm.pca.explained_variance_[:, np.newaxis]) @ self.pm.pca.components_

        numerator: float = cov[i, j] * delta_xi
        denominator: float = cov[i, i]

        if denominator > 0:
            return numerator / denominator
        return 0.0
