# -------------------------------------------------------
# Copyright 2023 Novartis Institutes for BioMedical Research, INC
# Author: Andac Demir
# Licensed under the GNUv3 License (the "License");
# you may not use this file except in compliance with the License.
# -------------------------------------------------------
import numpy as np
from scipy.stats import multivariate_normal, norm
import matplotlib.pyplot as plt


def generate_random_params(dimension: int) -> tuple:
    """
    Generate random means and covariances for a given dimension.

    Args:
    dimension (int): The dimensionality of the mean and covariance.

    Returns:
        tuple: A tuple containing the mean vector and covariance matrix.
    """
    mean = np.random.randn(dimension)
    random_matrix = np.random.randn(dimension, dimension)
    # Symmetric, positive-definite matrix
    cov = np.dot(random_matrix, random_matrix.T)
    # Small positive value to the diagonal for numerical stability
    cov += np.eye(dimension) * 1e-6
    return mean, cov


def generate_non_independent_data(
    mean1: np.ndarray,
    cov1: np.ndarray,
    mean2: np.ndarray,
    cov2: np.ndarray,
    num_samples: int,
    correlation_factor: float,
) -> tuple:
    base_samples = np.random.multivariate_normal(mean1, cov1, num_samples)
    noise = np.random.multivariate_normal(mean2, cov2, num_samples)
    correlated_samples = base_samples * correlation_factor + noise * (
        1 - correlation_factor
    )
    return base_samples, noise, correlated_samples


def get_non_zero_density_bounds(
    cov: np.ndarray,
    point: np.ndarray,
    direction: np.ndarray,
    confidence: float = 0.95,
) -> tuple:
    """
    Calculate the bounds for non-zero density in a specified direction.

    Args:
        cov (np.ndarray): Covariance matrix of the distribution.
        point (np.ndarray): Current point in the distribution.
        direction (np.ndarray): Direction vector for calculating bounds.
        confidence (float, optional): Confidence level for the bounds.
            Default is 0.95.

    Returns:
        tuple: A tuple of lower and upper bounds.
    """
    # Normalize the direction vector
    direction_normalized = direction / np.linalg.norm(direction)

    # Project the point onto the direction
    point_projection = np.dot(point, direction_normalized)

    # Calculate the variance of the projection
    variance_projection = np.dot(
        direction_normalized.T, np.dot(cov, direction_normalized)
    )
    # Calculate the standard deviation of the projection
    std_deviation = np.sqrt(variance_projection)
    # Find the z-scores for the specified confidence interval
    z_score = norm.ppf((1 + confidence) / 2)

    # Calculate the bounds
    lower_bound = point_projection - z_score * std_deviation
    upper_bound = point_projection + z_score * std_deviation
    return lower_bound, upper_bound


def estimate_cross_covariance(
    mean1: np.ndarray,
    cov1: np.ndarray,
    mean2: np.ndarray,
    cov2: np.ndarray,
    ground_truth: np.ndarray,
    num_iterations: int = 1000,
) -> tuple:
    """
    Estimate the cross-covariance matrix using hit-and-run Gibbs sampling and
    calculate the RMSE and Frobenius norm at each iteration.

    Args:
        mean1, mean2 (np.ndarray): Mean vectors of the first and
            second distributions.
        cov1, cov2 (np.ndarray): Covariance matrices of the first and
            second distributions.
        ground_truth (np.ndarray): Ground truth cross-covariance matrix.
        num_iterations (int, optional): Number of iterations for sampling.
            Default is 1000.

    Returns:
        tuple: Estimated cross-covariance matrix, list of RMSE values,
            and list of Frobenius norms.
    """
    rmse_values, frobenius_norms = [], []

    # Random initialization of cross-covariance matrix
    dimension = mean1.shape[0]
    random_matrix = np.random.randn(dimension, dimension)
    cross_covariance = (
        np.dot(random_matrix, random_matrix.T) + np.eye(dimension) * 1e-6
    )

    # Initialize starting points
    x_current = np.random.multivariate_normal(
        mean1, np.diag(np.ones(dimension))
    )
    y_current = np.random.multivariate_normal(
        mean2, np.diag(np.ones(dimension))
    )

    for i in range(1, num_iterations + 1):
        # Hit-and-Run sampling for X
        direction_x = np.random.randn(len(mean1))
        lb_x, ub_x = get_non_zero_density_bounds(
            cov1, x_current, direction_x
        )
        x_current = np.random.uniform(lb_x, ub_x)

        # Hit-and-Run sampling for Y
        direction_y = np.random.randn(len(mean2))
        lb_y, ub_y = get_non_zero_density_bounds(
            cov2, y_current, direction_y
        )
        y_current = np.random.uniform(lb_y, ub_y)

        # Sequential update on cross-covariance matrix
        deviation_x = x_current - mean1
        deviation_y = y_current - mean2
        cross_covariance = (i / (i + 1)) * cross_covariance + (
            i / (i + 1) ** 2
        ) * np.outer(deviation_x, deviation_y)

        # Update rmse values and frobenius_norms
        rmse = np.sqrt(np.mean((cross_covariance - ground_truth) ** 2))
        rmse_values.append(rmse)
        fro_norm = np.linalg.norm(cross_covariance - ground_truth, "fro")
        frobenius_norms.append(fro_norm)

    return cross_covariance, rmse_values, frobenius_norms


def main(
    dimension: int,
    num_iterations: int = 1000,
    correlation_factor: float = 0,
) -> list:
    """
    Process and estimate cross-covariance for a given dimension with an
        imposed correlation factor.

    Args:
        dimension (int): The dimensionality of the distribution.
        num_iterations (int, optional): Number of iterations for sampling.
            Default is 1000.
        correlation_factor (float, optional): Factor to impose correlation
            between the two distributions. Default is 0.

    Returns:
        list: List of RMSE values over iterations and list of Frobenius norms
            for iterations.
    """
    mean1, cov1 = generate_random_params(dimension)
    mean2, cov2 = generate_random_params(dimension)

    # Generate correlated samples
    _, _, correlated_samples = generate_non_independent_data(
        mean1, cov1, mean2, cov2, num_iterations, correlation_factor
    )

    # Calculate the ground truth cross-covariance matrix
    ground_truth = np.cov(correlated_samples, rowvar=False)

    _, rmse, frobenius_norms = estimate_cross_covariance(
        mean1, cov1, mean2, cov2, ground_truth, num_iterations
    )
    return rmse, frobenius_norms


if __name__ == "__main__":
    dimensions = [100]
    num_iterations = 1000
    correlation_factors = [0.0, 0.1, 0.2, 0.5]

    # Dictionary to store results for each correlation factor
    all_results = {}

    # Compute results for each correlation factor and dimension
    for factor in correlation_factors:
        all_results[factor] = {
            dim: main(dim, num_iterations, factor) for dim in dimensions
        }

    plt.figure(figsize=(16, 5))
    lines = []  # To store line objects for the legend
    labels = []  # To store label strings for the legend

    for i, factor in enumerate(correlation_factors):
        plt.subplot(2, 4, i + 1)

        # Retrieve results for this correlation factor
        results = all_results[factor]

        # Plotting RMSE values
        for dim, (rmse, _) in results.items():
            (line,) = plt.plot(rmse, label=f"d={dim}", linewidth=1)
            if i == 0:  # Only add to legend for the first subplot
                lines.append(line)
                labels.append(f"d={dim}")

        plt.title(f"$\\rho$: {factor}", fontsize=12)
        plt.xlabel("# of iterations", fontsize=10)
        plt.ylabel("RMSE", fontsize=10)
        plt.xlim(0, num_iterations)
        plt.grid(True)

    for i, factor in enumerate(correlation_factors):
        plt.subplot(2, 4, i + 5)

        # Retrieve results for this correlation factor
        results = all_results[factor]

        # Plotting Frobenius norms
        for dim, (_, fro_norm) in results.items():
            plt.plot(fro_norm, label=f"d={dim}", linewidth=1)

        plt.title(f"$\\rho$: {factor}", fontsize=12)
        plt.xlabel("# of iterations", fontsize=10)
        plt.ylabel("Frobenius norm", fontsize=10)
        plt.xlim(0, num_iterations)
        plt.grid(True)

    plt.figlegend(
        lines,
        labels,
        loc="upper center",
        ncol=len(dimensions),
        fontsize=10,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig("../../figures/hit_and_run_mcmc.pdf", format='pdf', dpi=300)
    plt.show()

