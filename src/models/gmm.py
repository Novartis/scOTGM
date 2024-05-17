# -------------------------------------------------------
# Copyright 2023 Novartis Institutes for BioMedical Research, INC
# Author: Andac Demir
# Licensed under the GNUv3 License (the "License");
# you may not use this file except in compliance with the License.
# -------------------------------------------------------
import io
import random
import warnings
import logging
from typing import Optional, Tuple
import imageio
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import LabelEncoder
from scipy.stats import norm
from scipy.linalg import sqrtm
from sklearn.svm import SVC

# This will ignore all FutureWarnings
warnings.simplefilter(action="ignore", category=FutureWarning)
# Ignore all DeprecationWarnings
warnings.simplefilter(action="ignore", category=DeprecationWarning)

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def frobenius_norm_distance(mean1, cov1, mean2, cov2, alpha=1e-4):
    """
    Compute a combined measure of distance between two multivariate Gaussian distributions, using
    the Frobenius norm of the difference between their means and covariance matrices.

    The Frobenius norm gives a measure of the difference in the shape of the distributions, and the
    difference in their barycenters.

    Args:
        mean1 (np.array): Mean of the first Gaussian distribution.
        cov1 (np.array): Covariance matrix of the first Gaussian distribution.
        mean2 (np.array): Mean of the second Gaussian distribution.
        cov2 (np.array): Covariance matrix of the second Gaussian distribution.

    Returns:
        float: The combined distance measure between the two distributions.
    """
    # Compute the Euclidean distance between the means
    mean_distance = np.linalg.norm(mean1 - mean2)
    
    # Compute the Frobenius norm of the difference between covariance matrices
    cov_distance = np.linalg.norm(cov1 - cov2, 'fro')
    
    # Combine the distances (simple addition or other methods could be used depending on the application)
    combined_distance = mean_distance + alpha * cov_distance
    
    return combined_distance


class GMMTrainer:
    def __init__(self, X: np.ndarray, cell_types: np.ndarray, donor_ids: np.ndarray, n_control: int = 1, n_pd: int = 1,
                 variance_threshold: float = 0.1, random_seed: Optional[int] = None) -> None:
        """
        Initialize the GMMTrainer and split data into training and testing based on the leave_out_ids.

        Args:
            X (np.array): Gene expression matrix of shape (number of cells x number of genes).
            cell_types (np.array): Cell type classification matrix of shape (number of cells x number of cell types).
            donor_ids (np.array): Patient id matrix of shape (number of cells x number of subjects).
            n_control (int, optional): Number of control samples to be used. Defaults to 1.
            n_pd (int, optional): Number of PD samples to be used. Defaults to 1.
            variance_threshold (float): minimum captured variance by the principal components

        Note:
            Use variance_threshold=0.1 for Smajic and Feleke Datasets.
            Use variance_threshold=0.4 for Kamath Dataset.
        """
        if random_seed is not None:
            random.seed(random_seed)
                            
        control = [i for i in donor_ids if i.startswith("C") or i.startswith("Control")]
        pd = [i for i in donor_ids if i.startswith("PD") or i.startswith("Parkinson's")]

        # self.leave_out_ids (List[str]): Arbitrarily selected donor IDs to be left out for testing, e.g., ['C6', 'PD5']
        self.leave_out_ids = random.sample(sorted(set(control)), n_control) + \
                             random.sample(sorted(set(pd)), n_pd)
        logging.info(f"Control IDs: {sorted(set(control))}")
        logging.info(f"PD IDs: {sorted(set(pd))}")
        logging.info(f"Leave out IDs: {self.leave_out_ids}")

        # Modify donor_ids
        donor_ids_modified = np.array(["C" if "C" in donor or "Control" in donor else "PD" for donor in donor_ids])

        # Create cross-product label y
        y = [f"{cell}_{donor}" for cell, donor in zip(cell_types, donor_ids_modified)]
        
        # Instantiate a LabelEncoder
        self.le = LabelEncoder()
        y_encoded = self.le.fit_transform(y)

        # Log the mapping from original labels to encoded labels
        label_mapping = {original: encoded for original, encoded in zip(y, y_encoded)}
        sorted_label_mapping = sorted(label_mapping.items(), key=lambda item: item[1])
        for original, encoded in sorted_label_mapping:
            logging.info(f"Original Label: '{original}' -> Encoded Label: {encoded}")

        # Get boolean masks for train and test splits
        train_mask = ~np.isin(donor_ids, self.leave_out_ids)
        test_mask = np.isin(donor_ids, self.leave_out_ids)
        
        # Use the boolean masks to split data
        self.split_data(X, y_encoded, train_mask, test_mask)

        # Apply PCA to the training data and then transform test data using the same components
        self.variance_threshold = variance_threshold
        self.apply_pca()

    def print_encoded_labels(self) -> None:
        """
        Print the original labels and their corresponding integer labels.
        """
        for original_label, encoded_label in zip(self.le.classes_, range(len(self.le.classes_))):
            logging.info(f"Original Label: {original_label} -> Encoded Label: {encoded_label}")
        logging.info(50 * "#")

    def split_data(self, X: np.ndarray, y_encoded: np.ndarray, train_mask: np.ndarray, test_mask: np.ndarray) -> None:
        """Use boolean masks to split X and y into train and test subsets."""
        self.X_train = X[train_mask]
        self.y_train = np.array(y_encoded)[train_mask]

        self.X_test = X[test_mask]
        self.y_test = np.array(y_encoded)[test_mask]
        
        unique_elements, counts = np.unique(self.y_train, return_counts=True)
        train_result = dict(zip(unique_elements, counts))
        logging.info(f"Train labels distribution: {train_result}")

        unique_elements, counts = np.unique(self.y_test, return_counts=True)
        test_result = dict(zip(unique_elements, counts))
        logging.info(f"Test labels distribution: {test_result}")

    def apply_pca(self) -> None:
        """
        Fit PCA to training data and transform both training and test data.
        """
        pca_full = PCA(svd_solver='full')
        pca_full.fit(self.X_train)

        # Determine number of components to explain the desired variance
        explained_variances = pca_full.explained_variance_ratio_
        cumulative_variances = np.cumsum(explained_variances)
        n_components = np.where(cumulative_variances >= self.variance_threshold)[0][0] + 1
        logging.info(f'Number of principal components to keep after PCA: {n_components}')

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(cumulative_variances, marker="o", linestyle="--")
        plt.title("Variance Explained vs. Number of Components")
        plt.xlabel("Number of Components")
        plt.ylabel("Cumulative Variance Explained")
        plt.axvline(x=n_components, color="r", linestyle="--", label=f"{n_components} components")
        plt.axhline(y=self.variance_threshold, color="g", linestyle="--", label=f"{self.variance_threshold * 100}% variance")
        plt.legend()
        plt.savefig("figures/pca_component_variance_plot.png", format='png', dpi=300)
        plt.show()

        # Apply PCA with desired number of components
        self.pca = PCA(n_components=n_components, svd_solver='full')

        self.X_train = self.pca.fit_transform(self.X_train)
        self.X_test = self.pca.transform(self.X_test)

    def decode_labels(self, encoded_labels: np.ndarray) -> np.ndarray:
        """Revert integer labels back to original string labels."""
        return self.le.inverse_transform(encoded_labels)

    def train(self, verbose: bool = False) -> None:
        """Train the GMM model."""
        self.n_clusters = len(np.unique(self.y_train))
                
        # Calculate initial means
        self.initial_means = np.array([self.X_train[self.y_train == i].mean(axis=0) for i in range(self.n_clusters)])
        # Calculate initial covariances. np.cov treats rows as features and columns as samples, transpose data matrix.
        self.initial_covs = np.array(
            [np.cov(self.X_train[self.y_train == i].T) for i in range(self.n_clusters)]
        )
        # Covariance regularization to ensure invertibility and stability
        self.epsilon = 1e-3
        for i in range(self.n_clusters):
            self.initial_covs[i] += np.eye(self.X_train.shape[1]) * self.epsilon
        # Calculate initial priors using the frequency of cell phenotypes and their pathology
        self.initial_weights = np.array([(self.y_train == i).mean() for i in range(self.n_clusters)])

        self.gmm = GaussianMixture(
            n_components=self.n_clusters,
            means_init=self.initial_means,
            precisions_init=np.linalg.inv(self.initial_covs),
            weights_init=self.initial_weights,
            max_iter=1,
        )
        self.gmm.fit(self.X_train)
        logging.info(f"Number of iterations until convergence: {self.gmm.n_iter_}")

        # Print means and covariances
        if verbose:
            for i in range(self.gmm.n_components):
                logging.info(f"Component {i + 1}:")
                logging.info("Component means:")
                logging.info(self.gmm.means_[i])
                logging.info("Determinant of the covariance matrix to interpret its volume:")
                logging.info(np.linalg.det(self.gmm.covariances_[i]))
                logging.info(50 * "#")

    def predict(self) -> np.ndarray:
        """
        Predict using the trained GMM model and match the clusters to labels between the initial and final GMM components.
        """
        if self.gmm is None:
            raise ValueError("GMM has not been trained. Call `train` first.")

        # Predict the cluster assignments
        gmm_labels = self.gmm.predict(self.X_test)

        matching = {}
        for i, (final_mean, final_cov) in enumerate(zip(self.gmm.means_, self.gmm.covariances_)):
            costs = [frobenius_norm_distance(self.initial_means[j], self.initial_covs[j], final_mean, final_cov)
                     for j in range(self.n_clusters)]
            closest_label = np.argmin(costs)
            matching[i] = closest_label

        # Map the predicted cluster assignments to the corresponding labels
        y_pred = np.array([matching[label] for label in gmm_labels])

        return y_pred

    def evaluate(self, y_pred: np.ndarray) -> float:
        """Evaluate the GMM model."""
        accuracy = accuracy_score(self.y_test, y_pred)
        f1 = f1_score(self.y_test, y_pred, average="weighted")

        # Decode labels
        decoded_y_test = self.le.inverse_transform(self.y_test)
        decoded_y_pred = self.le.inverse_transform(y_pred)

        # Confusion matrix with decoded labels
        cm = confusion_matrix(decoded_y_test, decoded_y_pred, labels=self.le.classes_)
        plt.figure(figsize=(12, 10))
        ax = sns.heatmap(
            cm, annot=True, fmt="d", cmap="Blues", xticklabels=self.le.classes_, yticklabels=self.le.classes_
        )

        # Rotate x-axis labels and adjust font size
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=6)

        # Rotate y-axis labels and adjust font size
        ax.set_yticklabels(ax.get_yticklabels(), rotation=45, fontsize=6)

        plt.ylabel("Actual", fontsize=8)
        plt.xlabel("Predicted", fontsize=8)
        plt.title(f"Confusion Matrix on Test Set (GMM)- Accuracy: {accuracy:.4f}, F1: {f1: .4f}")
        plt.tight_layout()  # Adjust the layout to accommodate rotated labels
        plt.savefig("figures/gmm_confusion_matrix.png", format='png', dpi=300)
        plt.show()
        
        # Logging the results
        logging.info(f"Before in silico perturbation")
        logging.info(f'Accuracy: {accuracy}')
        logging.info(f'F1 Score (Weighted): {f1}')
        logging.info('Confusion Matrix:')
        logging.info('\n' + str(cm))

        return f1


class Visualizer:
    @staticmethod
    def plot_gmm_pca(X_test: np.ndarray, y_test: np.ndarray, label_encoder: LabelEncoder) -> None:
        """Visualize GMM with PCA reduction."""
        pca_plot = PCA(n_components=4)
        X_pca = pca_plot.fit_transform(X_test)

        fig, axes = plt.subplots(3, 3, figsize=(15, 15))

        # Upper triangle combinations of PCAs
        combinations = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

        for index, (x, y) in enumerate(combinations):
            ax = axes[x, y - 1 if y != 0 else y]
            Visualizer._scatter_on_axes(X_pca, y_test, x, y, ax, label_encoder, show_legend=(index == 0))

        # Create the legend in the first subplot
        handles, labels = axes[0, 0].get_legend_handles_labels()
        # Add the legend to the empty subplot in the second row, first column
        axes[1, 0].legend(handles, labels, loc="center", fontsize="small")
        # Hide the axis for the subplot
        axes[1, 0].axis("off")

        # Hide the rest of the lower triangle subplots
        for i in range(3):
            for j in range(3):
                if i > j:
                    axes[i, j].axis("off")

        plt.tight_layout()
        plt.show()

    @staticmethod
    def _scatter_on_axes(X_pca: np.ndarray, y_test: np.ndarray, x: int, y: int, ax: plt.Axes,
                         label_encoder: LabelEncoder, show_legend: bool = True) -> None:
        unique_labels = np.unique(y_test)
        for label in unique_labels:
            mask = y_test == label
            decoded_label = label_encoder.inverse_transform([label])[0]
            ax.scatter(X_pca[mask, x], X_pca[mask, y], label=f"Cluster {decoded_label}", s=5, alpha=0.5)

        ax.set_xlabel(f"PCA {x + 1}")
        ax.set_ylabel(f"PCA {y + 1}")

        if show_legend:
            # Only collect the legend details but don't display it
            ax.get_legend_handles_labels()
            ax.legend().set_visible(False)


class PerturbationModel:
    def __init__(self, gmm: GaussianMixture, initial_means: np.ndarray, initial_covs: np.ndarray, pca: PCA,
                 le: LabelEncoder, unperturbed_cluster: int, perturbed_cluster: int, num_samples: int = 10000) -> None:
        """
        Initialize the PerturbationModel class and match clusters based on the proximity
        of the initial and final GMM mean vectors.

        Args:
            gmm (GaussianMixture): The fitted GMM model for training data.
            initial_means (np.array): Array containing the initial mean vectors for each cluster.
                Used for matching clusters to true labels.
            initial_covs (np.array): Array containing the initial covariance matrices for each cluster.
                Used for matching clusters to true labels.
            pca (PCA): The fitted PCA model to the training data.
            le (LabelEncoder): Encodes target labels into values between 0 and n_classes-1.
            unperturbed_cluster (int): Index of the starting cluster from which the perturbation originates.
            perturbed_cluster (int): Index of the target cluster where the perturbation is applied.
        """
        self.gmm = gmm
        self.pca = pca
        self.le = le

        self.unperturbed_cluster = unperturbed_cluster
        self.perturbed_cluster = perturbed_cluster

        self.unperturbed_cluster_matched = self.match_cluster_to_label(unperturbed_cluster, initial_means, initial_covs)
        self.perturbed_cluster_matched = self.match_cluster_to_label(perturbed_cluster, initial_means, initial_covs)
    
    def match_cluster_to_label(self, cluster_id: int, initial_means: np.ndarray,
                               initial_covs: np.ndarray) -> int:
        """
        This function aims to match a specific cluster from a trained GMM to an "initial" cluster label 
        based on the initial means and covariances provided. 

        Args:
            cluster_id (int): The index of the cluster in the trained GMM that you want to match with an initial label.
            initial_means (np.array): Array containing the mean vectors of the initial clusters, before the GMM was trained.
            initial_covs (np.array): Array containing the covariance matrices of the initial clusters.

        Returns:
            int: The index of the matched initial cluster. This index serves as the "label" that best matches 
                the specified cluster from the trained GMM, effectively resolving the permutation ambiguity by 
                providing a consistent way to refer to clusters across different instances of the model.
        """
        initial_mean = initial_means[cluster_id]
        initial_cov = initial_covs[cluster_id]

        costs = []  

        for i in range(len(self.gmm.means_)):
            gmm_mean = self.gmm.means_[i]
            gmm_covariance = self.gmm.covariances_[i]

            cost = frobenius_norm_distance(initial_mean, initial_cov, gmm_mean, gmm_covariance)
            costs.append(cost)

        matched_cluster = np.argmin(costs)
        return matched_cluster

    def get_cross_covariance_by_monte_carlo_sampling(self, num_iterations=10000, visualize: bool = False) -> np.ndarray:
        """
        Estimate the cross-covariance matrix using Monte Carlo sampling.
        This technique involves direct sampling from the 2 Gaussian components and estimating the cross-covariance
        matrix using the outer product of sample deviations from their means. -> Simple and straightforward approach,
        but requires that the 2 distributions have strong dependencies.
        
        Args:
            num_iterations (int): Number of samples to draw from GMM.

        Returns:
            numpy.ndarray: Cross-covariance matrix.
        """
        # Sample from the two clusters
        samples, labels = self.gmm.sample(num_iterations)
        samples1 = samples[labels == self.unperturbed_cluster_matched]
        samples2 = samples[labels == self.perturbed_cluster_matched]

        # Make sure samples1 and samples2 have the same length
        min_len = min(len(samples1), len(samples2))
        samples1 = samples1[:min_len]
        samples2 = samples2[:min_len]

        # Compute means of the samples
        mean1 = samples1.mean(axis=0)
        mean2 = samples2.mean(axis=0)

        if visualize is True:
            # For visualization purposes
            intermediate_covs = []

            for i in range(1, num_iterations + 1):
                cross_cov_i = (samples1[:i] - mean1).T @ (samples2[:i] - mean2) / i
                intermediate_covs.append(cross_cov_i)

            # Save the intermediate covariances for visualization
            self.intermediate_covs = intermediate_covs
            # Visualize the change in the cross-covariance matrix estimation
            self.visualize_covariance_evolution()
            return intermediate_covs[-1]

        else:
            return (samples1 - mean1).T @ (samples2 - mean2) / num_iterations

    def get_cross_covariance_by_hit_and_run_gibbs_sampling(self, num_iterations: int = 1000) -> np.ndarray:
        """
        Estimate the cross-covariance matrix using Hit-and-Run Gibbs Sampling.

        Args:
            num_iterations (int): Number of iterations to perform.

        Returns:
            np.ndarray: Estimated cross-covariance matrix.
        """
        # Get means and covariances for each cluster
        means_x = self.gmm.means_[self.unperturbed_cluster_matched]
        cov_x = self.gmm.covariances_[self.unperturbed_cluster_matched]
        means_y = self.gmm.means_[self.perturbed_cluster_matched]
        cov_y = self.gmm.covariances_[self.perturbed_cluster_matched]

        # Initialize starting points
        x_current = np.random.multivariate_normal(means_x, cov_x)
        y_current = np.random.multivariate_normal(means_y, cov_y)

        # Initialize cross_covariance
        cross_covariance = np.zeros_like(cov_x)
        np.fill_diagonal(cross_covariance, 1e-4)

        for i in range(1, num_iterations + 1):
            # Hit-and-Run sampling for X
            direction_x = np.random.randn(len(means_x))
            lb_x, ub_x = PerturbationModel.get_non_zero_density_bounds(cov_x, x_current, direction_x)
            x_current = np.random.uniform(lb_x, ub_x)

            # Hit-and-Run sampling for Y
            direction_y = np.random.randn(len(means_y))
            lb_y, ub_y = PerturbationModel.get_non_zero_density_bounds(cov_y, y_current, direction_y)
            y_current = np.random.uniform(lb_y, ub_y)

            # Sequential update on cross-covariance matrix
            deviation_x = x_current - means_x
            deviation_y = y_current - means_y
            cross_covariance = (i / (i + 1)) * cross_covariance + (
                    i / (i + 1) ** 2
            ) * np.outer(deviation_x, deviation_y)

        return cross_covariance

    @staticmethod
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

    def visualize_covariance_evolution(self) -> str:
        """
        Visualize the change in the covariance matrix estimation with each Monte Carlo sample.

        Returns:
            str: Path to the generated GIF.
        """
        images = []

        for cov in self.intermediate_covs:
            fig, ax = plt.subplots()
            sns.heatmap(cov, ax=ax, cmap="Blues")

            # Convert the figure to an image
            buf = io.BytesIO()
            plt.savefig(buf, format="png")
            buf.seek(0)
            images.append(imageio.imread(buf))

            plt.close(fig)

        # Create GIF from the list of images
        gif_path = "../eval/covariance_evolution.gif"
        imageio.mimsave(gif_path, images, duration=0.0005, loop=1)

        return gif_path

    @staticmethod
    def make_positive_semidefinite(matrix):
        """
        Converts a given square matrix into a positive semidefinite matrix.

        This method computes the eigenvalues and eigenvectors of the input matrix using
        the Hermitian eigenvalue decomposition. It then adjusts any negative eigenvalues to 
        a very small positive value (1e-99) to ensure that the resulting matrix is positive 
        semidefinite. The adjusted eigenvalues are used to reconstruct a positive semidefinite 
        matrix using the original eigenvectors.

        Args:
            matrix (numpy.ndarray): A square numpy array representing the matrix to be converted
                into a positive semidefinite matrix.
        
        Returns:
            np.ndarray: A positive semidefinite matrix derived from the input matrix.
        """
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        eigenvalues[eigenvalues < 0] = 1e-99  # Set any small negative eigenvalues to a small positive value
        return eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T
    
    def additive_parameters(self, method: str = "hit and run mcmc") -> Tuple[np.ndarray, np.ndarray]:
        """
        Find the parameters of the additive Gaussian noise between two clusters.

        Args:
            method (str): Method to estimate the cross-covariance matrix when the joint distribution is unknown.
                Available options: "monte carlo" and "hit and run mcmc". Defaults to "hit and run mcmc".

        Returns:
            tuple: Mean and covariance matrix of the additive Gaussian noise.
        """
        # Mean of the additive Gaussian noise
        perturb_mean = self.gmm.means_[self.perturbed_cluster_matched] - self.gmm.means_[self.unperturbed_cluster_matched]

        # Covariance of the additive Gaussian noise
        if method == "monte carlo":
            cross_cov = self.get_cross_covariance_by_monte_carlo_sampling(num_iterations=10000)
        elif method == "hit and run mcmc":
            np.set_printoptions(threshold=np.inf)
            cross_cov = self.get_cross_covariance_by_hit_and_run_gibbs_sampling(num_iterations=10000)
        else:
            raise ValueError(
                f"Method '{method}' is not supported. Choose 'monte carlo' or 'hit and run mcmc'.")

        perturb_cov = (
            self.gmm.covariances_[self.unperturbed_cluster_matched]
            + self.gmm.covariances_[self.perturbed_cluster_matched]
            - 2 * cross_cov
        )
        perturb_cov = PerturbationModel.make_positive_semidefinite(perturb_cov)

        return perturb_mean, perturb_cov
        
    def sample_from_perturbation(self, perturb_mean: np.ndarray, perturb_cov: np.ndarray, n_samples: int = 10000) -> Tuple[np.ndarray, np.ndarray]:
        """Sample from the perturbation distribution."""
        # Extract mean and covariance of the perturbation distribution
        samples = np.random.multivariate_normal(perturb_mean, perturb_cov, n_samples)
        samples = self.pca.inverse_transform(samples)

        # Return mean and standard deviation for each gene
        gene_perturbation_means = np.mean(samples, axis=0)  # magnitude of up/down regulation in perturbation
        gene_perturbation_std_devs = np.std(samples, axis=0)  # uncertainty (confidence) in estimates

        return gene_perturbation_means, gene_perturbation_std_devs

    def test_perturbation_with_gmm(self, X_test: np.ndarray, y_test: np.ndarray, perturb_mean: np.ndarray, perturb_cov: np.ndarray) -> None:
        """
        Adds a perturbation to samples belonging to a specific cluster, identified by self.unperturbed_cluster,
        using a Gaussian distribution defined by perturb_mean and perturb_cov.
        These adversarial examples are designed to be misclassified by a pre-trained GMM model, 
        aiming to alter the original samples subtly enough to cause misclassification, 
        while remaining as indistinguishable as possible from the original samples.
        
        Args:
            X_test: np.ndarray, feature matrix for test samples.
            y_labels: np.ndarray, labels for test samples.
            perturb_mean: np.ndarray, mean vector of the Gaussian perturbation.
            perturb_cov: np.ndarray, covariance matrix of the Gaussian perturbation.
        """
        # Extract samples belonging to the specified cluster
        cluster_samples = X_test[y_test == self.unperturbed_cluster]

        # Sample vectors from perturbation distribution for each sample and apply it
        n_samples = len(cluster_samples)
        noise_vectors = np.random.multivariate_normal(perturb_mean, perturb_cov, n_samples)
        perturbed_data = cluster_samples + noise_vectors
        # perturbed_data = cluster_samples + perturb_mean

        # Predict cluster memberships of the perturbed samples using GMM
        y_pred = self.gmm.predict(perturbed_data)

        decoded_y_test = self.le.inverse_transform(y_test[y_test == self.unperturbed_cluster])
        decoded_y_pred = self.le.inverse_transform(y_pred)

        # Show the confusion matrix and accuracy
        accuracy = accuracy_score(decoded_y_test, decoded_y_pred)
        f1 = f1_score(decoded_y_test, decoded_y_pred, average="weighted")
        cm = confusion_matrix(decoded_y_test, decoded_y_pred, labels=self.le.classes_)
        
        # Logging the results
        logging.info(f"After in silico perturbation (GMM evaluation)")
        logging.info(f'Accuracy: {accuracy}')
        logging.info(f'F1 Score (Weighted): {f1}')
        logging.info('Confusion Matrix:')
        logging.info('\n' + str(cm))

        plt.figure(figsize=(12, 10))
        ax = sns.heatmap(
            cm, annot=True, fmt="d", cmap="Blues", xticklabels=self.le.classes_, yticklabels=self.le.classes_
        )  # Decoded class names for x and y axis

        # Rotate x-axis labels and adjust font size
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=6)

        # Rotate y-axis labels and adjust font size
        ax.set_yticklabels(ax.get_yticklabels(), rotation=45, fontsize=6)

        plt.ylabel("Actual", fontsize=8)
        plt.xlabel("Predicted", fontsize=8)
        plt.title(f"Confusion Matrix on Test Set after Perturbation (GMM)- Accuracy: {accuracy:.4f}, F1: {f1:.4f}")
        plt.savefig("figures/gmm_confusion_matrix_after_perturbation.png", format='png', dpi=300)
        plt.show()

    def test_perturbation_with_svm(self, X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray,
                                   y_test: np.ndarray, perturb_mean: np.ndarray, perturb_cov: np.ndarray) -> None:
        """
        Adds a perturbation to samples belonging to a specific cluster, identified by self.unperturbed_cluster,
        using a Gaussian distribution defined by perturb_mean and perturb_cov.
        These adversarial examples are designed to be misclassified by a pre-trained SVM model, 
        aiming to alter the original samples subtly enough to cause misclassification, 
        while remaining as indistinguishable as possible from the original samples.
        
        Args:
            X_test: np.ndarray, feature matrix for test samples.
            y_labels: np.ndarray, labels for test samples.
            perturb_mean: np.ndarray, mean vector of the Gaussian perturbation.
            perturb_cov: np.ndarray, covariance matrix of the Gaussian perturbation.
        """
        svm = SVC(kernel="poly", degree=3, C=1)  # Support Vector Classifier with a degree 3 polynomial kernel
        svm.fit(X_train, y_train)

        # Show the confusion matrix and accuracy
        y_pred = svm.predict(X_test)

        # Assuming you have an instance of the label encoder called label_encoder
        decoded_y_test = self.le.inverse_transform(y_test)
        decoded_y_pred = self.le.inverse_transform(y_pred)

        # Show the confusion matrix and accuracy
        accuracy = accuracy_score(decoded_y_test, decoded_y_pred)
        f1 = f1_score(decoded_y_test, decoded_y_pred, average="weighted")
        cm = confusion_matrix(decoded_y_test, decoded_y_pred, labels=self.le.classes_)
        
        logging.info(f"After in silico perturbation (SVM evaluation)")
        logging.info(f'Accuracy: {accuracy}')
        logging.info(f'F1 Score (Weighted): {f1}')
        logging.info('Confusion Matrix:')
        logging.info('\n' + str(cm))

        plt.figure(figsize=(12, 10))
        ax = sns.heatmap(
            cm, annot=True, fmt="d", cmap="Blues", xticklabels=self.le.classes_, yticklabels=self.le.classes_
        )  # Decoded class names for x and y axis

        # Rotate x-axis labels and adjust font size
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=6)

        # Rotate y-axis labels and adjust font size
        ax.set_yticklabels(ax.get_yticklabels(), rotation=45, fontsize=6)

        plt.ylabel("Actual", fontsize=8)
        plt.xlabel("Predicted", fontsize=8)
        plt.title(
            f"Confusion Matrix on Test Set (SVC with Polynomial Kernel of degree 3)- Accuracy: {accuracy:.4f}, "
            f"F1: {f1: .4f}"
        )
        plt.savefig("figures/svm_confusion_matrix.png", format='png', dpi=300)
        plt.show()

        # Identify the samples from the unperturbed_cluster
        cluster_samples = self.pca.inverse_transform(X_test[y_test == self.unperturbed_cluster])

        # Add the perturbation mean to those samples
        perturbed_data = cluster_samples + perturb_mean

        # Predict using the SVM
        y_pred = svm.predict(self.pca.transform(perturbed_data))

        # Assuming you have an instance of the label encoder called label_encoder
        decoded_y_test = self.le.inverse_transform(y_test[y_test == self.unperturbed_cluster])
        decoded_y_pred = self.le.inverse_transform(y_pred)

        # Compute the accuracy and confusion matrix with decoded labels
        accuracy = accuracy_score(decoded_y_test, decoded_y_pred)
        f1 = f1_score(decoded_y_test, decoded_y_pred, average="weighted")
        cm = confusion_matrix(decoded_y_test, decoded_y_pred, labels=self.le.classes_)

        plt.figure(figsize=(12, 10))
        ax = sns.heatmap(
            cm, annot=True, fmt="d", cmap="Blues", xticklabels=self.le.classes_, yticklabels=self.le.classes_
        )  # Decoded class names for x and y axis

        # Rotate x-axis labels and adjust font size
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=6)

        # Rotate y-axis labels and adjust font size
        ax.set_yticklabels(ax.get_yticklabels(), rotation=45, fontsize=6)

        plt.ylabel("Actual", fontsize=8)
        plt.xlabel("Predicted", fontsize=8)
        plt.title(
            f"Confusion Matrix on Test Set after Perturbation (SVC with polynomial Kernel)- "
            f"Accuracy: {accuracy:.4f}, F1: {f1:.4f}"
        )
        plt.savefig("figures/svm_confusion_matrix_after_perturbation.png", format='png', dpi=300)
        plt.show()