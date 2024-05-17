.. GMM for Target ID documentation master file, created by
   sphinx-quickstart on Tue Dec  5 13:04:54 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

sc-OTGM: Single-Cell Perturbation Modeling by Solving Optimal Mass Transport on the Manifold of Gaussian Mixtures
=================================================================================================================

**sc-OTGM effectively performs**:

1. **Cell type & state annotation**: Identifies and classifies cell types within a dataset, enhancing the understanding of cellular diversity and function.

2. **Differential gene expression analysis**: Enables the comparison of gene expression levels across different conditions, cell types, or treatments, providing insights into biological processes and disease mechanisms.

3. **Prediction of single-gene perturbations on downstream gene regulation**: Offers insights into the effects of altering a single gene on the regulation of other genes, aiding in the comprehension of genetic networks and pathways.


**Overview**

Single-cell RNA sequencing (sc-RNAseq) has emerged as a powerful tool to dissect cellular heterogeneity and understand disease pathogenesis at the individual cell level. By profiling the transcriptomes of thousands to millions of cells in a single experiment, it provides a high-resolution lens into the molecular intricacies of various diseases. Given the vast dimensionality of gene expression profiles, understanding the inherent structure and distributions becomes crucial for effective analyses.

The core function of sc-OTGM is to create a probabilistic latent space utilizing a Gaussian mixture model (GMM) as its prior distribution and distinguish between distinct cell populations by learning their respective marginal probability density functions (PDFs). It then employs a Hit-and-Run Markov chain sampler to learn the optimal transport (OT) plan across these PDFs within the GMM manifold.

sc-OTGM is a probabilistic framework simulating the cellular responses to combinatorial gene perturbations in diseased cells. This enables a deeper understanding of the regulatory processes required to shift them towards a healthier state.

One fundamental assumption in our analysis is that gene expression data, particularly from sc-RNAseq, can be approximated by a multivariate Gaussian distribution. Let :math:`X` represent the gene expression matrix, where each row corresponds to a single cell and each column to a gene. The expression level of genes in cell :math:`i` is given by :math:`X_{i}`. The Gaussian assumption can be formally written as:

.. math::
   X_{i} \sim \mathcal{N}(\mu_k, \Sigma_k)

where :math:`\mu_k` and :math:`\Sigma_k` are the mean and covariance of the underlying distribution that generates the sc-RNA data for a particular cell phenotype :math:`k`. This assumption allows us to leverage statistical methods tailored for Gaussian distributions, facilitating robust modeling of the underlying distribution that generates sc-RNA data. Building on this assumption, we propose a computationally efficient, probabilistic framework to infer gene perturbations for drug target identification.

Evaluation on the CROP-seq Dataset
-------------------------------------
To evaluate the performance of sc-OTGM, we used the CROP-seq dataset from in-vitro experiments on human-induced pluripotent stem cell (iPSC)-derived neurons subjected to genetic perturbations. Raw published data is available from the Gene Expression Omnibus under accession code `GSE152988 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152988>`_ [1]. These perturbations were implemented using CRISPRi technology, which enables the selective knockdown of specific genes, reflecting a targeted approach to studying gene function and its impact on neuronal survival and oxidative stress responses among other processes. Using the `rank_genes_groups <https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html>`_ method [2] from the scanpy package for differential expression analysis, we scrutinized the effects of knocking down 185 genes identified as potentially relevant to neuronal health and disease states. Of these, only 57 genes met our significance threshold (adjusted p-value less than 0.05), indicating a significant alteration in expression levels post-perturbation. These findings are summarized in the table below, which includes the genes that presented significant differential expression, highlighting their potential roles in neuronal function and susceptibility to oxidative stress—a key factor in the pathogenesis of neurodegenerative diseases.