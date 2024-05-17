# sc-OTGM: Single-Cell Perturbation Modeling by Solving Optimal Mass Transport on the Manifold of Gaussian Mixtures

![Image](figures/hueless_a_stunning_coding_competition_poster_that_shows_DNA_PR_bd408e66-aeee-4a06-ae63-9fd87c705acf.png)

In single-cell genomics, the integration of machine learning techniques, particularly foundation models inspired by advancements in large language models, has significantly enhanced our ability to interpret complex biological data. However, the effectiveness of these large networks in scenarios characterized by sparse and noisy data is unclear. Single-cell data, unlike its NLP counterparts, often suffers from technical artifacts, dropout events, and batch effects, which complicates data analysis, especially under weakly supervised conditions.

Addressing these challenges, we introduce sc-OTGM, a streamlined model with fewer than 500,000 parameters. This model is roughly 100 times more compact than typical state-of-the-art foundation models, providing an efficient tool for genomic data analysis without the overhead of larger models. sc-OTGM adopts an unsupervised approach, leveraging a Gaussian mixture model (GMM) to form a probabilistic latent space that effectively distinguishes between different cell populations by learning their individual probability density functions. Through the use of a Hit-and-Run Markov chain sampler, sc-OTGM optimizes the transport plan across these densities within the GMM framework.

We have rigorously tested sc-OTGM on the CROP-seq dataset, which includes 57 single-gene perturbations. Our results demonstrate that sc-OTGM can 

1.  Perform cell phenotype annotation,
2.  Facilitate the analysis of differential gene expression, 
3.  Rank genes top to bottom for target identification (recommender system),
4.  Predict how single-gene perturbations affect the regulation of downstream genes, and
5.  Generate synthetic sc-RNA data conditioned on a specific subject, cellular phenotype or disease.

For more detailed information about the methodologies and experiments conducted in our study, please refer to our paper: [sc-OTGM: Single-Cell Perturbation Modeling by Solving Optimal Mass Transport on the Manifold of Gaussian Mixtures](https://openreview.net/pdf?id=qjFwY7cwXe)


## Setup and Installation
* It is highly suggested that you install all dependencies into a separate conda virtual environment for easy package management.
* The dependencies are in [`requirements.txt`](requirements.txt). You will need to install dependencies by running in the root directory:
    ```shell
    $conda create -n <myenv> --file requirements.txt python=3.10
    $conda activate <myenv>  
    $pip install -e .
    $python -m ipykernel install --user --name=<myenv>
    ```
* Please check you have the same versions of these dependencies.


## License
The code is licensed under GPLv3. For more details, see the [LICENSE](LICENSE) file.


## Getting benchmark dataset

#### CROP-seq
To evaluate the performance of sc-OTGM, we used the CROP-seq dataset from in-vitro experiments on human-induced pluripotent stem cell (iPSC)-derived neurons subjected to genetic perturbations [Tian2021]. Raw published data is available from the Gene Expression Omnibus under accession code [GSE152988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152988). These perturbations were implemented using CRISPRi technology, which enables the selective knockdown of specific genes, reflecting a targeted approach to studying gene function and its impact on neuronal survival and oxidative stress responses among other processes. Using the [`rank_genes_groups`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) method from the scanpy package [Wolf2018] for differential expression analysis, we scrutinized the effects of knocking down 185 genes identified as potentially relevant to neuronal health and disease states. Of these, only 57 genes met our significance threshold (adjusted p-value less than 0.05), indicating a significant alteration in expression levels post-perturbation.


## Code style
Perform these steps manually in the root directory:

```shell
source format_and_lint.sh
format-and-lint .
```


### References
[Tian2021]: Tian, R., Abarientos, A., Hong, J. et al. "Genome-Wide CRISPRi/a Screens in Human Neurons Link Lysosomal Failure to Ferroptosis". Nature Neuroscience, 24(7), 1020–1034, 2023, Nature Publishing Group.

[Wolf2018]: Wolf, F Alexander and Angerer, Philipp and Theis, Fabian J. "SCANPY: large-scale single-cell gene expression data analysis". Genome Biology, 19, 1-5, 2018, Springer.

Corresponding author: Andac Demir (andac.demir@novartis.com), AICS, BR 
