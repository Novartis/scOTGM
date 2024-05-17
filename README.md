# In Silico Prediction of Cellular Responses to Gene Perturbations with Distributional Shift Learning under Mixture of Gaussians

## Introduction  
sc-OTGM provides a scalable and unified framework for:
1. Cell state classification
2. Differential gene expression analysis
3. Gene ranking for target identification through a recommender system
4. Perturbation response prediction
5. Generation of synthetic scRNA-seq data 


## Setup and Installation
* It is highly suggested that you install all dependencies into a separate virtual environment for easy package management.
* The dependencies are in [`requirements.txt`](requirements.txt). You will need to install dependencies by running in the root directory:
    ```shell
    $conda create -n <myenv> python=3.10
    $conda activate <myenv>  
    $pip install -e .
    $python -m ipykernel install --user --name=<myenv>
    ```
* Please check you have the same versions of these dependencies.


#### Code style
Perform these steps manually in the root directory:

```shell
source format_and_lint.sh
format-and-lint .
```


## Loading benchmark datasets

#### CROP-seq
To evaluate the performance of sc-OTGM, we used the CROP-seq dataset from in-vitro experiments on human-induced 
pluripotent stem cell (iPSC)-derived neurons subjected to genetic perturbations Tianet al. (2021). 
These perturbations were executed via CRISPRi, enabling targeted gene knockdown to investigate its effects on neuronal 
survival and oxidative stress. Using the rank genes groups method from scanpy package Wolf et al. (2018) 
for differential expression analysis, we scrutinized the effects of knocking down 185 genes identified as potentially 
relevant to neuronal health anddisease states. Of these, only 57 genes met our significance threshold 
(adjusted p-value less than 0.05), indicating a significant alteration in expression levels post-perturbation. 
Raw published data is available from the Gene Expression Omnibus under accession code GSE152988. 


## Baseline comparison on CROP-seq 
In CRISPR interference (CRISPRi) experiments, not all cells receiving the CRISPR components achieve successful targeted
gene knockdown. Additionally, the extent of gene suppression can vary significantly among cells due to variations 
in Cas9 activity, guide RNA efficiency, and individual cellular responses. Following gene knockdown, cells often 
activate compensatory mechanisms that alter the expression profiles of other genes. This requires robust computational 
models that can accurately identify genes that have been knocked down and predict the subsequent changes in the expression of other genes.

To address this, we benchmarked statistical models for their efficacy in differential gene expression analysis 
post-CRISPRi. The techniques include the Mann-Whitney U test, t-test, and sc-OTGM, with the results detailed in Table 1. 
We assessed the ability of each method to rank the true knockeddown gene within the top-k results, with performance 
metrics based on how lower p-values from the Mann-Whitney U test and t-test correlate with higher rankings. 
sc-OTGM demonstrated superior performance, particularly in Top-1 accuracy, showing its effectiveness in identifying 
the most likely perturbed gene. The performance advantage of using sc-OTGM decreases as the ranking threshold increases.

#### Table 1: Benchmarking Statistical Differential Gene Expression Analysis Techniques
| Method             | Top-1 Accuracy | Top-5 Accuracy | Top-10 Accuracy | Top-50 Accuracy | Top-100 Accuracy |
|--------------------|----------------|----------------|-----------------|-----------------|------------------|
| Mann–Whitney U test| 0.37           | 0.40           | 0.42            | 0.58            | 0.60             |
| t-test             | 0.39           | 0.67           | 0.70            | 0.79            | 0.86             |
| sc-OTGM            | **0.56**       | 0.68           | 0.74            | 0.82            | 0.91             |

Table 2 shows the performance of sc-OTGM in identifying differentially expressed genes (DEGs) following targeted gene 
knockdowns in CRISPRi experiments. Based on sc-OTGM’s ranking, a cutoff of 100 is set to select the top predicted DEGs.
Predicted DEGs are obtained from the top of the ranked gene list, and compared against the list of known DEGs. 
We conducted Fisher’s exact tests which yielded p-Values, significantly below the standard threshold of 0.05. 
This confirms a strong statistical correlation between known DEGs after targeted gene knockdown and those predicted by 
sc-OTGM.  

The performance of sc-OTGM is also evaluated based on its accuracy to predict the direction of expression changes 
(upregulation or downregulation) in DEGs. The metrics include Accuracy (%) and F1-score, where the mean accuracy is 
75.2% with a standard deviation of 15.4%, and the mean F1-score is 0.79 with a standard deviation of 0.14. 

#### Table 2: Quantitative Analysis of Gene Perturbation Responses
| Gene    | p-Value        | Accuracy (%) | F1-score |
|---------|----------------|--------------|----------|
| TUBB4A  | 7.37 × 10^−6   | 100.0        | 1.00     |
| ATP1A3  | 3.12 × 10^−20  | 78.7         | 0.85     |
| KIFAP3  | 6.39 × 10^−10  | 87.5         | 0.67     |
| MAPT    | 9.01 × 10^−4   | 100.0        | 1.00     |
| CASP3   | 6.93 × 10^−6   | 100.0        | 1.00     |
| APEX1   | 2.61 × 10^−16  | 100.0        | 1.00     |
| COX10   | 1.43 × 10^−16  | 83.7         | 0.87     |
| NDUFS8  | 2.29 × 10^−36  | 82.7         | 0.83     |
| ZNF292  | 4.91 × 10^−16  | 65.4         | 0.72     |
| GSTA4   | 5.50 × 10^−21  | 89.5         | 0.92     |
| STX1B   | 4.95 × 10^−38  | 72.1         | 0.76     |
| OPTN    | 1.51 × 10^−6   | 100.0        | 1.00     |
| SOD1    | 1.16 × 10^−7   | 88.9         | 0.90     |
| NDUFV1  | 4.33 × 10^−30  | 73.9         | 0.79     |
| CALB1   | 1.38 × 10^−4   | 40.0         | 0.46     |
| EEF2    | 5.08 × 10^−29  | 92.2         | 0.92     |
| BIN1    | 6.39 × 10^−8   | 88.9         | 0.89     |
| SCFD1   | 1.59 × 10^−42  | 56.4         | 0.66     |
| PON2    | 8.41 × 10^−55  | 53.6         | 0.60     |
| BAX     | 1.27 × 10^−30  | 78.3         | 0.83     |
| SCAPER  | 1.48 × 10^−31  | 87.2         | 0.89     |
| CYB561  | 5.66 × 10^−33  | 60.4         | 0.66     |
| AKAP9   | 3.69 × 10^−14  | 100.0        | 1.00     |
| VPS35   | 9.36 × 10^−14  | 80.0         | 0.85     |
| PRNP    | 3.95 × 10^−53  | 59.4         | 0.72     |
| AP2A2   | 2.48 × 10^−57  | 75.4         | 0.82     |
| SOD2    | 2.01 × 10^−7   | 90.5         | 0.91     |
| BECN1   | 6.22 × 10^−4   | 73.1         | 0.75     |
| SNCB    | 6.79 × 10^−41  | 87.7         | 0.89     |
| CDH11   | 1.77 × 10^−14  | 66.7         | 0.70     |
| ELOVL5  | 2.28 × 10^−14  | 92.3         | 0.93     |
| NTRK2   | 4.50 × 10^−29  | 58.7         | 0.61     |
| DAP     | 1.27 × 10^−45  | 81.7         | 0.86     |
| EIF4G1  | 3.00 × 10^−24  | 76.4         | 0.84     |
| TRPM7   | 3.66 × 10^−14  | 66.7         | 0.80     |
| COASY   | 1.51 × 10^−6   | 83.3         | 0.84     |
| TRAP1   | 2.11 × 10^−35  | 80.0         | 0.89     |
| CYP46A1 | 4.56 × 10^−4   | 50.0         | 0.33     |
| PARP1   | 5.51 × 10^−24  | 57.7         | 0.67     |
| FOXRED1 | 1.49 × 10^−25  | 75.8         | 0.78     |
| AFG3L2  | 1.01 × 10^−15  | 75.0         | 0.81     |
| RAB7A   | 1.21 × 10^−12  | 83.3         | 0.84     |
| PPP2R2B | 3.10 × 10^−35  | 63.9         | 0.75     |
| RGS2    | 5.03 × 10^−30  | 63.0         | 0.68     |
| AMFR    | 4.06 × 10^−12  | 62.5         | 0.74     |
| MRPL10  | 1.11 × 10^−22  | 43.6         | 0.57     |
| ANO10   | 1.38 × 10^−4   | 40.0         | 0.46     |
| DMXL1   | 6.72 × 10^−33  | 66.7         | 0.76     |
| HYOU1   | 7.43 × 10^−35  | 55.3         | 0.65     |
| HTT     | 8.81 × 10^−24  | 61.1         | 0.66     |
| ECHS1   | 5.44 × 10^−9   | 71.4         | 0.79     |
| CYCS    | 7.80 × 10^−3   | 77.8         | 0.80     |
| CEP63   | 2.80 × 10^−14  | 75.0         | 0.77     |
| FARP1   | 1.35 × 10^−22  | 83.3         | 0.88     |
| FRMD4A  | 8.53 × 10^−31  | 83.1         | 0.84     |
| RPL6    | 2.74 × 10^−36  | 67.7         | 0.77     |
| PFN1    | 3.91 × 10^−16  | 80.0         | 0.85     |

  
- **Gene Ranking Mechanism**: sc-OTGM's recommendation system ranks genes across selected cell types within 
each dataset. This ranking approach prioritizes genes that are more likely to shift the distribution of the cell state
from healthy to diseased.

- **Evaluation**: Classification accuracies, F1-scores, and confusion matrices were computed to evaluate the 
predictive power of sc-OTGM in distinguishing between cell types and states.


## Citation
Single-Cell Perturbation Modeling by Solving Optimal Mass Transport on the Manifold of 
Gaussian Mixtures. ICLR 2024 Workshop on Machine Learning for Genomics Explorations.

Corresponding author: andac.demir@novartis.com
