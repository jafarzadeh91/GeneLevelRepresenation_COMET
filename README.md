## Estimation of Gene-Level Represenation for DNA Methtylation Data 
In this project, we infer the latent represenation of genes and the gene-gene interaction network in human genome using the DNA methylation profile of subjects in a cohort. 
## Development Environment 
[![Matlab(Base_Optimization toolbox)](https://img.shields.io/badge/Matlab_core+optimization_toolbox-2017b-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![QUIC](https://img.shields.io/badge/QUIC-1.1-blue.svg)](http://bigdata.ices.utexas.edu/software/1035/)
[![macOS Mojave](https://img.shields.io/badge/macOS_Mojave-10.14-blue.svg)](https://www.apple.com/ca/macos/mojave/)

## Installation
Add the project folder in additon to your `QUIC` installation folder to your Matlab path.

## Functions

### Data Simulation

#### _Generate Gamma Matrix_
```r
Gamma = 
    probe_gene_distance_matrix_to_gamma(probe_gene_distance_matrix, chunk_size)
```
**or**
```r
Gamma = 
    probe_gene_distance_matrix_to_gamma(probe_gene_distance_matrix)
```
**inputs**: 
- `probe_gene_distance_matrix`: a sparse `number_of_cpg_sites*number_of_genes` matrix representing the base-pair distance (signed or absolute value) between CpG sites and TSS of genes. 
- `chunk_size`: It shows the number of CpG sites for which we can have the full(not sparse) matrix of distances in the memory (default: 30000). Adjust it according to your memory. Bigger value results in faster execution of function but needs more available memory.

**output**:
- `Gamma`: a `number_of_cpg_sites*number_of_genes` sparse matrix, each row represents the distribution for the distal proximity of a CpG site to different genes (values are non-negative and each row sums up to 1.).

#### _Generate K Matrix_
```r
K = generate_sparse_psd_block_matrix(number_of_genes,number_of_blocks,value_cut_off,min_eigen_value )
```
or
```r
K = generate_sparse_psd_block_matrix(number_of_genes,number_of_blocks,value_cut_off)
```
or
```r
K = generate_sparse_psd_block_matrix(number_of_genes,number_of_blocks)
```
**description**:
This function generate positive definite matrix consisting of blocks(modules) of highly connected genes.

**inputs**: 
- `number_of_genes`
- `number_of_blocks`: number of blocks of highly-connected genes in the K matrix.
- `value_cut_off`: minimum absolute value for the non-zero elements of K (default: 0.5).
- `min_eigen_value`: the value of the smallest eigen-value of K matrix (default: 0.5).

**output**:
- `K`: a `number_of_genes*number_of_genes` symmetric positive definite matrix of gene-gene interactions.

#### _Generate W Matrix_
```r
W=W_initialiation(Gamma, beta)
```
**description**:
This function generate the W matrix from the generative model using the Gamma matrix and hyper-parameter of beta.

**inputs**: 
- `Gamma`
- `beta`: The [![beta](https://latex.codecogs.com/gif.latex?%5Cbeta)]() in the model to regularize the weights values.

**output**:
- `W`:a `number_of_cpg_sites*number_of_genes` sparse matrix representing the weights between pairs of CpG sites and genes in their window.

#### _Generate L Matrix_
```r
L = mvnrnd(zeros(number_of_subjects, number_of_genes), inv(K));
```
**description**:
`mvnrnd` is a part of matlab standard libraries. We mentioned here to explain the way we can generate L in the process of simulating data.

#### _Generation Data(i.e. X matrix)_
```r
X = 
    data_simulator(Gamma, W, L, sigma2, number_of_subjects, number_of_genes)
```

**description**:
X is the "dataset" tha be used as input to our inference algorithm or any other similar method. It represents the value of CpG sites for different subjects in the cohort.
**inputs**:
- `Gamma`
- `W`
- `L`
- `sigma2`: the [![](https://latex.codecogs.com/gif.latex?%5Csigma%5E2)]() hyper-parameter in the model.
- `number_of_subjects`
- `number_of_genes`

**output**:
- `X`: the `number_of_subjects*number_of_cpg_sites` matrix where each column shows the simulated value for a CpG site across all subjects.


### Optimization

#### _L Initialization_
```r
L_initial = L_initialization_value(Gamma,X,weighted)
```
or
```r
L_initial = L_initialization_value(Gamma,X)
```
**description**:
This function initializes the gene level representation for different genes (represented by columns of L). For each column `i` of L, it  is initialized by the first principle component of PCA (or weighted PCA where weighted are the inverse of distances) of the CpG sites in the window of gene `i`.
**inputs**:
- `Gamma`
- `X`
- `weighted`: set `true` for usage of weighted PCA algorithm and `false` for ordinary PCA to infer the initial representations for genes.

**output**:
- `L_initial`

#### _Inference_
```r
[L_in_diff_iters,  K_in_diff_iters, W_in_diff_iters, objective_in_diff_iters] = 
infer(X, Gamma, L_initial, beta, eta, lambda, number_of_iters, opt_func)
```

**inputs**:
- `X`
- `Gamma`
- `L_initial`
- `beta`: [![beta](https://latex.codecogs.com/gif.latex?%5Cbeta)]()
- `sigma2`: [![sigma2](https://latex.codecogs.com/gif.latex?%5Csigma%5E2)]()
- `lambda`: [![lambda](https://latex.codecogs.com/gif.latex?%5Clambda)]()
- `number_of_iters`: number of iterations in coordinate descent algorithm.
- `opt_func`: algorithm to use for covariance matrix estimation: `bigquic`, `quic`, `bcd`. 

**output**:
- `L_in_diff_iters`
- `K_in_diff_iters`
- `W_in_diff_iters`
- `objective_in_diff_iters`

## Examples
In this example, simulate DNA methylation on chromosome number 22 for 500 subjects.

```r
%% Data Generation

probe_gene_distance_matrix_KoborDNAm_100_kb=load('probe_gene_distance_matrix_KoborDNAm_100_kb.mat');
MSK_probes_subjects = load('MSK_probes_subjects.mat');

probe_gene_distance_matrix = probe_gene_distance_matrix_KoborDNAm_100_kb.probe_gene_distance_matrix;

probe_gene_distance_matrix=
    probe_gene_distance_matrix(probe_gene_distance_matrix_KoborDNAm_100_kb.rows_metadata.CHR=='22',:);

probe_gene_distance_matrix(:,sum(probe_gene_distance_matrix,1)==0)=[];
probe_gene_distance_matrix(sum(probe_gene_distance_matrix,2)==0,:)=[];

number_of_genes = size(probe_gene_distance_matrix,2);
number_cpg_sites = size(probe_gene_distance_matrix,1);

Gamma = probe_gene_distance_matrix_to_gamma(probe_gene_distance_matrix);
K = generate_sparse_psd_block_matrix(number_of_genes,10,0.5,0.5 );
imagesc(abs(K));
W=W_initialiation(Gamma, 2);
number_of_subjects=500;
L = mvnrnd(zeros(number_of_subjects, number_of_genes), inv(K));
X = data_simulator(Gamma, W, K, 0.01, number_of_subjects, number_of_genes);

%% Optimization
L_initial = L_initialization_value(Gamma,X,true);

[L_in_diff_iters,  K_in_diff_iters, W_in_diff_iters, objective_in_diff_iters] = 
    infer(X, Gamma, L_initial, 1,0.1, 0.1, 3, 'bcd');

```
## Authors' Info
Sina Jafarzadeh. SJafarzadeh@cmmt.ubc.ca
