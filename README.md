# STORM
A software package restores the spatial structure from single-cell transcriptome with spectral graph sparsification.
# Pre-requirements
* python3
* numpy
* pandas
* scipy
* umap-learn
* loess
* lshashpy3
# Installation
## Installation from the source code
```shell
cd STORM
python setup.py install
```
# Usages
## Reconstruction utilizing single-cell transcriptome coupled with spatial transcriptomics data
```shell
python ./script/STORM_ST.py -s spatial_exp.tsv -v spatial_decon.tsv -r sc_exp.tsv -m sc_meta.tsv -l lr_pair.txt
```
```
usage: STORM_ST.py -s ST_FILE -v W_FILE -r SC_FILE -m META_FILE -l LR_FILE [-o OUT_DIR] [-p NUM_PER_SPOT] [-a MODE] [-h] 

optional arguments:
  -s ST_FILE, --st-file ST_FILE
                        Spatial transcriptomics data
  -v W_FILE, --weight-file W_FILE
                        Deconvoluted ST data
  -r SC_FILE, --sc-file SC_FILE
                        Single-cell candidate library of the corresponding ST tissue
  -m META_FILE, --meta-file META_FILE
                        Cell-type annotation of the single-cell candidate library
  -l LR_FILE, --lr-file LR_FILE
                        Ligand-receptor pair data
  -o OUT_DIR, --out-dir OUT_DIR
                        Output file path， default is the current working dir
  -p NUM_PER_SPOT, --cell-num-per-spot NUM_PER_SPOT
                        Estimated cell number per spot. Default is 10
  -a MODE, --selection_mode MODE
                        The choice of either gather cells primarily from the same type (strict) 
                        or from all cells (wild) in the candidate library
  -h, --help            show this help message and exit                      
```
### Parameters
```python
embedding(sparse_A, path, left_range, right_range, steps, dim)
```
* left_range : int, optional, default: 0

* right_range : int, optional, default: 30

    The index range for neighbor number, the actual neighbor number is (i+1)*15
    
* steps : int, optional, default: 30

    The iteration number for each neighbor

* dim : int, optional, default: 2

    The embedding dimension of the reconstruction
### Input file format
1. Spatial expression matrix - spot x gene

|  | Gene 1 | Gene 2 | ... | Gene n|
| :---: | :---: | :---: | :---: | :---: |
| Spot 1 | 0 | 2 | ... | 5|
| Spot 2 | 1 | 30 | ... | 0|
| Spot ... | 0 | 0 | ... | 60|
| Spot m | 30 | 0 | ... | 20|

2. Cell-type deconvoluted spatial matrix - spot x cell-type

|  | Type 1 | Type 2 | ... | Type c|
| :---: | :---: | :---: | :---: | :---: |
| Spot 1 | 0 | 0.2 | ... | 0.001|
| Spot 2 | 0.6 | 0.13 | ... | 0|
| Spot ... | 0 | 0 | ... | 0.146|
| Spot m | 0.12 | 0 | ... | 0.13|

3. Spatial coordinates - spot x cell-type

|  | Type 1 | Type 2 | ... | Type c|
| :---: | :---: | :---: | :---: | :---: |
| Spot 1 | 0 | 0.2 | ... | 0.001|
| Spot 2 | 0.6 | 0.13 | ... | 0|
| Spot ... | 0 | 0 | ... | 0.146|
| Spot m | 0.12 | 0 | ... | 0.13|
## *de novo* reconstruction from the single-cell transcriptome
```shell
python ./scripts/STORM_SC.py -r tpm.txt -l LR_pairs_add.txt -o /home/wanwang6/scratch/5.UMAP/1.spatial/1.data/2.HCC/
```
### Parameters 
```python
optimizers.sparsify_embedding(Q, affinitymat, path, percent = per, left_range = 0, right_range = 30, steps = 30, rep = 5, dim = 3)
```
Parameters
* percent : float, optional, default: 0.2
  
  The percentage of edges to be preserved
  
* left_range : int, optional, default: 0

* right_range : int, optional, default: 30

    The index range for neighbor number, the actual neighbor number is (i+1)*15
    
* steps : int, optional, default: 30

    The iteration number for each neighbor
    
* rep : int, optional, default: 5

    The number for sparsification repeats

* dim : int, optional, default: 3

    The embedding dimension of the reconstruction

