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
## Installation with pip
```shell
pip install STORM-sc
```
## Installation from the source code
```shell
cd STORM
python setup.py install
```
# Usage
## Reconstruction utilizing single-cell transcriptome coupled with spatial transcriptomics data
```shell
python ./script/STORM_ST.py -s st-hippo2.tsv -v W.hippo2.tsv -r cnt_data.tsv -m mta_data.tsv -l mouse_lr_pair.txt
```
```
usage: STORM_ST.py -s ST_FILE -v W_FILE -r SC_FILE -m META_FILE -l LR_FILE [-o OUT_DIR] [-p NUM_PER_SPOT] [-a MODE] [-h] 

optional arguments:
  -s ST_FILE, --st-file ST_FILE
                        Spatial transcriptomics data
  -v W_FILE, --weight-file W_FILE
                        Deconvoluted ST data
  -r SC_FILE, --sc-file SC_FILE
                        Single-cell candidate library of the corresponding ST
                        tissue
  -m META_FILE, --meta-file META_FILE
                        Cell-type annotation of the single-cell candidate
                        library
  -l LR_FILE, --lr-file LR_FILE
                        Ligand-receptor pair data
  -o OUT_DIR, --out-dir OUT_DIR
                        Output file path， default is the current working dir
  -p NUM_PER_SPOT, --cell-num-per-spot NUM_PER_SPOT
                        Estimated cell number per spot. Default is 10
  -a MODE, --selection_mode MODE
                        The choice of either gather cells only from the same
                        type (strict) or from all cells (wild) in the
                        candidate library
  -h, --help            show this help message and exit                      
```
## *de novo* reconstruction from the single-cell transcriptome
```shell
python ./scripts/STORM_SC.py -r tpm.txt -l LR_pairs_add.txt -o /home/wanwang6/scratch/5.UMAP/1.spatial/1.data/2.HCC/
```
### parameters 
```python
embedding(aff, out_dir, left_range = 0, right_range = 30, rep=30,dim=3)
```
parameters
* left_range : int, optional, default: 0

* right_range : int, optional, default: 30

    Range for neighbor number
    
* rep : int, optional, default: 30

    Iteration number for each neighbor

* dim : int, optional, default: 3

    Embedding dimension of the reconstruction

