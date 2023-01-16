# SPROUT
The software implementation of the method in 
[SPROUT: spectral sparsification helps restore the spatial structure at single-cell resolution](https://academic.oup.com/nargab/article/4/3/lqac069/6700709).

A software package restores the spatial structure from single-cell transcriptome data with spectral graph sparsification.

# Pre-requirements
* python3, numpy, pandas, scipy
* umap-learn, loess, lshashpy3
# Installation
## Installation from the source code
```shell
wget https://github.com/deepomicslab/SPROUT/archive/refs/heads/main.zip
unzip main.zip
cd SPROUT-main
python setup.py install
```
# Usages
## 1. Reconstruction utilizing single-cell transcriptome coupled with spatial transcriptomics data
```python
from src import sprout
sprout_obj = sprout.SPROUT(st_exp = st_exp, st_coord = st_coord, weight = st_decon, 
                           sc_exp = sc_exp, meta_df = sc_meta, cell_type_key = 'celltype',lr_df = lr_df, 
                           verbose= 1, save_path = save_path)
spot_cor,picked_index_df = sprout_obj.select_sc(num_per_spot = 10, mode = 'strict', max_rep = 1, repeat_penalty = 10)
sc_coord = sprout_obj.spatial_recon(left_range = 0, right_range = 10, steps = 1, dim = 2,max_dist = 1)
```
### Input file format
* **Spatial Transcriptomics (ST) Count Data**
  * `st_exp` dataframe with spots as rows and genes as columns
 
* **Spatial coordinates**
  * `st_coord` dataframe with spot as rows, axis x and y as columns 

* **Cell-type deconvoluted spatial matrix**
  * `st_decon` dataframe with spot as rows and cell-type as columns


* **Single-cell RNA-seq Count Data**
  * `sc_exp` dataframe with cells as rows and genes as columns

* **Single-cell RNA-seq Metadata**
  * `sc_meta` dataframe with cells as rows and cell types as columns
  * `cell_type_key` column name of the celltype identity in `sc_meta`

* **Ligand and Receptor Data**
  * `lr_df` dataframe with ligand-receptor pairs as rows, ligand, receptor and its weight as columns

* **Verbose parameter**
  * `lr_df` dataframe with ligand-receptor pairs as rows, ligand, receptor and its weight as columns
### Parameters
```python
spatial_recon(left_range = 0, right_range = 10, steps = 1, dim = 2,max_dist = 1)
```
* verbose : Boolean, default: True

    True: save coordinate results of all parameters.
    
    False: only save the results with the highest shape correlation with the ST structure.
    
* left_range : int, optional, default: 0

* right_range : int, optional, default: 30

    The index range for neighbor number, the actual neighbor number is (i+1)*15
    
* steps : int, optional, default: 30

    The iteration number for each neighbor

* dim : int, optional, default: 2

    The embedding dimension of the reconstruction
    


### Output files
* **Selected single-cell profiels representing each spot**
  * `sc_agg.txt`  file with spots as rows and genes as columns

* **Metadata of the selected cells**
  * `sc_agg_index.txt`  file with cell as rows, cell type, and the spot its belongs to as columns 
 
* **Single-cell coordinates**
  * `coord_best_in_shape.csv` file with cells as rows, coordinates as columns 
