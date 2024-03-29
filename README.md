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
```shell
python ./SPROUT-main/scripts/STORM_ST.py -s spatial_exp.tsv -v spatial_decon.tsv -r sc_exp.tsv -m sc_meta.tsv -l lr_pair.txt
```
```
usage: STORM_ST.py -s ST_FILE [-c ST_COORD] -v W_FILE -r SC_FILE -m META_FILE -l LR_FILE [-o OUT_DIR] [-p NUM_PER_SPOT] [-a MODE] [-h] 

optional arguments:
  -s ST_FILE, --st-file ST_FILE
                        Spatial transcriptomics data
  -c ST_COORD, --st-coordinate ST_COORD
                        Spatial coordinates of the spatial transcriptomics data
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
embedding(sparse_A, path, verbose = True, left_range, right_range, steps, dim)
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
    
### Input file format
* **Spatial Transcriptomics (ST) Count Data**
  * `.tsv` file with spots as rows and genes as columns

* **Cell-type deconvoluted spatial matrix**
  * `.tsv` file with spot as rows and cell-type as columns
 
* **Spatial coordinates**
  * `.tsv` file with spot as rows axis x and y as columns 

* **Single-cell RNA-seq Count Data**
  * `.tsv` file with cells as rows and genes as columns

* **Single-cell RNA-seq Metadata**
  * `.tsv` file with cells as rows and cell types as columns


### Reproducing The PDAC Analysis
#### 1. Downloading the data
The PDAC data we used in the paper is downloaded from [GSE111672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672), the processed data is deposited in [Google drive](https://drive.google.com/file/d/1r4hiv9z0HgmnXNRyYHk7FMg69hhFL579/view?usp=sharing).
- ST.zip
  - spatial_exp.tsv 
  - spatial_decon.tsv
  - sc_exp.tsv
  - sc_meta.tsv
#### 2. Run the analysis
```shell
python ./SPROUT-main/scripts/STORM_ST.py -s spatial_exp.tsv -v spatial_decon.tsv -r sc_exp.tsv -m sc_meta.tsv -l ./SPROUT-main/LR/human_LR_pairs.txt
```
```console
file loaded.
Cell type in weight matrix is unequal to single-cell meta file.

Start to select single-cell aggregates...
- Cell num per spot is 10, mode set as strict.
- 19738 genes in spatial data, 19738 genes in single-cell data.
- 14210 shared and expressed genes has been kept.
- Given the user-defined parameter keep_lr_per, STORM kept 80%, that is, 721 highly variable LR genes.
- STORM selects 3472 feature genes.
initial solution: 0.021963083441482775 0.4169362767129037 0.9069221853765786
swapped solution: 0.045146409164253054 0.4890899710706522 0.9419340778225115
Single-cell aggregates selection completed.

Start embedding...

```
### Output files
* **Selected single-cell profiels representing each spot**
  * `sc_agg.txt`  file with spots as rows and genes as columns

* **Metadata of the selected cells**
  * `sc_agg_index.txt`  file with cell as rows, cell type, and the spot its belongs to as columns 
 
* **Single-cell coordinates**
  * `coord_best_in_shape.csv` file with cells as rows, coordinates as columns 

## 2. *de novo* reconstruction from the single-cell transcriptome
```shell
python ./scripts/STORM_SC.py -r sc_exp.tsv -l lr_pair.txt
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

### Input file format
* **Single-cell RNA Sequencing Count Data**
  * `.tsv` file with genes as rows and cells as columns

### Reproducing The heart Analysis
#### 1. Downloading the data
The heart dataset we used in the paper is downloaded from [Developmental_heart_filtered_scRNA-seq_and_meta_data.zip](https://data.mendeley.com/datasets/mbvhhf8m62/2/files/33fb42ae-7b40-4a70-b61d-676f44d68d4c).

[Alternative storage cite](https://drive.google.com/drive/folders/1AuTPbBLGwK8qsrPiXzloxQW8wzwnZkxi?usp=drive_link)
- Developmental_heart_filtered_scRNA-seq_and_meta_data.zip
  - all_cells_count_matrix_filtered.tsv
  - all_cells_meta_data_filtered.tsv
#### 2. Run the analysis
```shell
python ./SPROUT-main/scripts/STORM_SC.py -r all_cells_count_matrix_filtered.tsv -l ./SPROUT-main/LR/human_LR_pairs.txt
```
