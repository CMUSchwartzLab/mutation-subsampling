# Theoretical Estimates on the Expected Number of Mutations for Reconstructing Clonal Lineage Trees.

<img width="1156" height="249" alt="image" src="https://github.com/user-attachments/assets/190a9fca-3708-4ab1-abe8-06a5ef78cf9e" />


## Installation

  
## Instructions for Generating Perfect Binary Phylogeny.
### Input
- `-n` : number of leaves
- `-m` : number of mutations
- `-t` : random seed for tree generation
- `-s` : random seed for mutation placement
- `-o` : output_prefix

  
### Outputs
- output_prefix_tree.dot: dot format tree with mutations as edge labels
- output_prefix_matrix.dot: mutation matrix of size `(2n-1) * m`, where n is number of leaves and m is number of mutations. 


### Instructions for Running
```
python generate_perfect.py -n 3 -m 10 -t 90 -s 90 -o "perfect" 
```
This will create two files named `perfect_tree.dot` and `perfect_matrix.dot`. 

## Instructions for Generating K-Dollo Phylogeny from the Perfect Phylogeny. 
### Input
The first two arguments are for the perfect phylogeny cell (or clone) * mutation matrix and the perfect phylogeny in dot format. The remaining arguments are the following - 
- `-k` : K for K-Dollo loss.
- `--loss` : loss probability
- `--mut-base` : 0/1 (whether mutations are 0 or 1 indexed)
- `-A` : K-dollo helper matrix (K-dollo completion of B)
- `-B` : K-dollo mutation matrix (this is the input for reconstructing K-Dollo phylogenies)
- `--dot` : K-Dollo output tree in dot format
  
### Output
- output.A: K-dollo helper matrix (K-dollo completion of B)
- output.B: K-dollo mutation matrix (this is the input for reconstructing K-Dollo phylogenies)
- output.dot: K-Dollo output tree in dot format

### Instructions for Running

```
python generate_dollo_from_perfect.py perfect_matrix.tsv perfect_tree.dot -k 1 --loss 0.1 -s 90 --mut-base 0 -A onedollo.A -B onedollo.B --dot onedollo_tree.dot
```



