# Theoretical Estimates on the Expected Number of Mutations for Reconstructing Clonal Lineage Trees.

<img width="1156" height="249" alt="image" src="https://github.com/user-attachments/assets/190a9fca-3708-4ab1-abe8-06a5ef78cf9e" />


## Installation

  
## Instructions for Generating Perfect Binary Phylogeny
### Input
- `-n` : number of leaves
- `-m` : number of mutations
- `-t` : random seed for tree generation
- `-s` : random seed for mutation placement
- `-o` : output_prefix

  
### Outputs
- output_prefix_tree.dot: dot format tree with mutations as edge labels
- output_prefix_matrix.dot: mutation matrix of size (2n-1) * m, where n is number of leaves and m is number of mutations. 


### Instructions for Running
```
python generate_perfect.py -n 3 -m 10 -t 90 -s 90 -o "perfect" 
```
This will create two files named `perfect_tree.dot` and `perfect_matrix.dot`. 

