# theoretical_bounds_on_mutation_subsampling

## Installation

  
## Instructions for Generating Perfect Binary Phylogeny
### Input
- `-n` : number of leaves
- `-m` : number of mutations
- `-t` : random seed for tree generation
- `-s` : random seed for mutation placement
- `-o` : output_prefix
- 
### Outputs
- output_prefix_tree.dot: dot format tree with mutations as edge labels
- output_prefix_matrix.dot: mutation matrix of size (2n-1) * m, where n is number of leaves and m is number of mutations. 


### Instructions for Running
```
python generate_perfect.py -n 3 -m 10 -t 90 -s 90 -o "perfect" 
```
This will create two files named `perfect_tree.dot` and `perfect_matrix.dot`. 

