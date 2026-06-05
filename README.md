# Theoretical Estimates on the Expected Number of Mutations for Reconstructing Clonal Lineage Trees.

<img width="1156" height="249" alt="image" src="https://github.com/user-attachments/assets/190a9fca-3708-4ab1-abe8-06a5ef78cf9e" />


## Installation
The codes are written in Python 3. The required packages are -
- numpy
- pydot
  
We have provided a `requirements.txt` file for the required packages. 

To install the benchmarking methods, we refer to the repositories of [SphyR](https://github.com/elkebir-group/SPhyR), [CellCoal](https://github.com/dapogon/cellcoal), and [PAUP](https://paup.phylosolutions.com/). 
  
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
python src/generate_perfect.py -n 3 -m 10 -t 90 -s 90 -o "perfect" 
```
This will create two files named `perfect_tree.dot` and `perfect_matrix.dot`. 

## Instructions for Generating K-Dollo Phylogeny from the Perfect Phylogeny. 
### Input
The first two arguments are for the perfect phylogeny cell (or clone) * mutation matrix and the perfect phylogeny in dot format. The remaining arguments are the following - 
- `-k` : K for K-Dollo loss.
- `--loss` : loss probability
- `--mut-base` : 0/1 (whether mutations are 0 or 1 indexed)
- `-A` : K-dollo helper matrix file in tsv format (K-dollo completion of B)
- `-B` : K-dollo mutation matrix file in tsv format (this is the input for reconstructing K-Dollo phylogenies)
- `--dot` : K-Dollo output tree in dot format
  
### Output
- output.A: K-dollo helper matrix (K-dollo completion of B)
- output.B: K-dollo mutation matrix (this is the input for reconstructing K-Dollo phylogenies)
- output.dot: K-Dollo output tree in dot format

### Instructions for Running

```
python src/generate_dollo_from_perfect.py perfect_matrix.tsv perfect_tree.dot -k 1 --loss 0.1 -s 90 --mut-base 0 -A onedollo.A -B onedollo.B --dot onedollo_tree.dot
```

## Instructions for Generating Coalescent Trees.

We generate coalescent trees using [CellCoal](https://github.com/dapogon/cellcoal). The command for running the program we use is -

```
./cellcoal-1.2.0 -n10 \
    -s"num_leaves" \
    -l10000 \
    -e100000 \
    -g1.0e-05 \
    -j3000 \
    -k1 \
    -i1 \
    -b0 \
    -c0 \
    -C5 \
    -u1.0e-07 \
    -f0.3 0.2 0.2 0.3 \
    -r0.00 0.03 0.12 0.04 0.11 0.00 0.02 0.68 0.68 0.02 0.00 0.11 0.04 0.12 0.03 0.00 \
    -1 -2 -3 -4 -6 -v -x -W \
    -o"output_directory" \
    -#200011
```
Description of the parameters: 
- `-n`: number of replicates
- `-l` : total sites in the genome
- `-e` : population size
- `-g` : exponential growth rate
- `-j` : number of sampled site
- `-k` : root branch length ratio
- `i` : rate variation among branches
- `b` : alphabet (0 for binary)
- `c` : germline mutation rate
- `C` : sequencing coverage
- `u` : mutation rate per site per generation 
- `f` : base frequencies
- `r` : mutation matrix ACGT x ACGT
-1 -2 -3 -4 -6 -v -x -W \
- `o` : output_directory \
- `#`: random seed
