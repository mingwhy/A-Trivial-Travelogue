algorithms related to information theory


- Information: Information theory-based analyses to find a minimal non-redundant gene set for sex label classification
- Information: Calculate the information content of a GO term given a cell classificaiton scheme

## paper
Li, Hongjie, et al. "Classifying Drosophila olfactory projection neuron subtypes by single-cell RNA sequencing." Cell 171.5 (2017): 1206-1220.

github: https://github.com/felixhorns/FlyPN/tree/master/analysis

GH146_CombinatorialCode_LeaveOneOutValidation_Fig7.html

GH146_CombinatorialCode_Fig7.html


# information theory

r code demo

# a minimal code for sex label

# `normalized MI`

I've implemented a iterative PCA approach to remove multi-batch effects in analyzing metablome data from DAG (dog aging project). When accessing the 'contribution' of each covar to PCs, I noticed covars with a higher number of categories tend to have larger MI value with PCs. N

Now I realized it may be the inherent `bias` in calculating MI, inspired from a paper titled '[2015 Discovering What Dimensionality Reduction Really Tells Us About RNA-Seq Data](https://www.liebertpub.com/doi/10.1089/cmb.2015.0085)' from Prof Bonnie Berger.

Thus, one way to adjust this bias might be: MI.adj = MI(PCX,COVAR)/H(COVAR)
