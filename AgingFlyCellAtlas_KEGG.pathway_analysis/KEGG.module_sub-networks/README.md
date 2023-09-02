# Metabolica
A systems biology approach to prioritize METABOLItes in CAncer.

Here you can find the R script of the Metabolica tool and some toy examples of its usage. 
Metabolica is descendant of our previous tool Metabolizer [[1](http://cancerres.aacrjournals.org/content/78/21/6059), [2](https://www.nature.com/articles/s41540-019-0087-2), [3](http://metabolizer.babelomics.org/), [4](https://github.com/babelomics/metabolizer)] which was limited with the KEGG metabolic module activities.
Metabolica extends its application to whole metabolism. Dissects KEGG metabolic pathways into subpathways and calculates their activities which account for metabolite abundances.

Metabolica uses transcriptomic and/or genomic data as input to predict reaction activities and then propagates metabolic flux over metabolic hypergraph that accounts for the production of a metabolite. Metabolica returns individual-level results that can be used in wide-range of posterior analysis to explain state-of-the-art cell biology and complex cellular mechanisms of diseases.
