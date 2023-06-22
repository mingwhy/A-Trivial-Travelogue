
can you assign a sex label (female or male) to individual cells based on single-cell RNA-seq data alone?

# preclude

## adult fly single cell data
male: roX1, roX2
female: Yp1, Yp2, Yp3
Li, Hongjie, et al. "Fly Cell Atlas: A single-nucleus transcriptomic atlas of the adult fruit fly." Science 375.6584 (2022): eabk2432.

## human embryo single cell data
Petropoulos, Sophie, et al. "Single-cell RNA-seq reveals lineage and X chromosome dynamics in human preimplantation embryos." Cell 165.4 (2016): 1012-1026.
see Figure S1 and section “Calling Embryonic Sex”
2016_FigS1.png
chrY and chrX gene expression together 

## fly embryo single cell ATAC data
Reddington, James P., et al. "Lineage-resolved enhancer and promoter usage during a time course of embryogenesis." Developmental Cell 55.5 (2020): 648-664.
Extended Data Figure 5: Sex of individual cells identified by ratio of X-chromosome to autosomal reads.
2021_extended.figure5.png
X:non-X reads ratio

# main

## fly single embryo bulk data
see details in embryo.data_summary.xlsx

informative prior genes for sexing cells: Sxl, msl-2, roX1, roX2
use either SVM, RF, or TSP, could easily train a model（2019）and then apply it to the unseen 'holdout' test data to predict sex labels （other studies).

## fly single embryo single-cell data
see details in embryo.data_summary.xlsx

train a model using 2021 data with sex labels, then apply it to 2017 and Jay's datasets.

one problem: 2021 data was generated using embryos collected at 5-8hr.
training models of this stage output low confidence scores for cells in other embryo stages.

cehck out: Liu, Jialin, et al. "The hourglass model of evolutionary conservation during embryogenesis extends to developmental enhancers with signatures of positive selection." Genome Research 31.9 (2021): 1573-1581.

# recipe

## TSP
k–Top Scoring Pairs (kTSP) algorithms
use gene relationships (eg. gene A> gene B, gene C<gene D) to esblish a set of rules for classification.
Marzouka, Nour-Al-Dain, and Pontus Eriksson. "multiclassPairs: an R package to train multiclass pair-based classifier." Bioinformatics 37.18 (2021): 3043-3044.

used case: 
digital.sexing.cells/recipe/TSP/TSP_on_single.cell_embryo.cell_embryo
	00_embryo_integrate.data_chrY.R
	01_TSP_single.cell.R

## SVM
digital.sexing.cells/recipe/SVM 

## RF
digital.sexing.cells/recipe/RF
