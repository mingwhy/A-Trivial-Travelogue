# dataset
## sc or sn RNA-seq	dataset	stage	sex.label	#embryo	#cells	platform	technique	ref
- 2017	6	no	pool 100~200	1297	drop-seq	single-cell	Karaiskos, Nikos, et al. "The Drosophila embryo at single-cell transcriptome resolution." Science 358.6360 (2017): 194-199.
- 2021	11-12	yes	unknown	11001 female cells + 7222 male cells	10X	single-cell	Li, Yi-Ru, et al. "Trajectory mapping of the early Drosophila germline reveals controls of zygotic activation and sex differentiation." Genome research 31.6 (2021): 1011-1023.
- Jay	1~18	no	unknown	497921	sci-RNA-seq	single-nucleus	https://shendure-web.gs.washington.edu/content/members/DEAP_website/public/

# strategy
train a model using 2021 data with sex labels, then apply it to 2017 and Jay's datasets.
one problem: 2021 data was generated using embryos collected at 5-8hr.
training models of this stage may output low confidence scores for cells in other embryo stages.

1) date preprocessing
- 2021, /Users/mingyang/Documents/Data_fly_FCA/embryo_germline/0_*.R, three R scripts
- 2017, /Users/mingyang/Documents/Data_fly_FCA/embryo/0_*.R, two R scripts 
- Jay, /Users/mingyang/Documents/Data_Jay_fly_development/RNA_seurat_object/pred_windows/*.rds  

2) 

####################################################################################################################
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

Calderon, Diego, et al. "The continuum of Drosophila embryonic development at single-cell resolution." Science 377.6606 (2022): eabn5800.
SuppMethod: Inferring nuclei sex
2022_Jay_Inferring.nuclei.sex.png

# main

## fly single embryo single-cell data
see details in embryo.data_summary.xlsx

train a model using 2021 data with sex labels, then apply it to 2017 and Jay's datasets.

one problem: 2021 data was generated using embryos collected at 5-8hr.
training models of this stage output low confidence scores for cells in other embryo stages.

cehck out: Liu, Jialin, et al. "The hourglass model of evolutionary conservation during embryogenesis extends to developmental enhancers with signatures of positive selection." Genome Research 31.9 (2021): 1573-1581.

# recipe

## SVM
digital.sexing.cells/recipe/SVM 

## TSP
k–Top Scoring Pairs (kTSP) algorithms
use gene relationships (eg. gene A> gene B, gene C<gene D) to esblish a set of rules for classification.
Marzouka, Nour-Al-Dain, and Pontus Eriksson. "multiclassPairs: an R package to train multiclass pair-based classifier." Bioinformatics 37.18 (2021): 3043-3044.

used case: 
digital.sexing.cells/recipe/TSP/TSP_on_single.cell_embryo.cell_embryo
	00_embryo_integrate.data_chrY.R
	01_TSP_single.cell.R

## RF
digital.sexing.cells/recipe/RF
