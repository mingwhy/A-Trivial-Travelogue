
## 资料整理
### 基因组annotation (genomic range and genomic feature)
annotate genomic region 的整体思路：用自己的inquire 和 已知公共数据库genomic feature做交集。
- [AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html)
- [Handling genomic data using Bioconductor](http://www.haowulab.org/teaching/bioc/genomicRanges.pdf), more tutorials https://www.haowulab.org//pages/teaching.html
- [annotatr: Making sense of genomic regions](https://www.bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html)
- [Annotating Genomic Ranges, hg19 and mm10](https://bioconductor.riken.jp/packages/3.14/workflows/vignettes/annotation/inst/doc/Annotating_Genomic_Ranges.html#mouse-mm10)

### 甲基化
- [学一学DNA甲基化芯片分析流程](https://www.jieandze1314.com/post/cnposts/227/)
- [一个甲基化芯片信号值矩阵差异分析的标准代码](https://github.com/jmzeng1314/methy_array)
- [A cross-package Bioconductor workflow for analysing methylation array data](http://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html) and [F1000 manuscript](https://pubmed.ncbi.nlm.nih.gov/27347385/)
- [解读GEO数据存放规律及下载](https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486063&idx=1&sn=156bee5397e979722b36b78284188538&scene=21#wechat_redirect) 
- [DNA-methylation-analysis](https://github.com/crazyhottommy/DNA-methylation-analysis)
- [Quantitative comparison of within-sample heterogeneity scores for DNA methylation data](https://github.com/MPIIComputationalEpigenetics/WSHPackage) 

### 甲基化时钟
- [DNA methylation age and the epigenetic clock by Steve Horvath](https://horvath.genetics.ucla.edu/html/dnamage/)
- [R codes for preprocessing Illumina Infinium DNA Methylation array data (from IDAT files) and epigenetic clocks](https://github.com/wt2015-github/Methylation-array) 
- [Chronological and gestational DNAm age estimation using different methylation-based clocks](https://rpubs.com/jrgonzalezISGlobal/methylclock)
- [DNA methylation indices of exposure and phenotype (meffonym)](https://github.com/perishky/meffonym/)
- [Methylation clock resources](https://github.com/yiluyucheng/dnaMethyAge) and [R methylclockData package](https://github.com/isglobal-brge/methylclockData) https://bioconductor.org/packages/release/data/experiment/vignettes/methylclockData/inst/doc/methylcockData.html

### 甲基化芯片or测序数据分析实例
- [Gene set enrichment analysis for genome-wide DNA methylation data](https://github.com/Oshlack/methyl-geneset-testing)
- [Epigenomics Workshop 2022](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/index.html)
- [GBM 450k dataset](https://jokergoo.github.io/cola_examples/GBM_450K/)
- [Tutorial on mapping WGBS data using Bismark](http://statisticalrecipes.blogspot.com/2015/05/tutorial-on-mapping-wgbs-data-using.html), more tutorials https://github.com/genomicsclass/colonCancerWGBS/blob/master/scripts/createObject.Rmd

