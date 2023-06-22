https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68062

GPL13304	Illumina HiSeq 2000 (Drosophila melanogaster)
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68062/matrix/

We sequenced mRNAs from individual unfertilized eggs, precisely staged and sexed blastoderm embryos and compared levels between D. melanogaster, D. yakuba, D. pseudoobscura and D. virilis.

GSM1662140	mel.M.1
GSM1662141	mel.M.2
GSM1662142	mel.M.3
GSM1662143	mel.F.1
GSM1662144	mel.F.2
GSM1662145	mel.F.3
GSM1662164	mel.U.2
GSM1662165	mel.U.3

methods part from 2015 publication

Six normalization procedures were tested: 
(i) either no further normalization after eXpress, 
normalization after eXpress to the median (ii), 
to the 75th percentile (iii), 
to the 95th percentile (iv), 
TMM normalization (v) or 
quantile normalization (see [42] for review). 

The analyses pre- sented in the main text were conducted on the 75th percentile normalization but the different normalizations essentially gave similar results. A few representative plots similar to the main figures are given in S10 and S11 Figs for the 6 normalizations. In addition, tables created from the different normalized datasets are available on GEO (accession number GSE68062).

GSE68062 $wc -l *txt
    6004 GSE68062_Gene_abundances_after_FPKM_normalization.txt
    6080 GSE68062_Gene_abundances_after_Q_normalization.txt
    6030 GSE68062_Gene_abundances_after_TMM_normalization.txt
    5999 GSE68062_Gene_abundances_after_p50_normalization.txt
    5994 GSE68062_Gene_abundances_after_p75_normalization.txt
    5999 GSE68062_Gene_abundances_after_p95_normalization.txt
    