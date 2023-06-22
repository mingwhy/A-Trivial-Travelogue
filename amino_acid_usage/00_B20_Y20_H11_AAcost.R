# supp table 1 from https://www.nature.com/articles/s41467-018-06461-1
# import data from PDF
# https://crimebythenumbers.com/scrape-table.html
#BiocManager::install('pdftools')
library(pdftools)

df= pdf_text("suppTable1_raw.pdf")
head(df)
length(df)

sector_profile <- df[1]
sector_profile <- strsplit(sector_profile, "\n")
sector_profile <-sector_profile[[1]]
sector_profile[6:8]
header=sector_profile[6]

aa=sector_profile[9:29]
aa[12] 
aa=aa[-12]
length(aa)
aa.list=lapply(aa,function(x) unlist(strsplit(x,'\\s+')))
sapply(aa.list,length)
aa.df=as.data.frame(Reduce(`rbind`,aa.list))
head(aa.df)
aa.df$V1='EAA'
aa.df[12,]
aa.df[1:11,]$V1='NEAA'
colnames(aa.df)=c('essential','AA','AA.short','AA.abbreviation','Decay.rate','b20','y20','h11','B20','Y20','H11')
head(aa.df)
str(aa.df)
for(i in 5:11){
  aa.df[,i]=as.numeric(aa.df[,i])
}
aa.df
data.table::fwrite(aa.df,'cost_per_AA.txt')
