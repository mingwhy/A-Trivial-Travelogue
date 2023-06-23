# https://elifesciences.org/articles/40675
# https://github.com/gr-meer/WLMT
# https://github.com/kerepesi/MouseAgingClocks
library(readxl)
library(stringr)

# read in Meer sample age info (https://elifesciences.org/articles/40675/figures#supp1)
dat=read_excel('Meer_elife-40675-supp1-v2.xlsx',sheet='Results on the clock sets',skip=1)
table(dat$Paper)
dat=dat[dat$Paper=='This paper',]
dim(dat) #81 samples 
table(dat$`Age, days`) 
#180 300 360 600 900 
#20   2  20  20  19 
