#https://bernatgel.github.io/karyoploter_tutorial//Examples/CpGIslands/CpGIslands.html

library(AnnotationHub)
ahub <- AnnotationHub()
ahub=ahub[ahub$species=='Mus musculus',]
str(ahub) #1649
ahub[grep('CpG',ahub$title,ignore.case = T),] #

ahub["AH6117"]
cpgs <- ahub[["AH6117"]]
cpgs #16026 ranges