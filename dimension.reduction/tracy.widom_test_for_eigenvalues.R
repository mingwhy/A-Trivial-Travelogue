##################################################################
# R package `AssocTests` https://rdrr.io/cran/AssocTests/man/tw.html
library(AssocTests)
out=tw(eigenvalues = c(5, 3, 1, 0), eigenL = 4, criticalpoint = 2.0234)
out$statistic
out$method
out$SigntEigenL # choose PC1-4

#tw funciton: https://rdrr.io/cran/AssocTests/man/tw.html
#a numeric value corresponding to the significance level. If the significance level is 
# 0.05, 0.01, 0.005, or 0.001,
# the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly. The default is 2.0234.

##################################################################
# R package `LEA` https://www.rdocum entation.org/packages/LEA/versions/1.4.0/topics/tracy.widom
#BiocManager::install("LEA")
library(LEA)
# Creation of the genotype file "genotypes.lfmm"
# with 1000 SNPs for 165 individuals.
data("tutorial")
write.lfmm(tutorial.R,"genotypes.lfmm")

pc = pca("genotypes.lfmm", scale = TRUE)
pc

tw = tracy.widom(pc)
head(tw)

# Plot the percentage of variance explained by each component.
plot(tw$percentage)
cumsum(tw$percentage)
# Display the p-values for the Tracy-Widom tests. 
tw$pvalues
which(tw$pvalues<0.05)

# remove pca Project
remove.pcaProject("genotypes.pcaProject")
