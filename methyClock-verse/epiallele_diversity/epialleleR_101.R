#https://www.bioconductor.org/packages/release/bioc/html/epialleleR.html
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/epialleleR_1.6.1.tar.gz", repos = NULL, type="source")
library(epialleleR)

capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
bam.data    <- preprocessBam(capture.bam)
#> Reading BAM file [0.037s]
#> 
#> # data.table::data.table object for
# CpG VEF report
cg.vef.report <- generateCytosineReport(bam.data)
#> Already preprocessed BAM supplied as an input. Options 'min.mapq', 'min.baseq', 'skip.duplicates' and 'nthreads' will have no effect.
#> Thresholding reads [0.001s]
#> Preparing cytosine report [0.024s]
head(cg.vef.report[order(meth+unmeth, decreasing=TRUE)])

# Exploring DNA methylation patterns
# First, let's extract base methylation information for sequencing reads
# of 1:9 mix of methylated and non-methylated control DNA
patterns <- extractPatterns(
  bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
  bed=as("chr17:43125200-43125600","GRanges")
)
#> Reading BAM file [0.010s]
#> Extracting methylation patterns [0.027s]

# that many read pairs overlap genomic region of interest
nrow(patterns)
#> [1] 238

# these are positions of bases within relevant read pairs
base.positions <- grep("^[0-9]+$", colnames(patterns), value=TRUE)

# let's make a summary table with counts for every pattern
patterns.summary <- patterns[, c(lapply(.SD, unique), .N),
                             by=.(pattern, beta), .SDcols=base.positions]

# that many unique methylation patterns were found
nrow(patterns.summary)
#> [1] 45

# let's melt and plot patterns
plot.data <- data.table::melt.data.table(patterns.summary,
                                         measure.vars=base.positions, variable.name="pos", value.name="base")

# categorical positions, all patterns sorted by beta, with counts on the right
if (require("ggplot2", quietly=TRUE) & require("ggstance", quietly=TRUE)) {
  ggplot(na.omit(plot.data),
         aes(x=pos, y=reorder(pattern,beta),
             group=pattern, color=factor(base))) +
    geom_line(color="grey", position=position_dodgev(height=0.5)) +
    geom_point(position=position_dodgev(height=0.5)) +
    scale_colour_grey(start=0.8, end=0) +
    theme_light() +
    theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1),
          axis.text.y=element_blank()) +
    labs(x="position", y="pattern", title="epialleles", color="base") +
    scale_x_discrete(expand=c(0.05,0)) +
    geom_text(inherit.aes=FALSE, data=patterns.summary,
              mapping=aes(x="count", y=pattern, label=N), size=3)
}


# continuous positions, nonunique patterns according to their counts
if (require("ggplot2", quietly=TRUE) & require("ggstance", quietly=TRUE)) {
  ggplot(na.omit(plot.data)[N>1],
         aes(x=as.numeric(as.character(pos)), y=factor(N),
             group=pattern, color=factor(base))) +
    geom_line(color="grey", position=position_dodgev(height=0.5)) +
    geom_point(position=position_dodgev(height=0.5)) +
    scale_colour_grey(start=0.8, end=0) +
    theme_light() +
    labs(x="position", y="count", title="epialleles", color="base")
}

# upset-like plot of all patterns, continuous positions, sorted by counts
if (require("ggplot2", quietly=TRUE) & require("gridExtra", quietly=TRUE)) {
  grid.arrange(
    ggplot(na.omit(plot.data),
           aes(x=as.numeric(as.character(pos)), y=reorder(pattern,N),
               color=factor(base))) +
      geom_line(color="grey") +
      geom_point() +
      scale_colour_grey(start=0.8, end=0) +
      theme_light() +
      theme(axis.text.y=element_blank(), legend.position="none") +
      labs(x="position", y=NULL, title="epialleles", color="base"),
    
    ggplot(unique(na.omit(plot.data)[, .(pattern, N, beta)]),
           aes(x=N+0.5, y=reorder(pattern,N), alpha=beta, label=N)) +
      geom_col() +
      geom_text(alpha=0.5, nudge_x=0.2, size=3) +
      scale_x_log10() +
      theme_minimal() +
      theme(axis.text.y=element_blank(), legend.position="none") +
      labs(x="count", y=NULL, title=""),
    ncol=2, widths=c(0.75, 0.25)
  )
}


# now let's explore methylation patterns in RAD51C gene promoter using
# methylation capture data
capture.patterns <- extractPatterns(
  bam=system.file("extdata", "capture.bam", package="epialleleR"),
  bed=as("chr17:58691673-58693108", "GRanges"),
  verbose=FALSE
)
capture.positions <- grep("^[0-9]+$", colnames(capture.patterns), value=TRUE)
capture.summary <-
  capture.patterns[, c(lapply(.SD, unique), .N),
                   by=.(pattern, beta), .SDcols=capture.positions]
capture.data <- data.table::melt.data.table(capture.summary,
                                            measure.vars=capture.positions, variable.name="pos", value.name="base")
if (require("ggplot2", quietly=TRUE) & require("ggstance", quietly=TRUE)) {
  ggplot(na.omit(capture.data),
         aes(x=as.numeric(as.character(pos)), y=factor(N),
             group=pattern, color=factor(base))) +
    geom_line(color="grey", position=position_dodgev(height=0.9)) +
    geom_point(position=position_dodgev(height=0.9)) +
    scale_colour_grey(start=0.8, end=0) +
    theme_light() +
    labs(x="position", y="count", title="RAD51C", color="base")
}


