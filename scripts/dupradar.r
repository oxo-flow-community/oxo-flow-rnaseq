#!/usr/bin/env Rscript

# DupRadar duplication-rate QC for nf-core/rnaseq (port of
# modules/nf-core/dupradar/templates/dupradar.r @ 3.26.0).
#
# Written by Phil Ewels and released under the MIT license.
# Ported to nf-core/modules with template by Jonathan Manning.
# oxo-flow port: the Nextflow template variables are supplied via
# environment variables (BAM, GTF, PREFIX, THREADS, STRANDED, PAIRED).
# feature_type is fixed to 'exon' (the upstream default); the nf-core
# versions.yml output is not written (software versions are pinned in
# envs/dupradar.yaml and reported by the workflow-level versions file).

################################################
################################################
## Pull in module inputs                      ##
################################################
################################################

input_bam <- Sys.getenv("BAM")
output_prefix <- Sys.getenv("PREFIX")
annotation_gtf <- Sys.getenv("GTF")
threads <- as.integer(Sys.getenv("THREADS"))
feature_type <- 'exon'

stranded <- as.integer(Sys.getenv("STRANDED"))

paired_end <- Sys.getenv("PAIRED") == "true"

# Debug messages (stderr)
message("Input bam      : ", input_bam)
message("Input gtf      : ", annotation_gtf)
message("Strandness     : ", c("unstranded", "forward", "reverse")[stranded+1])
message("paired/single  : ", ifelse(paired_end, 'paired', 'single'))
message("feature type   : ", feature_type)
message("Nb threads     : ", threads)
message("Output basename: ", output_prefix)

# Load / install packages
library("dupRadar")
library("parallel")

# Duplicate stats
dm <- analyzeDuprates(input_bam, annotation_gtf, stranded, paired_end, threads, GTF.featureType = feature_type, verbose = TRUE)
write.table(dm, file=paste(output_prefix, "_dupMatrix.txt", sep=""), quote=F, row.name=F, sep="\t")

# 2D density scatter plot
pdf(paste0(output_prefix, "_duprateExpDens.pdf"))
duprateExpDensPlot(DupMat=dm)
title("Density scatter plot")
mtext(output_prefix, side=3)
dev.off()
fit <- duprateExpFit(DupMat=dm)
cat(
    paste("- dupRadar Int (duprate at low read counts):", fit$intercept),
    paste("- dupRadar Sl (progression of the duplication rate):", fit$slope),
    fill=TRUE, labels=output_prefix,
    file=paste0(output_prefix, "_intercept_slope.txt"), append=FALSE
)

# Create a multiqc file dupInt
sample_name <- gsub("Aligned.sortedByCoord.out.markDups", "", output_prefix)
line="#id: DupInt
#plot_type: 'generalstats'
#pconfig:
#    dupRadar_intercept:
#        title: 'dupInt'
#        namespace: 'DupRadar'
#        description: 'Intercept value from DupRadar'
#        max: 100
#        min: 0
#        scale: 'RdYlGn-rev'
Sample dupRadar_intercept"

write(line,file=paste0(output_prefix, "_dup_intercept_mqc.txt"),append=TRUE)
write(paste(sample_name, fit$intercept),file=paste0(output_prefix, "_dup_intercept_mqc.txt"),append=TRUE)

# Get numbers from dupRadar GLM
curve_x <- sort(log10(dm$RPK))
curve_y = 100*predict(fit$glm, data.frame(x=curve_x), type="response")
# Remove all of the infinite values
infs = which(curve_x %in% c(-Inf,Inf))
curve_x = curve_x[-infs]
curve_y = curve_y[-infs]
# Reduce number of data points
curve_x <- curve_x[seq(1, length(curve_x), 10)]
curve_y <- curve_y[seq(1, length(curve_y), 10)]
# Convert x values back to real counts
curve_x = 10^curve_x
# Write to file
line="#id: dupradar
#plot_type: 'linegraph'
#section_name: 'DupRadar'
#section_href: 'bioconductor.org/packages/release/bioc/html/dupRadar.html'
#description: \"provides duplication rate quality control for RNA-Seq datasets. Highly expressed genes can be expected to have a lot of duplicate reads, but high numbers of duplicates at low read counts can indicate low library complexity with technical duplication.
#    This plot shows the general linear models - a summary of the gene duplication distributions. \"
#pconfig:
#    title: 'DupRadar General Linear Model'
#    xlog: True
#    xlab: 'expression (reads/kbp)'
#    ylab: '% duplicate reads'
#    ymax: 100
#    ymin: 0
#    tt_label: '<b>{point.x:.1f} reads/kbp</b>: {point.y:,.2f}% duplicates'
#    x_lines:
#        - color: 'green'
#          dash: 'LongDash'
#          label:
#                text: '0.5 RPKM'
#          value: 0.5
#          width: 1
#        - color: 'red'
#          dash: 'LongDash'
#          label:
#                text: '1 read/bp'
#          value: 1000
#          width: 1"

write(line,file=paste0(output_prefix, "_duprateExpDensCurve_mqc.txt"),append=TRUE)
write.table(
    cbind(curve_x, curve_y),
    file=paste0(output_prefix, "_duprateExpDensCurve_mqc.txt"),
    quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE,
)

# Distribution of expression box plot
pdf(paste0(output_prefix, "_duprateExpBoxplot.pdf"))
duprateExpBoxplot(DupMat=dm)
title("Percent Duplication by Expression")
mtext(output_prefix, side=3)
dev.off()

# Distribution of RPK values per gene
pdf(paste0(output_prefix, "_expressionHist.pdf"))
expressionHist(DupMat=dm)
title("Distribution of RPK values per gene")
mtext(output_prefix, side=3)
dev.off()

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()
