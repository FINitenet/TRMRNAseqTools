library(ReportingTools)
library(RNAseqQC)
library(DESeq2)
library(ensembldb)
library(ggplot2)
library(purrr)
library(tibble)
library(magrittr)
library(AnnotationHub)
library(ggplotify)

# library stats
htmlRep <- HTMLReport(shortName = "librart_status", reportDirectory = "./reports")
publish("library stats summary", htmlRep)
publish(meta, htmlRep)

total_counts <- plot_total_counts(dds)
publish("total counts", htmlRep)
publish(total_counts, htmlRep)


library_complexity <- plot_library_complexity(dds)
publish("library complexity", htmlRep)
publish(library_complexity, htmlRep)

gene_detection <- plot_gene_detection(dds)
publish("gene detection", htmlRep)
publish(gene_detection, htmlRep)


# not recommand
vsd <- vst(dds)
mean_sd_plot <- mean_sd_plot(vsd)
publish("mean sd plot", htmlRep)
publish(mean_sd_plot, htmlRep)

# map(c("1", "5"), ~plot_chromosome(vsd, .x))

# colData(vsd)$trt_mut <- paste0(colData(vsd)$treatment, "_", colData(vsd)$mutation)
# ma_plots <- plot_sample_MAs(vsd, group = "trt_mut")
# cowplot::plot_grid(plotlist = ma_plots[17:24], ncol = 2)

set.seed(1)
sample_clustering <- as.ggplot(plot_sample_clustering(vsd, anno_vars = c("group", "ecotype", "rep")))
publish("sample clustering", htmlRep)
publish(sample_clustering, htmlRep)

pca <- plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "group")
publish("sample pca", htmlRep)
publish(pca, htmlRep)

finish(htmlRep)


des2Report <- HTMLReport(shortName = 'sRNAseq_analysis_with_DESeq2',
                         title = 'sRNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./reports")

publish(dds,des2Report, pvalueCutoff=0.05, factor = colData(dds)$group,
        reportDir="./reports")
finish(des2Report)









