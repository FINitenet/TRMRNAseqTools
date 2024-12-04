#!/usr/local/bin/R
library(tidyverse)
library(patchwork)

#### ARGS ####
args <- commandArgs(trailingOnly = T)
input_path <- args[1]
# input_path <- '/bios-store1/chenyc/Project/ZHB/Project_zhangbin_240108N/Annotation-type_dis'

#### defalut parameter ####
rnatype_colour <- c(spikein = '#1ae4d3',
                    unassigned = "#a6cee3",
                    snoRNA = "#33a02c",snRNA = "#fb9a99",rRNA = "#1f78b4",tRNA = "#b2df8a",
                    lncRNA = "#e31a1c",protein_coding = "#fdbf6f",
                    transposable_element = "#ff7f00",`hc-siRNA` = "#cab2d6",`phasi_tasiRNA` = "#6a3d9a",miRNA_primary_transcript = "#ffff99")
rnatype_level <- c("spikein",
                   "unassigned",
                   "snoRNA","snRNA","rRNA","tRNA",
                   "lncRNA","protein_coding",
                   "transposable_element","hc-siRNA","phasi_tasiRNA","miRNA_primary_transcript")

filelist <- list.files(path = input_path, pattern = "_priority_spikein.summary")

#### plot function ####
plot_fuc <- function(mt){
  df <- read.table(file.path(input_path,mt), header = TRUE, sep = '\t')
  df$feature <- gsub("otherRNA", "unassigned", df$feature)
  df$feature <- factor(df$feature, levels = rnatype_level)
  max <- max(df$length)
  min <- min(df$length)
  if (max > 50) {
    max <- 50
  }
  
  title <- gsub("_priority_spikein.summary","",mt)

  p <- ggplot(df) +
    aes(x = length, fill = feature, weight = Percentage) +
    geom_bar() +
    scale_fill_manual(values = rnatype_colour) +
    theme_classic()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          #axis.line.x = element_line(size = 1),#轴线粗细
          axis.text.y = element_text(size = 12),
          #axis.line.y = element_line(size = 1),#轴线粗细
          plot.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))+
    labs(x = "",y = "Relative percentage(%)",fill = "RNA type",title = title)+
    scale_x_continuous(limits = c(min-1,max+1), breaks = seq(min, max, 1))
}

#####

plots <- list()
for (file in filelist) {
  p <- plot_fuc(file)
  plots[[file]] <- p
}

plots_combined <- wrap_plots(plots, ncol = 2) + plot_layout(guides = 'collect')

ggsave(file.path(input_path, "rna_type_dis_spikein.pdf"), plots_combined, width = 20, height = 10)

