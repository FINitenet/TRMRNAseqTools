library(tidyverse)
library(ggplot2)
library(data.table)
library(openxlsx)

options(scipen = 999)  # 设置为极大值,禁用科学计数法

count_norm <- function(file, basename, xlab, norm_size) {
  # read count file
  # count <- fread("2_srna_analysis_1828/03.Quant/Annotation-pri-miRNA/miRNA.hairpin.annotation")
  count <- fread(file)
  count <- as.data.frame(count)
  
  for (i in 7:length(count)) {
    col_name <- paste0(names(count)[i], "_rpm")  # 构造新列名
    # count[col_name] <- count[, names(count)[i]] / sample[match(names(count)[i], sample$Sample), "Without_trsnoRNA_reads"] * 10^6
    count[col_name] <- count[, names(count)[i]] / sample[match(names(count)[i], sample$Sample), norm_size] * 10^6
  }
  
  # without reple
  count$log2FoldChange <- log2((count[,10]+1)/(count[,9]+1))
  count <- count %>% mutate(
    Sig = case_when(
      (.[, 7] > 10 | .[, 8] > 10) & log2FoldChange > 1 ~ "up",
      (.[, 7] > 10 | .[, 8] > 10) & log2FoldChange < -1 ~ "down",
      TRUE ~ "normal"
    ))
  
  
  p <- ggplot(count) +
    aes(
      x = get(xlab),
      y = log2FoldChange,
      colour = Sig
    ) +
    geom_point() +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed")+
    scale_color_brewer(palette = "PiYG", direction = -1) +
    scale_x_continuous(trans = "log10", name = xlab) +
    ggthemes::theme_few()+ 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
  ggsave(filename = paste0(file,".pdf"), plot = p)
  
  return(count)
}


summary <- fread("2_srna_analysis/03.Quant/Annotation-all.reads/summary.txt", sep = '\t')
summary <- as.data.frame(t(summary))
colnames(summary) <- summary[1, ]
summary <- summary[-1, ]
summary[] <- lapply(summary, as.numeric)
summary$sample <- gsub("_trimmed.bam", "", row.names(summary))

sample <- read.csv("2_srna_analysis/02.Mapping/mapping_results_ShortStack.csv",sep = '\t')
sample$type <- c(rep("mutant",6),"wildtype")
sample <- left_join(sample, summary, by = c("Sample" = "sample")) %>% 
  mutate(Without_trsnoRNA_reads = Mapped_reads - tRNA - snoRNA - snRNA - rRNA) %>%
  select(1, 7, 2:6, 19, everything())
# normalization
sample$sizefactors <- 10^6/sample$Mapped_reads

wb <- createWorkbook()
addWorksheet(wb, "library_stats")
writeData(wb, "library_stats", sample)

counts <- list.files(path = "2_srna_analysis/03.Quant/", pattern = "annotation$", recursive = TRUE, full.names = TRUE)
counts <- counts[!grepl("all.reads", counts)]
for (file in counts) {
  basename <- str_extract(file, "[^/]*$")
  print(file)
  addWorksheet(wb, basename)
  if (grepl("miRNA", file)) {
    tmp <- count_norm(file, basename, "hen1_8_240510N_R1", "Without_trsnoRNA_reads")
  }else{
    tmp <- count_norm(file, basename, "hen1_8_240510N_R1", "Mapped_reads")
  }
  writeData(wb, basename, tmp)
}
saveWorkbook(wb, "annotation_quant.xlsx", overwrite = TRUE)
