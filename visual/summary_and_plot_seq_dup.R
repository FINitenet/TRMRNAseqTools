library(tidyverse)
library(ggplot2)
library(data.table)
library(openxlsx)
library(bbum)

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
  
  # with reple
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


summary <- fread("Annotation-all.reads/summary.txt", sep = '\t')
summary <- as.data.frame(t(summary))
colnames(summary) <- summary[1, ]
summary <- summary[-1, ]
summary[] <- lapply(summary, as.numeric)
summary$sample <- row.names(summary)

sample <- read.csv("3_extract_umi/mapping_results_bowtie_3_extract_umi.csv",sep = '\t')
sample <- left_join(sample, summary, by = c("Sample" = "sample"))%>% 
  mutate(Without_trsnoRNA_reads = Reported - tRNA - snoRNA - snRNA - rRNA)
spkein <- read.csv("3_extract_umi/mapping_results_bowtie_4_map2spikein.csv", sep = '\t')
spkein$Sample <- gsub("_umi","",spkein$Sample)
sample <- left_join(sample, spkein, by = "Sample")



# normalization
sample$sizefactors <- 10^6/sample$Reported
sample$Sample <- gsub("_240716N","",sample$Sample)

wb <- createWorkbook()
addWorksheet(wb, "library_stats")
writeData(wb, "library_stats", sample)

counts <- fread("Annotation_ncRNA_gene/ncRNA_gene.annotation")

g <- counts[,-c(10,19,20,22)]
row_sums <- apply(g[, 7:ncol(g)], 1, sum)
g <- g[row_sums > 100, ]

for (i in 7:length(g)) {
  col_name <- paste0(names(g)[i], "_rpm")  # 构造新列名
  rpm_values <- g[[names(g)[i]]] / as.numeric(sample[match(names(g)[i], sample$Sample), 5]) * 10^6
  set(g, j = col_name, value = rpm_values)  # 通过引用添加新列
}


# 定义一个函数，生成分组列表并计算每组的行均值
generate_and_calculate_groups <- function(df, group_names, patterns) {
  # 生成分组列表
  groups <- list()
  
  for (group in group_names) {
    # 找到与当前组名匹配的所有模式
    matching_patterns <- grep(group, patterns, value = TRUE)
    
    # 使用这些模式匹配数据框中的列名
    matching_columns <- grep(paste(matching_patterns, collapse="|"), colnames(df), value = TRUE)
    
    # 将匹配的列名分配给相应的组
    groups[[group]] <- matching_columns
  }
  
  df_filtered <- df
  
  # 计算每个组的行均值，并添加到数据框中
  for (group in group_names) {
    df_filtered <- df_filtered %>%
      dplyr::mutate("{group}_avg" := rowMeans(select(., all_of(groups[[group]])), na.rm = TRUE))
  }
  
  return(df_filtered)
}

# 示例使用
group_names <- c("g_leaf", "h_leaf", "L_leaf", "g_root", "h_root", "L_root")
patterns <- c("g_leaf1_rpm","h_leaf1_rpm","L_leaf1_rpm", "g_root1_rpm", "h_root1_rpm", "L_root1_rpm", 
              "g_leaf2_rpm","h_leaf2_rpm","L_leaf2_rpm", "g_root2_rpm", "h_root2_rpm", "L_root2_rpm")

# 假设 df 是你的数据框
g_filtered <- generate_and_calculate_groups(g, group_names, patterns)

# 查看生成的结果
print(g_filtered)



g_filtered$leaf_g_h <- log2((g_filtered$g_leaf_avg+1)/(g_filtered$h_leaf_avg+1))
g_filtered$root_g_h <- log2((g_filtered$g_root_avg+1)/(g_filtered$h_root_avg+1))
g_filtered <- g_filtered %>% filter(!grepl("^gene:", Geneid))

addWorksheet(wb, "normal")
writeData(wb, "normal", g_filtered)

colorpalettes <- c("#5E4FA2","#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "red", "#9E0142")


p1 <- g_filtered %>% filter(leaf_g_h > 1) %>% 
ggplot(aes(x=log10(g_leaf_avg+1),y=log10(L_leaf_avg+1))) +
  geom_point(aes(color =leaf_g_h,size=leaf_g_h)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  scale_color_gradientn(colors =colorpalettes)+
  scale_x_continuous(breaks = seq(0,6,1),expand = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0.2,0))+
  scale_radius(range=c(1,4),guide=NULL)+
  labs(x = "log10(Mean of graft sample)", 
       y = "log10(Mean of Ler sample)",
       title = "leaf") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5))

p2 <- g_filtered %>% filter(root_g_h > 2) %>% 
  ggplot(aes(x=log10(g_root_avg+1),y=log10(L_root_avg+1))) +
  geom_point(aes(color =root_g_h,size=root_g_h)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  scale_color_gradientn(colors =colorpalettes)+
  scale_x_continuous(breaks = seq(0,6,1),expand = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0,0.5))+
  scale_radius(range=c(1,4),guide=NULL)+
  labs(x = "log10(Mean of graft sample)", 
       y = "log10(Mean of Ler sample)",
       title = "root") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5))

p <- p1+p2
p

P <- counts[,c(1:6,10,19,20,22)]
row_sums <- apply(P[, 7:ncol(P)], 1, sum)
P <- P[row_sums > 100, ]

for (i in 7:length(P)) {
  col_name <- paste0(names(P)[i], "_rpm")  # 构造新列名
  rpm_values <- P[[names(P)[i]]] / as.numeric(sample[match(names(P)[i], sample$Sample), 5]) * 10^6
  set(P, j = col_name, value = rpm_values)  # 通过引用添加新列
}
P <- P %>% filter(!grepl("^gene:", Geneid))
addWorksheet(wb, "P")
writeData(wb, "P", P)

namelist <- read.table("/bios-store1/chenyc/scripts/TRMRNAseqTools/list/MIRNA_AGI_MIRBASE.list")
df_filtered <- P[grepl("^MI", P$Geneid), ]
df_filtered <- left_join(df_filtered, namelist, by = c("Geneid"="V2"))

test <- df_filtered %>% 
  dplyr::mutate(label = ifelse(Ph_rpm %in% sort(Ph_rpm, decreasing = TRUE)[1:30], V3, "")) %>% 
  ggplot() +
    aes(x = log10(Ph_rpm), y = log10(hL_rpm)) +
    geom_point(colour = "#112446") +
    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
    theme_classic(base_size = 20) +
    labs(x = "log10(Ph)", 
         y = "log10(hL)") +
  ggrepel::geom_text_repel(aes(label = label), 
                           max.overlaps =30*10,
                           color = "#3288BD", 
                           vjust = 0.8, size = 4, alpha = 0.7)


saveWorkbook(wb, "annotation_quant.xlsx", overwrite = TRUE)
