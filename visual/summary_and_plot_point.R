library(tidyverse)
library(patchwork)

list <- read.table("/bios-store1/chenyc/scripts/TRMRNAseqTools/list/MIRNA_AGI_MIRBASE.list",sep = '\t')

df <- read.csv("res_degs/cleandata_normalized_counts.csv")
subdf <- df %>% select(1, ends_with("RPM"))
subdf <- left_join(list, subdf, by = c("V2"="Geneid"))

p1 <- ggplot(subdf) +
  aes(x = At_H_R_RPM, y = At_Za_R_RPM) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # 对角线
  geom_point(colour = "#112446", size = 0.5) +

  scale_y_continuous(trans = "log10", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = function(x) parse(text = paste0("10^", log10(x)))) +
  scale_x_continuous(trans = "log10", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = function(x) parse(text = paste0("10^", log10(x)))) +
  labs(x = "At_L_R_RPM")+
  theme_bw()+
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 9),
        panel.grid=element_blank())

subdf$'log2(L_R/Za_R)' <- log2((subdf$At_H_R_RPM+1)/(subdf$At_Za_R_RPM+1))
subdf <- subdf %>% mutate(flag = if_else(abs(`log2(L_R/Za_R)`) > 1, "a", "b"),
                          id = if_else(flag == "a", V3, ""))

p2 <- ggplot(subdf) +
  aes(x = At_H_R_RPM, y = `log2(L_R/Za_R)`, color = flag) +
  geom_hline(yintercept = c(1,-1), linetype=2, size=0.25)+ 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("red","grey"))+
  scale_x_continuous(trans = "log10", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = function(x) parse(text = paste0("10^", log10(x)))) +
  labs(x = "At_L_R_RPM")+
  theme_bw()+
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 9),
        panel.grid=element_blank())

p <- p1+p2

