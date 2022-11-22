
df <- data.frame(genes   = c("ubiquitin-dependent protein catabolic process", "regulation of pinocytosis", "endosome", "endoplasmic reticulum", "regulation of amyloid beta formation", "autophagy", "endocytosis"),
                 UpRegulated = c(83, 2, 74, 73, 10 ,55,52),
                 DownRegulated = c(116, 5, 128, 114,8,134,102 ))

library(tidyr)
library(ggplot2)

df %>% 
  pivot_longer(cols = c("UpRegulated", "DownRegulated")) %>%
  ggplot(aes(name, value, fill = genes, alpha = name)) +
  geom_col(position = position_dodge(), color = "black") +
  scale_alpha_manual(values = c(0.5, 1), guide = guide_none()) +
  facet_grid(~genes, scales = "free_x", switch = "x") +
  theme(strip.placement  = "outside",
        panel.spacing    = unit(0, "points"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 12)) +
  labs(x = "Up and down regulated genes", y = "Number of genes")