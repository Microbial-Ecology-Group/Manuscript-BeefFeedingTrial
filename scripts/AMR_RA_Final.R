#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################

beta_data_noSNP <- data_noSNP
beta_data_noSNP <- prune_taxa(taxa_sums(beta_data_noSNP) > 0, data_noSNP)
beta_data_noSNP.css <- phyloseq_transform_css(beta_data_noSNP, log = F)
rel_abund <- transform_sample_counts(beta_data_noSNP.css, function(x) {x/sum(x)}*100)

beta_ra_group <- tax_glom(rel_abund, taxrank = "group")
beta_ra_class <- tax_glom(rel_abund, taxrank = "class")
beta_ra_class_melt_temp <- psmelt(beta_ra_class)

# ------------------------------------------------------------------------------------------
# COLLAPSE REDUNDANT COMBINED-RESISTANCE CLASSES SAFELY (regex-based)
# ------------------------------------------------------------------------------------------

beta_ra_class_melt_temp$class <- as.character(beta_ra_class_melt_temp$class)

# TRUE if class name contains:
# "biocide", "biocides", "metal", "multi-metal", "copper"
collapse_idx <- grepl("biocide",  beta_ra_class_melt_temp$class, ignore.case = TRUE) |
  grepl("metal",    beta_ra_class_melt_temp$class, ignore.case = TRUE) |
  grepl("copper",   beta_ra_class_melt_temp$class, ignore.case = TRUE)

beta_ra_class_melt_temp$class[collapse_idx] <- "Drug/Biocide/Metal resistance"

# --------- LOW ABUNDANCE FILTER NOW <0.1% -----------------------------------------------
medians <- ddply(beta_ra_class_melt_temp, ~class, function(x) c(median = median(x$Abundance)))
remainder <- medians[medians$median <= 0.1,]$class
beta_ra_class_melt_temp[beta_ra_class_melt_temp$class %in% remainder, ]$class <- 
  "low abundance class (<0.1%)"
# ------------------------------------------------------------------------------------------

factor_by_abund <- beta_ra_class_melt_temp %>%
  dplyr::group_by(class) %>%
  dplyr::summarize(median_class = median(Abundance)) %>%
  arrange(-median_class)

ps_AMR.dist <- vegdist(t(otu_table(beta_ra_group)), method = "bray")
ps_AMR.hclust <- hclust(ps_AMR.dist)
ps_AMR.dendro <- as.dendrogram(ps_AMR.hclust)
ps_AMR.dendro.data <- dendro_data(ps_AMR.dendro, type = "rectangle")
sample_names_ps_AMR <- ps_AMR.dendro.data$labels$label

AMR_mapfile <- as(sample_data(data_noSNP), "data.frame")
rownames(AMR_mapfile) <- gsub("-", ".", rownames(AMR_mapfile))
AMR_mapfile$Sample_ID <- rownames(AMR_mapfile)

ps_AMR.dendro.data$labels <- left_join(ps_AMR.dendro.data$labels,
                                       AMR_mapfile,
                                       by = c("label" = "Sample_ID"))

segment_data_ps_AMR <- with(
  segment(ps_AMR.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend)
)

gene_pos_table_ps_AMR <- with(ps_AMR.dendro.data$labels,
                              data.frame(y_center = x, gene = as.character(label), x = y,
                                         trt = as.character(Trial), order = as.character(Combined_order),
                                         height = 1))

sample_pos_table_ps_AMR <- data.frame(sample = sample_names_ps_AMR) %>%
  dplyr::mutate(x_center = (1:dplyr::n()), width = 1)

joined_ra_class_ps_class_melt <- beta_ra_class_melt_temp %>%
  left_join(gene_pos_table_ps_AMR, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_ps_AMR, by = c("Sample" = "sample"))

joined_ra_class_ps_class_melt$class <- factor(
  joined_ra_class_ps_class_melt$class,
  levels = as.character(factor_by_abund$class)
)

bar_totals <- joined_ra_class_ps_class_melt %>%
  dplyr::group_by(x_center) %>%
  dplyr::summarize(Total = sum(Abundance), .groups = "drop")

carto_bold <- c(
  "#06D6A0", "#FF61CC",  "#3333CC", "#00FFFF","#FF073A",
  "#5F0F40","#FFFF33", "#9B5DE5", "#6A1B9A",
  "#117733","#e0e0e0"
)

plt_rel_class_ps_class <- ggplot(joined_ra_class_ps_class_melt,
                                 aes(x = x_center, y = Abundance, fill = class)) +
  coord_flip() +
  geom_bar(stat = "identity", alpha = 0.95, width = 1, colour = NA) +
  geom_col(data = bar_totals, aes(x = x_center, y = Total),
           fill = NA, colour = "black", width = 1, linewidth = 0.5) +
  scale_fill_manual(values = carto_bold) +
  scale_x_continuous(breaks = sample_pos_table_ps_AMR$x_center,
                     labels = sample_pos_table_ps_AMR$sample,
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Relative Abundance", fill = "Antimicrobial Class") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black", face = "bold"),
    legend.key.size = unit(0.6, "cm"),
    legend.box.margin = margin(0, 10, 0, 10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.75),
    axis.text.x = element_text(size = 28, colour = "black"),
    axis.title.x = element_text(size = 28),
    axis.ticks.x = element_line(linewidth = 0.75),
    plot.margin = unit(c(-0.9, 0.5, 0.3, 0), "cm")
  )

plt_rel_class_ps_class

# (Dendrogram + save steps unchanged)

plt_rel_class_ps_class

legend_colors <- c("Start_No_Claim", "Midpoint_No_Claim", "End_No_Claim",
                   "Start_RWA", "Midpoint_RWA", "End_RWA",
                   "No Claim Beef", "RWA Beef")

c_colors <- c("Start_No_Claim"="#660000", "Midpoint_No_Claim"="#CC0066",
              "End_No_Claim"="#ff9999",
              "Start_RWA"="#000999", "Midpoint_RWA"="#006699", "End_RWA"="#66CCFF",
              "No Claim Beef"="#FF9900", "RWA Beef"="#117744")

plt_dendr_class_ps_class <- ggplot(segment_data_ps_AMR) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
               lineend = "round", linejoin = "round") +
  geom_point(data = gene_pos_table_ps_AMR,
             aes(x, y_center, color = order),
             size = 12, shape = 15, stroke = 0,
             position = position_nudge(x = 0.05)) +
  labs(x = "Ward's Distance", y = "", colour = "Sample Clusters") +
  scale_y_discrete(expand = c(0, 0, 0, 0)) +
  scale_x_reverse() +
  scale_color_manual(values = c_colors, breaks = legend_colors) +
  guides(color = guide_legend(
    nrow = 8,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  )) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black", face = "bold"),
    legend.key.size = unit(0.6, "cm"),
    legend.box.margin = margin(0, 10, 0, 10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.75),
    axis.text.x = element_text(size = 28, colour = "black"),
    axis.title.x = element_text(size = 28),
    axis.ticks.length = unit(0, "cm"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  )

plt_dendr_class_ps_class

combined_main <- plot_grid(
  plt_dendr_class_ps_class,
  plt_rel_class_ps_class,
  axis = "tb", align = "h",
  nrow = 1, rel_widths = c(0.5, 1.5)
)

combined_main

padded_combined <- ggdraw(combined_main) +
  theme(plot.margin = margin(20, 20, 20, 35, "pt"))

ggsave("Figures/AMR_RA_Final.png", padded_combined,
       width = 24, height = 34, units = "in", dpi = 600, bg = "white")

#############################################################################################
##############################   EXTRACT AND SAVE MATCHED LEGENDS (AMR)   ###################
#############################################################################################
library(cowplot)
library(grid)
library(gridExtra)


legend_dendr_amr <- get_legend(
  plt_dendr_class_ps_class +
    theme(
      legend.position = "right",
      legend.text  = element_text(size = 24, colour = "black"),
      legend.title = element_text(size = 26, face = "bold", colour = "black"),
      legend.key.size = unit(1, "cm")
    )
)

legend_rel_amr <- get_legend(
  plt_rel_class_ps_class +
    theme(
      legend.position = "right",
      legend.text  = element_text(size = 24, colour = "black"),
      legend.title = element_text(size = 26, face = "bold", colour = "black"),
      legend.key.size = unit(1, "cm")
    )
)

png("AAALegend_Dendrogram_AMR.png", width = 6, height = 8, units = "in", res = 300, bg = "white")
grid.newpage()
grid.draw(legend_dendr_amr)
dev.off()

png("AAALegend_RelAbundance_AMR.png", width = 6, height = 8, units = "in", res = 300, bg = "white")
grid.newpage
