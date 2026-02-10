#############################################################################################
##############################   MICROBIOME RELATIVE ABUNDANCE + DENDROGRAM (MATCHED TO AMR)  ##############################
#############################################################################################

#############################################################################################
####################################   CSS TRANSFORM    ###################################
#############################################################################################

# INCLUDE BOTH FECAL + BEEF SAMPLES
# SUBSET: remove Final_washout samples only
micro_subset <- subset_samples(data_micro.ps, Combined_order_ext != "Final_Washout")
micro.ps <- micro_subset

#############################################################################################
##############################   RELATIVE ABUNDANCE + DENDROGRAM   #########################
#############################################################################################

# CSS Transform
micro.ps.css <- phyloseq_transform_css(micro.ps, log = FALSE)
micro.ps.css.df <- as(sample_data(micro.ps.css), "data.frame")

# Remove unclassifieds
to_remove <- apply(tax_table(micro.ps.css), 1, function(x) any(x %in% c("unclassified Unassigned", "unclassified Bacteria")))
micro.ps.css.filtered <- prune_taxa(!to_remove, micro.ps.css)

# Relative abundance
rel_abund <- transform_sample_counts(micro.ps.css.filtered, function(x) {x/sum(x)} * 100)

# Agglomerate at Order level
ra_order <- tax_glom(rel_abund, taxrank = "Order")
ra_order_melt_temp <- psmelt(ra_order)

# Collapse low abundance taxa (<1%)
ra_order_melt_temp$Order <- as.character(ra_order_melt_temp$Order)
medians <- ddply(ra_order_melt_temp, ~Order, function(x) c(median = median(x$Abundance, na.rm = TRUE)))
remainder <- medians[medians$median <= 0.01, ]$Order
ra_order_melt_temp$Order[ra_order_melt_temp$Order %in% remainder] <- "Low abundance Orders (<1%)"

# Reorder class levels by median abundance
factor_by_abund <- ra_order_melt_temp %>%
  dplyr::group_by(Order) %>%
  dplyr::summarize(median_order = median(Abundance)) %>%
  arrange(-median_order)

#############################################################################################
###################################   CLUSTERING   #########################################
#############################################################################################

ps_micro.dist <- vegdist(t(otu_table(micro.ps.css.filtered)), method = "bray")
ps_micro.hclust <- hclust(ps_micro.dist)
ps_micro.dendro <- as.dendrogram(ps_micro.hclust)
ps_micro.dendro.data <- dendro_data(ps_micro.dendro, type = "rectangle")

sample_names_ps_micro <- ps_micro.dendro.data$labels$label

# Add metadata
mapfile_micro <- as(sample_data(micro.ps), "data.frame")
mapfile_micro$Sample_ID <- rownames(mapfile_micro)

ps_micro.dendro.data$labels <- left_join(ps_micro.dendro.data$labels,
                                         mapfile_micro,
                                         by = c("label" = "Sample_ID"))

segment_data_ps_micro <- with(segment(ps_micro.dendro.data),
                              data.frame(x = y, y = x, xend = yend, yend = xend))

pos_table_ps_micro <- with(ps_micro.dendro.data$labels,
                           data.frame(y_center = x,
                                      sample = as.character(label),
                                      x = y,
                                      Combined_order_ext = as.character(Combined_order_ext),
                                      height = 1))

sample_pos_table_ps_micro <- data.frame(sample = sample_names_ps_micro) %>%
  dplyr::mutate(x_center = (1:dplyr::n()), width = 1)

#############################################################################################
##############################   RELATIVE ABUNDANCE PLOT   #################################
#############################################################################################

joined_ra_order_melt <- ra_order_melt_temp %>%
  left_join(pos_table_ps_micro,  by = c("Sample" = "sample")) %>%
  left_join(sample_pos_table_ps_micro, by = c("Sample" = "sample"))

joined_ra_order_melt$Order <- as.character(joined_ra_order_melt$Order)

factor_by_abund <- joined_ra_order_melt %>%
  dplyr::group_by(Order) %>%
  dplyr::summarize(median_order = median(Abundance)) %>%
  arrange(-median_order)

joined_ra_order_melt$Order <- factor(joined_ra_order_melt$Order,
                                     levels = as.character(factor_by_abund$Order))

bar_totals_micro <- joined_ra_order_melt %>%
  dplyr::group_by(x_center) %>%
  dplyr::summarize(Total = sum(Abundance), .groups = "drop")




MM_colors <- c("#77CCCC","#E73f74","#11a579","#ccccff","#0047ab","#117777",
               "#cf1c90","#771155","#00ffff","#764E9F","#ffea00","#0a1172",
               "#006400","#66ff66","#ff5e00","#ff1493","#e0e0e0","#ADFF2F",
               "#8A2BE2","#00BFFF","#7B68EE")

# --- Relative Abundance Plot ---
plt_rel_ra_order <- ggplot(joined_ra_order_melt,
                           aes(x = x_center, y = Abundance, fill = Order)) +
  coord_flip() +
  geom_bar(stat = "identity", alpha = 0.95, width = 1, colour = NA) +
  geom_col(data = bar_totals_micro, aes(x = x_center, y = Total),
           fill = NA, colour = "black", width = 1, linewidth = 0.5) +
  scale_fill_manual(values = MM_colors, name = "Order") +
  scale_x_continuous(breaks = sample_pos_table_ps_micro$x_center,
                     labels = sample_pos_table_ps_micro$sample,
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "",
    y = "Relative Abundance"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black", face = "bold"),
    legend.key.size = unit(0.6, "cm"),
    legend.box.margin = margin(0, 10, 0, 10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.75),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.75, lineend = "square", colour = "black"),
    axis.title.x = element_text(size = 28),
    axis.text.x = element_text(size = 28, colour = "black"),
    plot.margin = unit(c(-0.9, 0.5, 0.3, 0), "cm")
  )

plt_rel_ra_order


#############################################################################################
##############################   DENDROGRAM PLOT   #########################################
#############################################################################################

legend_colors_micro <- c("Start No Claim Fecal", "Midpoint No Claim Fecal", "End No Claim Fecal",
                         "Start RWA Fecal", "Midpoint RWA Fecal", "End RWA Fecal",
                         "No Claim Beef", "RWA Beef")

c_colors_micro <- c("Start No Claim Fecal" = "#660000",
                    "Midpoint No Claim Fecal" = "#CC0066",
                    "End No Claim Fecal" = "#ff9999",
                    "Start RWA Fecal" = "#000999",
                    "Midpoint RWA Fecal" = "#006699",
                    "End RWA Fecal" = "#66CCFF",
                    "No Claim Beef" = "#ff9900",
                    "RWA Beef" = "#117744")

plt_dendr_order <- ggplot(segment_data_ps_micro) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
               lineend = "round", linejoin = "round") +
  geom_point(data = pos_table_ps_micro,
             aes(x, y_center, color = Combined_order_ext),
             size = 12, shape = 15, stroke = 0,
             position = position_nudge(x = 0.05)) +
  labs(
    x = "Ward's Distance",
    y = "",
    colour = "Sample Clusters"
  ) +
  scale_y_discrete(expand = c(0, 0, 0, 0)) +
  scale_x_reverse() +
  scale_color_manual(values = c_colors_micro, breaks = legend_colors_micro) +
  guides(color = guide_legend(
    nrow = 8,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  )) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 12, colour = "black", face = "bold"),
    legend.key.size = unit(0.6, "cm"),
    legend.box.margin = margin(0, 10, 0, 10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.75),
    axis.title.x = element_text(size = 28),
    axis.text.x = element_text(size = 28, colour = "black"),
    axis.ticks.length = unit(0, "cm"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  )

plt_dendr_order


#############################################################################################
##############################   COMBINE + SAVE (MATCHED)   ################################
#############################################################################################

combined_micro <- plot_grid(
  plt_dendr_order,
  plt_rel_ra_order,
  axis = "tb",
  align = "h",
  nrow = 1,
  rel_widths = c(0.5, 1.5)
)
combined_micro

padded_combined_micro <- ggdraw(combined_micro) +
  theme(plot.margin = margin(20, 20, 20, 35, "pt"))

ggsave("Microbiome_RA_Final.png", padded_combined_micro,
       width = 28, height = 34, units = "in", dpi = 300, bg = "white")



#Legends only
library(cowplot)
library(grid)

legend_dendr <- get_legend(
  plt_dendr_order +
    theme(
      legend.text  = element_text(size = 24, colour = "black"),
      legend.title = element_text(size = 26, face = "bold", colour = "black"),
      legend.key.size = unit(1, "cm")
    )
)

legend_rel <- get_legend(
  plt_rel_ra_order +
    theme(
      legend.position = "right",
      legend.text  = element_text(size = 24, colour = "black"),
      legend.title = element_text(size = 26, face = "bold", colour = "black"),
      legend.key.size = unit(1, "cm")
    )
)

library(cowplot)
library(grid)

legend_dendr <- get_legend(
  plt_dendr_order +
    theme(
      legend.text  = element_text(size = 24, colour = "black"),
      legend.title = element_text(size = 26, face = "bold", colour = "black"),
      legend.key.size = unit(1, "cm")
    )
)

legend_rel <- get_legend(
  plt_rel_ra_order +
    theme(
      legend.position = "right",
      legend.text  = element_text(size = 24, colour = "black"),
      legend.title = element_text(size = 26, face = "bold", colour = "black"),
      legend.key.size = unit(1, "cm")
    )
)

# Save dendrogram legend
png("Legend_Dendrogram_Microbiome.png", width = 6, height = 8, units = "in", res = 300, bg = "white")
grid.newpage()
grid.draw(legend_dendr)
dev.off()

# Save relative abundance legend
png("Legend_RelAbundance_Microbiome.png", width = 6, height = 8, units = "in", res = 300, bg = "white")
grid.newpage()
grid.draw(legend_rel)
dev.off()

