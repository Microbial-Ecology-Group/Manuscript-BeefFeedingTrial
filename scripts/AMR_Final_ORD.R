#############################################################################################
##############################   AMR ORDINATION (GENE GROUP LEVEL, FINAL CLEAN + BIG LEGEND, UNBOLDED AXES)   ##############################
#############################################################################################

# Make new phyloseq object
beta_data_noSNP <- data_noSNP
beta_data_noSNP <- prune_taxa(taxa_sums(beta_data_noSNP) > 0, beta_data_noSNP)

# Normalize with CSS
beta_data_noSNP.css <- phyloseq_transform_css(beta_data_noSNP, log = F)

# Relative abundance
rel_abund_amr <- transform_sample_counts(beta_data_noSNP.css, function(x) {x/sum(x)}*100)

# Agglomerate at the gene GROUP level
beta_ra_group <- tax_glom(rel_abund_amr, taxrank = "group")

# Relabel fecal groups with "Fecal" suffix
sample_data(beta_ra_group)$Combined_order <- dplyr::recode(
  as.character(sample_data(beta_ra_group)$Combined_order),
  "Start_No_Claim"    = "Start No Claim Fecal",
  "Midpoint_No_Claim" = "Midpoint No Claim Fecal",
  "End_No_Claim"      = "End No Claim Fecal",
  "Start_RWA"         = "Start RWA Fecal",
  "Midpoint_RWA"      = "Midpoint RWA Fecal",
  "End_RWA"           = "End RWA Fecal",
  "No Claim Beef"     = "No Claim Beef",
  "RWA Beef"          = "RWA Beef"
)

# Set factor levels (Beef first, then Fecal)
sample_data(beta_ra_group)$Combined_order <- factor(
  as.character(sample_data(beta_ra_group)$Combined_order),
  levels = c("No Claim Beef","RWA Beef",
             "Start No Claim Fecal","Midpoint No Claim Fecal","End No Claim Fecal",
             "Start RWA Fecal","Midpoint RWA Fecal","End RWA Fecal")
)

# Define colors
amr_colors <- c(
  "No Claim Beef"          = "#FFD700",
  "RWA Beef"               = "#117744",
  "Start No Claim Fecal"   = "#660000",
  "Midpoint No Claim Fecal"= "#CC0066",
  "End No Claim Fecal"     = "#FF9999",
  "Start RWA Fecal"        = "#000999",
  "Midpoint RWA Fecal"     = "#006699",
  "End RWA Fecal"          = "#66CCFF"
)

# ----- Subset for beef and fecal -----
fecal_levels <- c("Start No Claim Fecal","Midpoint No Claim Fecal","End No Claim Fecal",
                  "Start RWA Fecal","Midpoint RWA Fecal","End RWA Fecal")
beef_levels  <- c("No Claim Beef","RWA Beef")

beta_ra_group_fecal <- subset_samples(beta_ra_group, Combined_order %in% fecal_levels)
beta_ra_group_fecal <- prune_taxa(taxa_sums(beta_ra_group_fecal) > 0, beta_ra_group_fecal)

beta_ra_group_beef <- subset_samples(beta_ra_group, Combined_order %in% beef_levels)
beta_ra_group_beef <- prune_taxa(taxa_sums(beta_ra_group_beef) > 0, beta_ra_group_beef)

# ----- Ordinations -----
ord_amr_fecal <- metaMDS(t(otu_table(beta_ra_group_fecal)), distance = "bray",
                         try = 10, trymax = 999, autotransform = F)
ord_amr_beef  <- metaMDS(t(otu_table(beta_ra_group_beef)),  distance = "bray",
                         try = 10, trymax = 999, autotransform = F)

ord_points_amr_fecal <- plot_ordination(beta_ra_group_fecal, ord_amr_fecal, color = "Combined_order")
ord_points_amr_beef  <- plot_ordination(beta_ra_group_beef,  ord_amr_beef,  color = "Combined_order")

# Flip NMDS1 to match microbiome orientation
ord_points_amr_fecal$data$NMDS1 <- -ord_points_amr_fecal$data$NMDS1
ord_points_amr_beef$data$NMDS1  <- -ord_points_amr_beef$data$NMDS1

# ---- MATCH AXIS LIMITS WITHIN AMR PLOTS ----
x_limits_amr <- range(c(ord_points_amr_fecal$data$NMDS1, ord_points_amr_beef$data$NMDS1))
y_limits_amr <- range(c(ord_points_amr_fecal$data$NMDS2, ord_points_amr_beef$data$NMDS2))

# ----- FECAL PLOT -----
plt_ord_amr_fecal <- ggplot() +
  geom_point(
    data = ord_points_amr_fecal$data,
    aes(x = NMDS1, y = NMDS2, fill = Combined_order),
    shape = 21, colour = "black", size = 3.5, stroke = 0.8
  ) +
  stat_ellipse(
    data = ord_points_amr_fecal$data,
    aes(x = NMDS1, y = NMDS2, colour = Combined_order),
    lty = 1, size = 1
  ) +
  scale_fill_manual(values = amr_colors, breaks = fecal_levels, drop = FALSE) +
  scale_color_manual(values = amr_colors, breaks = fecal_levels, drop = FALSE) +
  guides(
    fill = guide_legend(
      ncol = 2,
      override.aes = list(shape = 21, size = 9, colour = "black"),
      title.theme = element_text(size = 24, face = "bold"),
      label.theme = element_text(size = 24)
    ),
    color = "none"
  ) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 10, b = 10),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, colour = "black"),
    legend.title = element_blank(),
    plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 16),  # unbolded NMDS
    axis.ticks = element_line(colour = "black", size = 0.8),
    panel.border = element_rect(colour = "black", size = 1.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ----- BEEF PLOT -----
plt_ord_amr_beef <- ggplot() +
  geom_point(
    data = ord_points_amr_beef$data,
    aes(x = NMDS1, y = NMDS2, fill = Combined_order),
    shape = 21, colour = "black", size = 3.5, stroke = 0.8
  ) +
  stat_ellipse(
    data = ord_points_amr_beef$data,
    aes(x = NMDS1, y = NMDS2, colour = Combined_order),
    lty = 1, size = 1
  ) +
  scale_fill_manual(values = amr_colors, breaks = beef_levels, drop = FALSE) +
  scale_color_manual(values = amr_colors, breaks = beef_levels, drop = FALSE) +
  guides(
    fill = guide_legend(
      ncol = 2,
      override.aes = list(shape = 21, size = 9, colour = "black"),
      title.theme = element_text(size = 24, face = "bold"),
      label.theme = element_text(size = 24)
    ),
    color = "none"
  ) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 10, b = 10),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, colour = "black"),
    legend.title = element_blank(),
    plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 16),  # unbolded NMDS
    axis.ticks = element_line(colour = "black", size = 0.8),
    panel.border = element_rect(colour = "black", size = 1.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Show plots
plt_ord_amr_fecal
plt_ord_amr_beef



# === Save ordination ===
#ggsave("AMR_Ord_Final.png", plt_ord_by_combined,
       width = 8, height = 6, dpi = 300, bg = "white")
