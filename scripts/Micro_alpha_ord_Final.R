#############################################################################################
##############################   MICROBIOME ORDINATION (FINAL CLEAN + BIG LEGEND, UNBOLDED AXES)   ##############################
#############################################################################################

# Remove Final_Washout
beta_micro.ps <- subset_samples(data_micro.ps, Combined_order_ext != "Final_Washout")
beta_micro.ps <- prune_taxa(taxa_sums(beta_micro.ps) > 0, beta_micro.ps)

# Normalize with CSS
beta_micro.ps.css <- phyloseq_transform_css(beta_micro.ps, log = F)


# Relative abundance
rel_abund_micro <- transform_sample_counts(beta_micro.ps.css, function(x) {x/sum(x)}*100)

# Set Combined_order_ext factor levels (Beef first, then Fecal)
sample_data(rel_abund_micro)$Combined_order_ext <- factor(
  as.character(sample_data(rel_abund_micro)$Combined_order_ext),
  levels = c("No Claim Beef","RWA Beef",
             "Start No Claim Fecal","Midpoint No Claim Fecal","End No Claim Fecal",
             "Start RWA Fecal","Midpoint RWA Fecal","End RWA Fecal")
)

# ----- Define level sets -----
fecal_levels <- c("Start No Claim Fecal","Midpoint No Claim Fecal","End No Claim Fecal",
                  "Start RWA Fecal","Midpoint RWA Fecal","End RWA Fecal")
beef_levels  <- c("No Claim Beef","RWA Beef")

# ----- Define matching color palettes -----
fecal_colors <- c(
  "Start No Claim Fecal"    = "#660066",
  "Midpoint No Claim Fecal" = "#990099",
  "End No Claim Fecal"      = "#cc00ff",
  "Start RWA Fecal"         = "#006699",
  "Midpoint RWA Fecal"      = "#00ccdd",
  "End RWA Fecal"           = "#99ffff"
)

beef_colors <- c(
  "No Claim Beef" = "#cc0033",
  "RWA Beef"      = "#000099"
)

# ----- Subset data -----
rel_abund_micro_fecal <- subset_samples(rel_abund_micro, Combined_order_ext %in% fecal_levels)
rel_abund_micro_fecal <- prune_taxa(taxa_sums(rel_abund_micro_fecal) > 0, rel_abund_micro_fecal)

rel_abund_micro_beef  <- subset_samples(rel_abund_micro, Combined_order_ext %in% beef_levels)
rel_abund_micro_beef  <- prune_taxa(taxa_sums(rel_abund_micro_beef) > 0, rel_abund_micro_beef)

# ----- Ordinations -----
data.ord.micro_fecal <- metaMDS(t(otu_table(rel_abund_micro_fecal)), distance = "bray",
                                try = 10, trymax = 999, autotransform = F)
data.ord.micro_beef  <- metaMDS(t(otu_table(rel_abund_micro_beef)),  distance = "bray",
                                try = 10, trymax = 999, autotransform = F)

# Extract ordination data
ord_points_micro_fecal <- plot_ordination(rel_abund_micro_fecal, data.ord.micro_fecal, color = "Combined_order_ext")
ord_points_micro_beef  <- plot_ordination(rel_abund_micro_beef,  data.ord.micro_beef,  color = "Combined_order_ext")

# Flip NMDS1 to match AMR
ord_points_micro_fecal$data$NMDS1 <- -ord_points_micro_fecal$data$NMDS1
ord_points_micro_beef$data$NMDS1  <- -ord_points_micro_beef$data$NMDS1

# ----- FECAL-ONLY PLOT -----
plt_ord_micro_fecal <- ggplot() +
  geom_point(
    data = ord_points_micro_fecal$data,
    aes(x = NMDS1, y = NMDS2, fill = Combined_order_ext),
    shape = 21, colour = "black", size = 3.5, stroke = 0.8
  ) +
  stat_ellipse(
    data = ord_points_micro_fecal$data,
    aes(x = NMDS1, y = NMDS2, colour = Combined_order_ext),
    lty = 1, size = 1
  ) +
  scale_fill_manual(values = fecal_colors, breaks = fecal_levels, drop = FALSE) +
  scale_color_manual(values = fecal_colors, breaks = fecal_levels, drop = FALSE) +
  guides(
    fill = guide_legend(
      ncol = 2,
      override.aes = list(shape = 21, size = 9, colour = "black"),  # bigger symbols
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
    legend.key.size = unit(1.5, "cm"),            # bigger legend keys
    legend.text = element_text(size = 24, colour = "black"),  # bigger legend text
    legend.title = element_blank(),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 16),  # keep same
    axis.ticks = element_line(colour = "black", size = 0.8),
    panel.border = element_rect(colour = "black", size = 1.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ----- BEEF-ONLY PLOT -----
plt_ord_micro_beef <- ggplot() +
  geom_point(
    data = ord_points_micro_beef$data,
    aes(x = NMDS1, y = NMDS2, fill = Combined_order_ext),
    shape = 21, colour = "black", size = 3.5, stroke = 0.8
  ) +
  stat_ellipse(
    data = ord_points_micro_beef$data,
    aes(x = NMDS1, y = NMDS2, colour = Combined_order_ext),
    lty = 1, size = 1
  ) +
  scale_fill_manual(values = beef_colors, breaks = beef_levels, drop = FALSE) +
  scale_color_manual(values = beef_colors, breaks = beef_levels, drop = FALSE) +
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
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 16),  # keep same
    axis.ticks = element_line(colour = "black", size = 0.8),
    panel.border = element_rect(colour = "black", size = 1.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Show plots
plt_ord_micro_fecal
plt_ord_micro_beef


############################################################################################
####################################   MICROBIOME PERMANOVA (FECAL ONLY)   ####################################
############################################################################################

# --- Subset to FECAL ONLY ---
beta_micro_fecal <- subset_samples(rel_abund_micro, sample_type == "Fecal")
beta_micro_fecal <- prune_taxa(taxa_sums(beta_micro_fecal) > 0, beta_micro_fecal)

# --- Distance matrix and metadata (FECAL ONLY) ---
micro_dist_fecal <- phyloseq::distance(beta_micro_fecal, method = "bray")
micro_meta_fecal <- as(sample_data(beta_micro_fecal), "data.frame")

# --- Clean up missing values (remove washout or NA treatment samples) ---
micro_meta_fecal <- micro_meta_fecal[
  complete.cases(micro_meta_fecal[, c("pool_timepoint", "trt", "treatment_order")]),
]

# --- Convert distance object to full matrix and match sample order ---
micro_dist_fecal <- as.matrix(micro_dist_fecal)
micro_dist_fecal <- micro_dist_fecal[rownames(micro_meta_fecal), rownames(micro_meta_fecal)]

# --- Run PERMANOVA with participant-level strata ---
perm_micro_fecal <- adonis2(
  micro_dist_fecal ~ pool_timepoint * trt * treatment_order,
  data = micro_meta_fecal,
  permutations = 9999,
  strata = micro_meta_fecal$participant_id
)

# --- View results ---
print(perm_micro_fecal)
