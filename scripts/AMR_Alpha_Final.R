#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################

# Make new phyloseq object
alpha_data_noSNP <- data_noSNP


# Remove Final Washout but keep ALL sample types (Fecal + Beef)
alpha_data_noSNP_noWashout.ps <- subset_samples(alpha_data_noSNP, Combined_order != "Final_washout")
alpha_data_noSNP_noWashout.ps <- prune_taxa(taxa_sums(alpha_data_noSNP_noWashout.ps) > 0, alpha_data_noSNP_noWashout.ps)

# Agglomerate for alpha diversity, NOT normalized
alpha_ra_group <- tax_glom(alpha_data_noSNP_noWashout.ps, taxrank = "group")

# Calculate richness/diversity
AMR_alpha_div <- estimate_richness(alpha_ra_group, measures = c("Observed","Shannon","Simpson","InvSimpson"))
AMR_alpha_div.df <- as(sample_data(alpha_ra_group), "data.frame")
AMR_alpha_div_meta <- cbind(AMR_alpha_div, AMR_alpha_div.df)

# Relabel fecal groups with "Fecal"
AMR_alpha_div_meta$Combined_order <- dplyr::recode(as.character(AMR_alpha_div_meta$Combined_order),
                                                   "Start_No_Claim"    = "Start_No_Claim Fecal",
                                                   "Midpoint_No_Claim" = "Midpoint_No_Claim Fecal",
                                                   "End_No_Claim"      = "End_No_Claim Fecal",
                                                   "Start_RWA"         = "Start_RWA Fecal",
                                                   "Midpoint_RWA"      = "Midpoint_RWA Fecal",
                                                   "End_RWA"           = "End_RWA Fecal",
                                                   "No Claim Beef"     = "No Claim Beef",
                                                   "RWA Beef"          = "RWA Beef"
)

# Order factor levels (Beef FIRST)
AMR_alpha_div_meta$Combined_order <- factor(
  AMR_alpha_div_meta$Combined_order,
  levels = c("No Claim Beef","RWA Beef",
             "Start_No_Claim Fecal","Midpoint_No_Claim Fecal","End_No_Claim Fecal",
             "Start_RWA Fecal","Midpoint_RWA Fecal","End_RWA Fecal")
)

# Colors
custom_colors2 <- c(
  "No Claim Beef" = "#FFD700", "RWA Beef" = "#117744",
  "Start_RWA Fecal" = "#0033ff", "Midpoint_RWA Fecal" = "#0033ff", "End_RWA Fecal" = "#0033ff",
  "Start_No_Claim Fecal" = "#ff0000", "Midpoint_No_Claim Fecal" = "#ff0000", "End_No_Claim Fecal" = "#ff0000"
)



##### Observed richness plot ####
alpha_richness_all <- ggplot(AMR_alpha_div_meta,
                             aes(x = Combined_order, y = Observed,
                                 fill = Combined_order)) +
  theme_bw() + 
  labs(y = "Observed ARGs") +
  geom_boxplot(alpha = 0.4, aes(colour = Combined_order), show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.15), size = 2,
             shape = 21, stroke = 0.7, colour = "black") +
  scale_fill_manual(values = custom_colors2) +
  scale_colour_manual(values = custom_colors2) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_rect(fill = "black", size = 1.0),
    strip.text = element_text(size = 24, colour = "white"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),  # <-- add this
    axis.ticks.x = element_line(colour = "black", size = 0.7),                       # <-- and this
    axis.ticks.y = element_line(colour = "black", size = 0.7),
    plot.title = element_text(size = 28),
    panel.border = element_rect(colour = "black", size = 1.0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Stats
#Observed.pw <- pairwise.wilcox.test(AMR_alpha_div_meta$Observed, AMR_alpha_div_meta$Combined_order,  p.adjust.method = "BH")
#Observed.pw


#### Shannon diversity plot ####
#### Shannon diversity plot ####
alpha_shannons_all <- ggplot(AMR_alpha_div_meta,
                             aes(x = Combined_order, y = Shannon,
                                 fill = Combined_order)) +
  theme_bw() + 
  labs(y = "Shannon's diversity") +
  geom_boxplot(alpha = 0.4, aes(colour = Combined_order), show.legend = FALSE) +  # no legend from boxplot
  geom_point(position = position_jitter(width = 0.15), size = 2,
             shape = 21, stroke = 0.7, colour = "black") +  # outlined dots
  scale_fill_manual(values = custom_colors2) +
  scale_colour_manual(values = custom_colors2) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_rect(fill = "black", size = 1.0),
    strip.text = element_text(size = 24, colour = "white"),
    axis.text = element_text(size = 14, colour = "black"),        # unified axis text
    axis.title = element_text(size = 18),                         # unified axis title
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),  # rotated x labels
    axis.ticks.x = element_line(colour = "black", size = 0.7),                       # show x ticks
    axis.ticks.y = element_line(colour = "black", size = 0.7),
    plot.title = element_text(size = 28),
    panel.border = element_rect(colour = "black", size = 1.0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Stats
#Shannon.pw <- pairwise.wilcox.test(AMR_alpha_div_meta$Shannon, AMR_alpha_div_meta$Combined_order,               p.adjust.method = "BH")
#Shannon.pw









#Combine into one panel
alpha_ord_panel <- plot_grid(
  alpha_richness_all,
  alpha_shannons_all,
  plt_ord_by_combined,
  ncol = 3,
  rel_widths = c(1, 1, 1.6),  # ord gets more space to balance width
  label_size = 20,
  label_fontface = "bold"
)

alpha_ord_panel

ggsave("AMR_Alpha_Ord_Final.png", alpha_ord_panel,
       width = 16, height = 8, dpi = 300, bg = "white")

# Save final
ggsave("AMR_Alpha_Final.png", alpha_div_all_panel,
       width = 16, height = 8, dpi = 300, bg = "white")
