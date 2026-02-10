#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

#Changed this to agglomerate earlier on in the analysis per PM/ED Oct 2024

#Make new phyloseq object
alpha_data_noSNP <- data_noSNP

#1-24-25 Why are you removing the Final Washout for alpha? 
#alpha_data_noSNP_noWashout.ps <- subset_samples(alpha_data_noSNP, Combined_order != "Final_washout")
#alpha_data_noSNP_noWashout.ps <- prune_taxa(taxa_sums(alpha_data_noSNP_noWashout.ps) > 0, alpha_data_noSNP_noWashout.ps)

#Agglomerate for alpha diversity, NOT normalized
alpha_data_noSNP.ps <- alpha_data_noSNP
alpha_ra_group <- tax_glom(alpha_data_noSNP.ps, taxrank = "group")
alpha_ra_group # 615 groups
alpha_ra_group_melt <- psmelt(alpha_ra_group)


AMR_alpha_div <- estimate_richness(alpha_ra_group, measures = c("Observed","Shannon","Simpson","InvSimpson"))
AMR_alpha_div
AMR_alpha_div.df <- as(sample_data(alpha_ra_group), "data.frame")
AMR_alpha_div_meta <- cbind(AMR_alpha_div, AMR_alpha_div.df)


#Make custom colors
cb_colors <- c("#0072B2","#882255")

alpha_richness_sampletype<- ggplot(AMR_alpha_div_meta, aes(x= Sample_Type, y= Observed, fill = Sample_Type, colour = Sample_Type)) +
  theme_bw() + 
  labs(y= "Observed ARGs") +
  geom_violin(alpha=0.4) +
  geom_point() +
  scale_fill_manual(values = cb_colors)+
  scale_color_manual(values = cb_colors)+
  theme(legend.position = "bottom", 
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
alpha_richness_sampletype

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
Sample_Type.richness.pw <- pairwise.wilcox.test(AMR_alpha_div_meta$Observed, AMR_alpha_div_meta$Sample_Type, p.adjust.method = "BH")
Sample_Type.richness.pw # p-values in the matrix


alpha_shannons_sampletype <- ggplot(AMR_alpha_div_meta, aes(x= Sample_Type, y= Shannon, fill = Sample_Type, colour = Sample_Type)) +
  theme_bw() + 
  labs(y= "Shannon's diversity") +
  geom_violin(alpha=0.4) +
  geom_point() +
  scale_fill_manual(values = cb_colors)+
  scale_color_manual(values = cb_colors)+
  theme(legend.position = "bottom",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
alpha_shannons_sampletype

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
Sample_Type.shannon.pw <- pairwise.wilcox.test(AMR_alpha_div_meta$Shannon, AMR_alpha_div_meta$Sample_Type, p.adjust.method = "BH")
Sample_Type.shannon.pw # p-values in the matrix

#Combine plots into panel for publication


#1-24-25 Added custom colors, now the extract legend is not working
# Extract the legend from the alpha_richness_sampletype plot
legend <- cowplot::get_legend(alpha_richness_sampletype)

# Combine plots into a single panel without the legend
Sample_Type_combined <- plot_grid(alpha_richness_sampletype + theme(legend.position = "bottom"), 
                                  alpha_shannons_sampletype + theme(legend.position = "bottom"), 
                                  ncol = 2)

# Combine the title, the combined plots, and the legend
Sample_Type_final <- plot_grid(
  ggdraw() + draw_label("Alpha Diversity by Sample Type", size = 18), 
  Sample_Type_combined, 
  legend, 
  ncol = 1, 
  rel_heights = c(0.1, 1, 0.1))
Sample_Type_final


##
#### Now, only fecal samples #####
##

# Extract fecal samples
alpha_div_fecal <- AMR_alpha_div_meta[AMR_alpha_div_meta$Sample_Type == "Fecal", ]
alpha_data_noSNP_noWashout <- alpha_div_fecal[alpha_div_fecal$Combined_order != "Final_washout", ]
alpha_div_fecal <- prune_taxa(taxa_sums(alpha_data_noSNP_noWashout.ps) > 0, alpha_data_noSNP_noWashout.ps)
if (exists("alpha_ra_group")) {alpha_data_noSNP_noWashout.ps <- prune_taxa(taxa_sums(alpha_ra_group) > 0, 
    subset_samples(alpha_ra_group, Sample_Type == "Fecal" & Combined_order != "Final_washout"))}

# Change ordre of "Combined_order"
alpha_div_fecal$Combined_order <- factor(alpha_div_fecal$Combined_order, levels = c("Start_CONV","Midpoint_CONV","End_CONV", "Start_RWA",
                                                                                    "Midpoint_RWA","End_RWA","Final_washout"))

#Define Custom Colors do you can split RWA from Conv
#Undecided on what to do with Final Washout at this point
custom_colors <- c(
  "Start_RWA" = "#000066",   
  "Midpoint_RWA" = "#0033ff",     
  "End_RWA" = "#6699ff",     
  "Start_CONV" = "#990000", 
  "Midpoint_CONV" = "#ff0000",  
  "End_CONV" = "#ff3333",
  "Final_washout" = "#00ffcc")

#Do these all one color for each group instead of timepoints and group
custom_colors2 <- c(
  "Start_RWA" = "#0033ff",   
  "Midpoint_RWA" = "#0033ff",     
  "End_RWA" = "#0033ff",     
  "Start_CONV" = "#ff0000", 
  "Midpoint_CONV" = "#ff0000",  
  "End_CONV" = "#ff0000",
  "Final_washout" = "#999999")




alpha_richness_fecal <- ggplot(alpha_div_fecal, aes(x= Combined_order, y= Observed, fill = Combined_order, colour = Combined_order)) +
  theme_bw() + 
  labs(y= "Observed ARGs") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  scale_fill_manual(values = custom_colors2) +
  scale_colour_manual(values = custom_colors2) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
alpha_richness_fecal

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
Combined_order.richness.pw <- pairwise.wilcox.test(alpha_div_fecal$Observed, alpha_div_fecal$Combined_order, p.adjust.method = "BH")
Combined_order.richness.pw # p-values in the matrix


## Shannon's #
alpha_shannons_fecal <- ggplot(alpha_div_fecal, aes(x= Combined_order, y= Shannon, fill = Combined_order, colour = Combined_order)) +
  theme_bw() + 
  labs(y= "Shannon's diversity") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  scale_fill_manual(values = custom_colors2) +
  scale_colour_manual(values = custom_colors2) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
alpha_shannons_fecal

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
Combined_order.shannon.pw <- pairwise.wilcox.test(alpha_div_fecal$Shannon, alpha_div_fecal$Combined_order, p.adjust.method = "BH")
Combined_order.shannon.pw # p-values in the matrix

#Combine plots into panel for publication
# Extract the legend from one plot
legend <- get_legend(alpha_richness_fecal)

# Combine plots into a single panel without the legend
fecal_combined <- plot_grid(alpha_richness_fecal + theme(legend.position = "none"), 
                                  alpha_shannons_fecal + theme(legend.position = "none"), 
                                  ncol = 2)

# Combine the title, the combined plots, and the legend
fecal_final <- plot_grid(
  ggdraw() + draw_label("Alpha Diversity of Fecal Samples by Time", size = 18), 
  fecal_combined, 
  legend, 
  ncol = 1, 
  rel_heights = c(0.1, 1, 0.1))
fecal_final
