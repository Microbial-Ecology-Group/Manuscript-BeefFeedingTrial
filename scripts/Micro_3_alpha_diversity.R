#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################
#Make custom colors
cb_colors <- c("#0072B2","#882255")

micro.alpha_div <- estimate_richness(data_micro.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))
micro.alpha_div

micro.alpha_div.df <- as(sample_data(data_micro.ps), "data.frame")
micro.alpha_div_meta <- cbind(micro.alpha_div, micro.alpha_div.df)

colnames(micro.alpha_div_meta)[colnames(micro.alpha_div_meta) == "sample_type"] <- "Sample_type"

micro_alpha_richness_sample_type <- ggplot(micro.alpha_div_meta, aes(x= Sample_type, y= Observed, fill = Sample_type, colour = Sample_type)) +
  theme_bw() + 
  labs(y= "Observed ASVs") +
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
micro_alpha_richness_sample_type

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.Combined_order.richness.pw <- pairwise.wilcox.test(micro.alpha_div_meta$Observed, micro.alpha_div_meta$Combined_order, p.adjust.method = "BH")
micro.Combined_order.richness.pw # p-values in the matrix


micro_alpha_shannons_sample_type <- ggplot(micro.alpha_div_meta, aes(x= Sample_type, y= Shannon, fill = Sample_type, colour = Sample_type)) +
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
micro_alpha_shannons_sample_type

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.Combined_order.shannon.pw <- pairwise.wilcox.test(micro.alpha_div_meta$Shannon, micro.alpha_div_meta$Combined_order, p.adjust.method = "BH")
micro.Combined_order.shannon.pw # p-values in the matrix

#Combine plots into panel for publication
# Extract the legend from one plot
#legend <- get_legend(micro_alpha_richness_sample_type)

# Combine plots into a single panel without the legend
Sample_Type_combined <- plot_grid(micro_alpha_richness_sample_type + theme(legend.position = "bottom"), 
                                  micro_alpha_shannons_sample_type + theme(legend.position = "bottom"), 
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
micro.alpha_div_fecal <- micro.alpha_div_meta[which(micro.alpha_div_meta$Sample_type=="Fecal"),]



# Change ordre of "Combined_order"
micro.alpha_div_fecal$Combined_order <- factor(micro.alpha_div_fecal$Combined_order, levels = c("Start_CONV","Midpoint_CONV","End_CONV", "Start_RWA",
                                                                                    "Midpoint_RWA","End_RWA","Final_washout"))
#Define Custom Colors do you can split RWA from Conv
#Do these all one color for each group instead of timepoints and group
custom_colors2 <- c(
  "Start_RWA" = "#0033ff",   
  "Midpoint_RWA" = "#0033ff",     
  "End_RWA" = "#0033ff",     
  "Start_CONV" = "#ff0000", 
  "Midpoint_CONV" = "#ff0000",  
  "End_CONV" = "#ff0000",
  "Final_washout" = "#999999")



fecal_micro_richness <- ggplot(micro.alpha_div_fecal, aes(x= Combined_order, y= Observed, fill = Combined_order, colour = Combined_order)) +
  theme_bw() + 
  labs(y= "Observed ASVs") +
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
fecal_micro_richness

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.Combined_order.richness.pw <- pairwise.wilcox.test(micro.alpha_div_fecal$Observed, micro.alpha_div_fecal$Combined_order, p.adjust.method = "BH")
micro.Combined_order.richness.pw # p =.97


## Shannon's #
fecal_micro_shannons <- ggplot(micro.alpha_div_fecal, aes(x= Combined_order, y= Shannon, fill = Combined_order, colour = Combined_order)) +
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
fecal_micro_shannons

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.Combined_order.shannon.pw <- pairwise.wilcox.test(micro.alpha_div_fecal$Shannon, micro.alpha_div_fecal$Combined_order, p.adjust.method = "BH")
micro.Combined_order.shannon.pw # p-values in the matrix



#Combine plots into panel for publication
# Extract the legend from one plot
legend_fecal <- get_legend(fecal_micro_richness)

# Combine plots into a single panel without the legend
fecal_micro_combined <- plot_grid(fecal_micro_richness + theme(legend.position = "bottom"), 
                                    fecal_micro_shannons + theme(legend.position = "bottom"), 
                                    ncol = 2)

# Combine the title, the combined plots, and the legend
fecal_micro_final <- plot_grid(
  ggdraw() + draw_label("Alpha Diversity of Fecal Samples by Time", size = 18), 
  fecal_micro_combined, 
  legend, 
  ncol = 1, 
  rel_heights = c(0.1, 1, 0.1))
fecal_micro_final






####Same thing without the Final Washout

NWmicro.alpha_div_fecal <- micro.alpha_div_meta[which(micro.alpha_div_fecal$Combined_order != "Final_washout"), ]



# Change ordre of "Combined_order"
NWmicro.alpha_div_fecal$Combined_order <- factor(NWmicro.alpha_div_fecal$Combined_order, levels = c("Start_CONV","Midpoint_CONV","End_CONV", "Start_RWA",
                                                                                                "Midpoint_RWA","End_RWA"))

#Do these all one color for each group instead of timepoints and group
NWcustom_colors <- c(
  "Start_RWA" = "#0033ff",   
  "Midpoint_RWA" = "#0033ff",     
  "End_RWA" = "#0033ff",     
  "Start_CONV" = "#ff0000", 
  "Midpoint_CONV" = "#ff0000",  
  "End_CONV" = "#ff0000")



NWfecal_micro_richness <- ggplot(NWmicro.alpha_div_fecal, aes(x= Combined_order, y= Observed, fill = Combined_order, colour = Combined_order)) +
  theme_bw() + 
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  scale_fill_manual(values = NWcustom_colors) +
  scale_colour_manual(values = NWcustom_colors) +
  theme(legend.position = "right",
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
NWfecal_micro_richness

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.Combined_order.richness.pw <- pairwise.wilcox.test(micro.alpha_div_fecal$Observed, micro.alpha_div_fecal$Combined_order, p.adjust.method = "BH")
micro.Combined_order.richness.pw # p-values in the matrix


## Shannon's #
NWfecal_micro_shannons <- ggplot(NWmicro.alpha_div_fecal, aes(x= Combined_order, y= Shannon, fill = Combined_order, colour = Combined_order)) +
  theme_bw() + 
  labs(y= "Shannon's diversity") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  scale_fill_manual(values = custom_colors2) +
  scale_colour_manual(values = custom_colors2) +
  theme(legend.position = "right",
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
NWfecal_micro_shannons

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
NWmicro.Combined_order.shannon.pw <- pairwise.wilcox.test(NWmicro.alpha_div_fecal$Shannon, NWmicro.alpha_div_fecal$Combined_order, p.adjust.method = "BH")
NWmicro.Combined_order.shannon.pw # p-values in the matrix

#Combine plots into panel for publication
# Extract the legend from one plot
NWlegend_fecal <- get_legend(NWfecal_micro_richness)

# Combine plots into a single panel without the legend
NWfecal_micro_combined <- plot_grid(NWfecal_micro_richness + theme(legend.position = "bottom"), 
                                  NWfecal_micro_shannons + theme(legend.position = "bottom"), 
                                  ncol = 2)

# Combine the title, the combined plots, and the legend
NWfecal_micro_final <- plot_grid(
  ggdraw() + draw_label("Alpha Diversity of Fecal Samples by Time", size = 18), 
  NWfecal_micro_combined, 
  legend, 
  ncol = 1, 
  rel_heights = c(0.1, 1, 0.1))
NWfecal_micro_final

###############################################
###############   ALPHA R2 STATS   ############
###############################################

# Helper function to extract R2 from a linear model
get_r2 <- function(model) summary(model)$r.squared

### ---- 1. MICROBIOME: SAMPLE TYPE COMPARISONS ---- ###
# Observed richness
lm_obs_st <- lm(Observed ~ Sample_type, data = micro.alpha_div_meta)
R2_obs_sampletype <- get_r2(lm_obs_st)

# Shannon
lm_shan_st <- lm(Shannon ~ Sample_type, data = micro.alpha_div_meta)
R2_shannon_sampletype <- get_r2(lm_shan_st)

R2_sampletype_results <- data.frame(
  Measure = c("Observed", "Shannon"),
  R2 = c(R2_obs_sampletype, R2_shannon_sampletype)
)
R2_sampletype_results


### ---- 2. FECAL ONLY: COMBINED ORDER ---- ###
fec <- micro.alpha_div_fecal

# Observed richness
lm_obs_fec <- lm(Observed ~ Combined_order, data = fec)
R2_obs_fecal <- get_r2(lm_obs_fec)

# Shannon
lm_shan_fec <- lm(Shannon ~ Combined_order, data = fec)
R2_shannon_fecal <- get_r2(lm_shan_fec)

R2_fecal_results <- data.frame(
  Measure = c("Observed", "Shannon"),
  R2 = c(R2_obs_fecal, R2_shannon_fecal)
)
R2_fecal_results


### ---- 3. FECAL (NO WASHOUT) ---- ###
fecNW <- NWmicro.alpha_div_fecal

# Observed richness
lm_obs_fecNW <- lm(Observed ~ Combined_order, data = fecNW)
R2_obs_fecalNW <- get_r2(lm_obs_fecNW)

# Shannon
lm_shan_fecNW <- lm(Shannon ~ Combined_order, data = fecNW)
R2_shannon_fecalNW <- get_r2(lm_shan_fecNW)

R2_fecalNW_results <- data.frame(
  Measure = c("Observed", "Shannon"),
  R2 = c(R2_obs_fecalNW, R2_shannon_fecalNW)
)
R2_fecalNW_results


### ---- PRINT ALL RESULTS CLEANLY ---- ###
cat("\n\n=========== R2 RESULTS FOR ALPHA DIVERSITY ===========\n")

cat("\n-- Sample Type (Microbiome) --\n")
print(R2_sampletype_results)

cat("\n-- Fecal Only (All Timepoints) --\n")
print(R2_fecal_results)

cat("\n-- Fecal Only (No Washout) --\n")
print(R2_fecalNW_results)



