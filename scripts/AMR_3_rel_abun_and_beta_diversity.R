#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################
#############################################################################################
#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################
#############################################################################################
#Make new phyloseq object
beta_data_noSNP <- data_noSNP

#NORMALIZE the beta_data_noSNP for beta diversity analysis
beta_data_noSNP <- prune_taxa(taxa_sums(beta_data_noSNP) > 0, data_noSNP)
any(taxa_sums(beta_data_noSNP)==0) # QUADRUPLE CHECKING - nope good.

beta_data_noSNP.css <- phyloseq_transform_css(beta_data_noSNP, log = F)
beta_data_noSNP.css.df <- as(sample_data(beta_data_noSNP.css), "data.frame")

rel_abund <- transform_sample_counts(beta_data_noSNP.css, function(x) {x/sum(x)}*100)

#Agglomerate
beta_ra_group <- tax_glom(rel_abund, taxrank = "group")
beta_ra_group # 502 groups
beta_ra_group_melt <- psmelt(beta_ra_group)

beta_ra_class <- tax_glom(rel_abund, taxrank = "class")
beta_ra_class # 40 classes, 216 samples
beta_ra_class_melt_temp <- psmelt(beta_ra_class)


######## Re-label low abundance taxa into single category
# convert Phylum to a character vector from a factor because R
beta_ra_class_melt_temp$class <- as.character(beta_ra_class_melt_temp$class)

# group dataframe by Class, calculate median rel. abundance
medians <- ddply(beta_ra_class_melt_temp, ~class, function(x) c(median=median(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
remainder <- medians[medians$median <= 1,]$class

# change the name of taxa to whose rel. abund. is less than 1% to "Low abundance phyla (<1%)"
beta_ra_class_melt_temp[beta_ra_class_melt_temp$class %in% remainder,]$class <- 'low abundance class (<1%)'

# Determine the number of taxa
length(unique(beta_ra_class_melt_temp$class))
#8



library(cowplot)

##
### Clustering at group level
##

# Use beta_ra_group to cluster at the group level
ps_AMR.dist <- vegdist(t(otu_table(beta_ra_group)), method = "bray")  
ps_AMR.hclust <- hclust(ps_AMR.dist)
plot(ps_AMR.hclust) # example plot

# Extract data as dendrogram
ps_AMR.dendro <- as.dendrogram(ps_AMR.hclust)
ps_AMR.dendro.data <- dendro_data(ps_AMR.dendro, type = "rectangle")

# Sample names in order based on clustering
sample_names_ps_AMR <- ps_AMR.dendro.data$labels$label

# Standardize metadata sample IDs
AMR_mapfile <- as(sample_data(data_noSNP), "data.frame")
rownames(AMR_mapfile) <- gsub("-", ".", rownames(AMR_mapfile))
AMR_mapfile$Sample_ID <- rownames(AMR_mapfile)

# Add metadata
ps_AMR.dendro.data$labels <- left_join(ps_AMR.dendro.data$labels,
                                       AMR_mapfile,
                                       by = c("label" = "Sample_ID"))

# Setup the data, so that the layout is inverted
segment_data_ps_AMR <- with(
  segment(ps_AMR.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend)
)

# Use the dendrogram label data to position the gene labels
gene_pos_table_ps_AMR <- with(ps_AMR.dendro.data$labels,
                              data.frame(y_center = x, gene = as.character(label), x = y,
                                         trt = as.character(Trial), order = as.character(Combined_order),
                                         height = 1)
)

# Table to position the samples
sample_pos_table_ps_AMR <- data.frame(sample = sample_names_ps_AMR) %>%
  dplyr::mutate(x_center = (1:dplyr::n()), width = 1)

##
######## Relative abundance bar plot #########
##

ra_class_melt <- beta_ra_class_melt_temp

# Use class melted data and add gene locations
joined_ra_class_ps_class_melt <- ra_class_melt %>%
  left_join(gene_pos_table_ps_AMR, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_ps_AMR, by = c("Sample" = "sample")) 

# Calculate the mean relative abundance of each class taxa, sort by most abundant to least
ra_class_melt$class <- as.character(ra_class_melt$class)

factor_by_abund <- ra_class_melt %>%
  dplyr::group_by(class) %>%
  dplyr::summarize(median_class = median(Abundance)) %>%
  arrange(-median_class)

# Use the sorted class levels to clean up the order of taxa in the relative abundance plots
joined_ra_class_ps_class_melt$class <- factor(joined_ra_class_ps_class_melt$class,
                                              levels = as.character(factor_by_abund$class))

# Limits for the vertical axes
gene_axis_limits_ps_class <- with(
  gene_pos_table_ps_AMR, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1)

# Color palette for color blind friendly
tol21 <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
           "#117777", "#DDDD77","#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
           "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
           "#771122", "#AA4455", "#DD7788")

# Compute totals per bar for an outer outline
bar_totals <- joined_ra_class_ps_class_melt %>%
  dplyr::group_by(x_center) %>%
  dplyr::summarize(Total = sum(Abundance), .groups = "drop")










# Relative abundance plot
plt_rel_class_ps_class <- ggplot(joined_ra_class_ps_class_melt,
                                 aes(x = x_center, y = Abundance, fill = class)) + 
  coord_flip() +
  geom_bar(stat = "identity", alpha = 0.95, width = 1, colour = NA) +
  geom_col(data = bar_totals, aes(x = x_center, y = Total),
           fill = NA, colour = "black", width = 1, linewidth = 0.5) +
  scale_fill_manual(values = tol21) +
  scale_x_continuous(breaks = sample_pos_table_ps_AMR$x_center, 
                     labels = sample_pos_table_ps_AMR$sample, 
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0,0)) +   # remove extra white space left/right
  labs(x = "", y = "Relative Abundance", fill = "Antimicrobial Class") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.box.margin = margin(0, 10, 0, 10),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(-0.9, 0.5, 0.3, 0), "cm"))

##
######## Dendrogram #########
##

# Force order of legend labels exactly how you want them
legend_colors <- c("Start_No_Claim", "Midpoint_No_Claim", "End_No_Claim",
                   "Start_RWA", "Midpoint_RWA", "End_RWA",
                   "No Claim Beef", "RWA Beef")

c_colors <- c("Start_No_Claim"="#660000", "Midpoint_No_Claim"="#CC0066",
              "End_No_Claim" = "#ff9999", 
              "Start_RWA" = "#000999", "Midpoint_RWA" = "#006699", "End_RWA" = "#66CCFF",
              "No Claim Beef"= "#DDDD77", "RWA Beef" = "#117744")

plt_dendr_class_ps_class <- ggplot(segment_data_ps_AMR) + 
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),
               lineend = "round", linejoin = "round") +
  geom_point(data = gene_pos_table_ps_AMR,
             aes(x, y_center, color = order),
             size = 6, shape = 15, stroke = 0,
             position = position_nudge(x = 0.04)) +
  labs(x = "Ward's Distance", y = "", colour = NULL, title = "") +
  scale_y_discrete(expand = c(0,0,0,0)) +
  scale_x_reverse() + 
  scale_color_manual(values= c_colors, breaks = legend_colors) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +  # wrap into 3 rows in order
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.box.margin = margin(0, 10, 0, 10),
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 14),
        axis.ticks.length = unit(0, "cm"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
        axis.text.x = element_text(size = 12, colour = "black"))

##
######## Combine dendrogram + RA #########
##

combined_main <- plot_grid(plt_dendr_class_ps_class,
                           plt_rel_class_ps_class,
                           axis = "tb",
                           align = 'h',
                           nrow = 1,
                           rel_widths = c(0.5, 1.5))

# Final object ready to embed in a multi-panel figure
class_ps_class_plots <- combined_main
class_ps_class_plots

# Save final version (adjust dimensions as needed)
ggsave("AMR_RA_Final.png", class_ps_class_plots,
       width = 18, height = 10, units = "in", dpi = 300)
