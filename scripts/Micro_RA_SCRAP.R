#Scrap trial. Still does not work


library(metagMisc)
library(cowplot)

#############################################################################################
####################################   CSS TRANSFORM    ###################################
#############################################################################################

# SUBSET: remove Final_washout samples only
micro_subset <- subset_samples(data_micro.ps, Combined_order != "Final_Washout")
micro.ps <- micro_subset


#############################################################################################
##############################   RELATIVE ABUNDANCE + DENDROGRAM   #########################
#############################################################################################

# CSS Transform
micro.ps.css <- phyloseq_transform_css(micro.ps, log = FALSE)
micro.ps.css.df <- as(sample_data(micro.ps.css), "data.frame")


################################################################################################
# Remove unclassified taxa
to_remove <- apply(tax_table(micro.ps.css), 1, function(x)
  any(x %in% c("unclassified Unassigned", "unclassified Bacteria")))
micro.ps.css.filtered <- prune_taxa(!to_remove, micro.ps.css)


#############################################################################################
##############################   CLEAN AGGLOMERATION & MELT (FIXED)  ########################
#############################################################################################

# --- Step 1: Double-check taxonomy integrity ---
taxa_tab <- tax_table(micro.ps.css.filtered) %>% as.data.frame()
if (!"Order" %in% colnames(taxa_tab)) stop("❌ No 'Order' column found in taxonomy table.")
taxa_tab$Order[is.na(taxa_tab$Order) | taxa_tab$Order == ""] <- "Unclassified"
tax_table(micro.ps.css.filtered) <- as.matrix(taxa_tab)

# --- Step 2: Recalculate relative abundance ---
rel_abund <- transform_sample_counts(micro.ps.css.filtered, function(x) (x / sum(x)) * 100)

# --- Step 3: Agglomerate at Order level ---
ra_Order <- tax_glom(rel_abund, taxrank = "Order")

# --- Step 4: Melt to long format ---
ra_Order_melt_temp <- psmelt(ra_Order)

# --- Step 5: Verify ---
message("✅ Unique Orders detected right after psmelt: ",
        length(unique(ra_Order_melt_temp$Order)))
print(head(unique(ra_Order_melt_temp$Order)))


#############################################################################################
###########################   LOW-ABUNDANCE FIX + SAFE NORMALIZATION   ######################
#############################################################################################

# Fix psmelt renaming issue
if ("sample_Sample" %in% names(ra_Order_melt_temp)) {
  ra_Order_melt_temp$Sample <- ra_Order_melt_temp$sample_Sample
  ra_Order_melt_temp <- ra_Order_melt_temp[, !(names(ra_Order_melt_temp) == "sample_Sample")]
}

# Standardize column names
sample_col <- names(ra_Order_melt_temp)[grepl("Sample", names(ra_Order_melt_temp), ignore.case = TRUE)][1]
order_col  <- names(ra_Order_melt_temp)[grepl("Order",  names(ra_Order_melt_temp), ignore.case = FALSE)][1]
names(ra_Order_melt_temp)[names(ra_Order_melt_temp) == sample_col] <- "Sample"
names(ra_Order_melt_temp)[names(ra_Order_melt_temp) == order_col]  <- "Order"

# Identify low-abundance Orders (<1%)
ra_Order_melt_temp$Order <- as.character(ra_Order_melt_temp$Order)
medians <- ddply(ra_Order_melt_temp, ~Order, function(x) c(median = median(x$Abundance)))
remainder <- medians[medians$median <= 0.01,]$Order

# Remove duplicates created by psmelt
dupes <- names(ra_Order_melt_temp)[duplicated(names(ra_Order_melt_temp))]
if (length(dupes) > 0) {
  message("Removing duplicate columns: ", paste(dupes, collapse = ", "))
  ra_Order_melt_temp <- ra_Order_melt_temp[, !duplicated(names(ra_Order_melt_temp))]
}

# --- Collapse and normalize ---
ra_Order_melt_temp <- ra_Order_melt_temp %>%
  dplyr::mutate(Order = ifelse(Order %in% remainder, "Low abundance Orders (<1%)", Order)) %>%
  dplyr::group_by(Sample, Order) %>%
  dplyr::summarise(Abundance = sum(as.numeric(Abundance)), .groups = "drop") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(Abundance = Abundance / sum(Abundance) * 100) %>%
  dplyr::ungroup()

# ✅ Verify totals and number of Orders
ra_Order_melt_temp %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(total = sum(Abundance)) %>%
  dplyr::summarise(range_total = range(total)) %>%
  print()
message("✅ Orders after collapse: ", length(unique(ra_Order_melt_temp$Order)))


#############################################################################################
###################################   CLUSTERING   #########################################
#############################################################################################

ps_micro.dist <- vegdist(t(otu_table(micro.ps.css.filtered)), method = "bray")
ps_micro.hclust <- hclust(ps_micro.dist)
ps_micro.dendro <- as.dendrogram(ps_micro.hclust)
ps_micro.dendro.data <- dendro_data(ps_micro.dendro, type = "rectangle")

dendro_order <- ps_micro.dendro.data$labels$label
mapfile_micro <- as(sample_data(micro.ps.css.filtered), "data.frame")
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
                                      Combined_order = as.character(Combined_order),
                                      height = 1))
sample_pos_table_ps_micro <- data.frame(sample = dendro_order) %>%
  mutate(x_center = 1:nrow(.), width = 1)


#############################################################################################
##############################   RELATIVE ABUNDANCE PLOT   #################################
#############################################################################################

ra_Order_melt_temp <- ra_Order_melt_temp %>%
  dplyr::mutate(Sample = as.character(Sample)) %>%
  dplyr::filter(Sample %in% dendro_order)

joined_ra_Order_melt <- ra_Order_melt_temp %>%
  dplyr::mutate(Sample = factor(Sample, levels = dendro_order)) %>%
  dplyr::left_join(pos_table_ps_micro, by = c("Sample" = "sample")) %>%
  dplyr::left_join(sample_pos_table_ps_micro, by = c("Sample" = "sample")) %>%
  dplyr::arrange(Sample)

joined_ra_Order_melt$Order <- as.character(joined_ra_Order_melt$Order)
factor_by_abund <- joined_ra_Order_melt %>%
  dplyr::group_by(Order) %>%
  dplyr::summarize(median_Order = median(Abundance)) %>%
  dplyr::arrange(-median_Order)
joined_ra_Order_melt$Order <- factor(joined_ra_Order_melt$Order,
                                     levels = as.character(factor_by_abund$Order))

bar_totals_micro <- joined_ra_Order_melt %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarize(Total = sum(Abundance), .groups = "drop") %>%
  dplyr::mutate(Sample = factor(Sample, levels = dendro_order))

joined_ra_Order_melt <- joined_ra_Order_melt %>%
  dplyr::rename(sample = Sample)
bar_totals_micro <- bar_totals_micro %>%
  dplyr::rename(sample = Sample)

joined_ra_Order_melt$sample <- factor(joined_ra_Order_melt$sample, levels = dendro_order)
bar_totals_micro$sample <- factor(bar_totals_micro$sample, levels = dendro_order)

joined_ra_Order_melt <- joined_ra_Order_melt %>%
  dplyr::filter(sample %in% dendro_order) %>%
  dplyr::mutate(sample = factor(sample, levels = dendro_order))
bar_totals_micro <- bar_totals_micro %>%
  dplyr::filter(sample %in% dendro_order) %>%
  dplyr::mutate(sample = factor(sample, levels = dendro_order))

# ✅ Color palette (21 colors)
MMColors <- c("#77CCCC","#E73f74","#11A579","#ffcccc","#F2B701","#117777",
              "#00bfff","#3969AC","#00ffff","#764E9F","#ED645A","#000fff",
              "#32cd32","#ff4500","#ff1493","#1E90FF","#771155","#ADFF2F",
              "#8A2BE2","#00BFFF","#7B68EE")

# Plot
plt_rel_ra_Order <- ggplot(joined_ra_Order_melt,
                           aes(x = factor(sample, levels = dendro_order),
                               y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", width = 1, colour = NA) +
  geom_col(data = bar_totals_micro,
           aes(x = factor(sample, levels = dendro_order), y = Total),
           fill = NA, colour = "black", width = 1, linewidth = 0.5,
           inherit.aes = FALSE) +
  coord_flip() +
  scale_fill_manual(values = MMColors, na.value = "grey80") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Relative Abundance", fill = "Order") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.75),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, colour = "black"),
    plot.margin = unit(c(-0.9, 0.5, 0.3, 0), "cm")
  )

plt_rel_ra_Order




#############################################################################################
##############################   DENDROGRAM PLOT   #########################################
#############################################################################################

# Legend order and colors (8 categories, same as AMR)
legend_colors_micro <- c("Start No Claim Fecal", "Midpoint No Claim Fecal", "End No Claim Fecal",
                         "Start RWA Fecal", "Midpoint RWA Fecal", "End RWA Fecal",
                         "No Claim Beef", "RWA Beef")

c_colors_micro <- c("Start No Claim Fecal" = "#660000",
                    "Midpoint No Claim Fecal" = "#CC0066",
                    "End No Claim Fecal" = "#ff9999",
                    "Start RWA Fecal" = "#000999",
                    "Midpoint RWA Fecal" = "#006699",
                    "End RWA Fecal" = "#66CCFF",
                    "No Claim Beef" = "#DDDD77",
                    "RWA Beef" = "#117744")

plt_dendr_Order <- ggplot(segment_data_ps_micro) +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),
               lineend = "round", linejoin = "round") +
  geom_point(data = pos_table_ps_micro,
             aes(x, y_center, color = Combined_order),
             size = 6, shape = 15, stroke = 0,
             position = position_nudge(x = 0.03)) +
  labs(x = "Ward's Distance", y = "", colour = NULL, title = "") +
  scale_y_discrete(expand = c(0,0,0,0)) +
  scale_x_reverse() +
  scale_color_manual(values = c_colors_micro, breaks = legend_colors_micro) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
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
plt_dendr_Order



#############################################################################################
##############################   COMBINE + SAVE (FIXED)   ##################################
#############################################################################################

# Ensure both plots share identical vertical coordinate space before combining
aligned_plots <- cowplot::align_plots(
  plt_dendr_Order,
  plt_rel_ra_Order,
  align = "hv",      # align both horizontally and vertically
  axis = "tblr"      # align on all axes, not just top/bottom
)

# Combine with fixed relative widths
combined_micro <- cowplot::plot_grid(
  aligned_plots[[1]],
  aligned_plots[[2]],
  nrow = 1,
  rel_widths = c(0.5, 1.5),
  align = "h",
  axis = "tb"
)

# Ensure margins are even so no drift occurs
combined_micro <- combined_micro + theme(plot.margin = margin(10, 15, 10, 15))

# View the final figure
combined_micro

