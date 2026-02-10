####Ordination####

#OCT 7 2024, per PM, Ordination done at group level!

cb_colors <- c("#0072B2","#882255","#332288", "#117733","#D55E00","#E69F00", "#56B4E9", "#009E73", "#F0E442","#CC79A7")
trialcolors <- c("#e69f00", "#117733", "#FF0000", "#0000CC")
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

beta_ra_group_css <- beta_ra_group_melt
beta_ra_group_css_df <- data.frame(beta_ra_group_css)


###
####
##### All samples- Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
#data.dist <- vegdist(decostand(t(otu_table(beta_ra_group_css)), "hell"), "euclidean") 
data.dist <- vegdist(t(otu_table(beta_ra_group_css)), method = "bray")
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "bray", try = 10, trymax = 999, autotransform = F)
#plot_ordination(data.css, data.ord, color = "Sample_Type") 
plt_ord_by_sampletype <- plot_ordination(beta_ra_group_css, data.ord, color = "Sample_Type") +
  theme_bw() +
  scale_fill_manual(values = cb_colors)+
  scale_color_manual(values = cb_colors)+
  labs(title = "ARGs by Sample Type") +
  stat_ellipse(aes(color = Sample_Type), linetype = 1, size = 0.75) +
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"),
    strip.background = element_rect(fill = "none", size = 1.0),
    strip.text = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(colour = "black", size = 0.7),
    plot.title = element_text(size = 18),
    panel.border = element_rect(colour = "black", size = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
plt_ord_by_sampletype


# permanova test of beta diversity, by End_trial
data.adonis <- adonis2(data.dist ~ Sample_Type, data = beta_ra_group_css_df, permutations = 9999)
data.adonis # 1e-04, R2 = 0.47731 

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
sample_type.disper <- betadisper(data.dist, beta_ra_group_css_df$Sample_Type)
sample_type.permdisp <- permutest(sample_type.disper, permutations = 9999, pairwise = T)
sample_type.permdisp # 


###
####
##### Trial comparison + meat #######
####
### 

# Make new phyloseq object so that you can subset out the NA values in Trial column
noNA_beta_ra_group.css <- beta_ra_group_css

# Subset out the NA values
noNA_beta_ra_group.css <- subset_samples(noNA_beta_ra_group.css, !is.na(Trial))

# Below step was necessary to get ordination to work
# Remove taxa with zero counts across all samples
noNA_beta_ra_group.css <- prune_taxa(taxa_sums(noNA_beta_ra_group.css) > 0, noNA_beta_ra_group.css)

#MM added the below two lines 6/21
#subset_vars <- c("RWA Fecal", "Conventional Fecal", "RWA Beef", "Conventional Beef")
noNA_beta_ra_group.css.df <- as(sample_data(noNA_beta_ra_group.css),"data.frame")
#subset_df <- noNA_beta_ra_group.css.df[noNA_beta_ra_group.css.df$Trial %in% subset_vars & !is.na(noNA_beta_ra_group.css.df$Trial), "Trial"]


sample_data(noNA_beta_ra_group.css)$Trial <- factor(sample_data(noNA_beta_ra_group.css)$Trial, levels = c("Conventional Beef","RWA Beef","RWA Fecal", "Conventional Fecal"))

noNA.dist <- vegdist(t(otu_table(noNA_beta_ra_group.css)), method = "bray")
noNA.ord <- vegan::metaMDS(comm = t(otu_table(noNA_beta_ra_group.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plt_ord_ps_class <- plot_ordination(noNA_beta_ra_group.css, noNA.ord, color = "Trial") +
  theme_bw() +
  scale_fill_manual(values = trialcolors)+
  scale_color_manual(values = trialcolors)+
  labs(title ="AMR by Trial") +
  stat_ellipse(aes(color= Trial), lty = 1, size = 0.75) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =22, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_ps_class 

# permanova test of beta diversity, by End_trial
noNA.adonis <- adonis2(noNA.dist ~Trial, data = noNA_beta_ra_group.css.df, permutations = 9999)
noNA.adonis #

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
trial.disper <- betadisper(noNA.dist, noNA_beta_ra_group.css.df$Trial)
trial.permdisp <- permutest(trial.disper, permutations = 9999, pairwise = T)
trial.permdisp # 






####
####### Combined order NMDS ####
####
###


ps_group.css <- tax_glom(beta_data_noSNP.css, taxrank = "group")
ps_group.css <- subset_samples(ps_group.css , Combined_order != "")
ps_group.css <- subset_samples(ps_group.css, Sample_Type != "Beef sample")

ps_group.css <- prune_taxa(taxa_sums(ps_group.css) > 0, ps_group.css)
any(taxa_sums(ps_group.css)==0) # QUADRUPLE CHECKING - nope good.


ps_group.css.df <- as(sample_data(ps_group.css), "data.frame")

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
ps_group.dist <- vegdist(t(otu_table(ps_group.css)), method = "bray")
ps_group.ord <- vegan::metaMDS(comm = t(ps_group.dist), distance = "bray", try = 10, trymax = 999, autotransform = F)


#Reorder the levels in the legend
sample_data(ps_group.css)$Combined_order <- factor(sample_data(ps_group.css)$Combined_order,
  levels = c("Start_CONV", "Midpoint_CONV", "End_CONV", "Start_RWA", "Midpoint_RWA", "End_RWA", "Final_washout"))

#Custom Colors
comb_colors <- c(
  "Start_RWA" = "#99ccff",   
  "Midpoint_RWA" = "#0033ff",     
  "End_RWA" = "#000099",     
  "Start_CONV" = "#ff9999", 
  "Midpoint_CONV" = "#ff0000",  
  "End_CONV" = "#cc0033",
  "Final_washout" = "#999999")



plt_ord_by_sampletype <- plot_ordination(ps_group.css, ps_group.ord, color = "Combined_order") +
  theme_bw() +
  labs(title ="AMR Fecal by Time") +
  stat_ellipse(aes(color = Combined_order), lty = 1, size = 0.7) +
  scale_fill_manual(values = comb_colors)+
  scale_color_manual(values = comb_colors)+
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_by_sampletype 

# permanova test of beta diversity, by Combined_order
ps_group.adonis <- adonis2(ps_group.dist ~ Combined_order, data = ps_group.css.df, permutations = 9999)
ps_group.adonis #


## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
combined_order.permanova <- pairwise.adonis(ps_group.dist, ps_group.css.df$Combined_order, perm = 9999, p.adjust.m = "BH")
combined_order.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
combined_order.disper <- betadisper(ps_group.dist, ps_group.css.df$Combined_order)
combined_order.permdisp <- permutest(combined_order.disper, permutations = 9999, pairwise = T)
combined_order.permdisp # looks like a few are significant



