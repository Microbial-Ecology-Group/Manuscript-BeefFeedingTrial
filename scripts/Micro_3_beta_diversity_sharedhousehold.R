# AMR beta diversity by shared household status
#############################################################################################
####################################   CSS TRANSFORM    ###################################
data_micro.ps <- prune_taxa(taxa_sums(data_micro.ps) > 0, data_micro.ps)
any(taxa_sums(data_micro.ps)==0) # QUADRUPLE CHECKING - nope good.

data_micro.ps.css <- phyloseq_transform_css(data_micro.ps, log = F)
data_micro.ps.css.df <- as(sample_data(data_micro.ps.css), "data.frame")




###
####
# Only subset "Habitating" ######
####
### 

Habitating_data_micro.ps.css <- subset_samples(data_micro.ps.css, Habitating!="")

Habitating_data_micro.ps.css <- prune_taxa(taxa_sums(Habitating_data_micro.ps.css) > 0, Habitating_data_micro.ps.css)


Habitating_data_micro.ps.css.df <- as(sample_data(Habitating_data_micro.ps.css),"data.frame")
Habitating_data_micro.ps.css.dist <- vegdist(t(otu_table(Habitating_data_micro.ps.css)), method = "bray")
Habitating_data_micro.ps.css.ord <- vegan::metaMDS(comm = t(otu_table(Habitating_data_micro.ps.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plt_Habitating <- plot_ordination(Habitating_data_micro.ps.css, Habitating_data_micro.ps.css.ord, color = "Habitating", shape = "Habitating") +
  theme_bw() +
  labs(title ="Habitating") +
  stat_ellipse(aes(fill= Habitating), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =22, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 24),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_Habitating 

# permanova test of beta diversity, by Habitating
Habitating_data_micro.ps.css.adonis <- adonis2(Habitating_data_micro.ps.css.dist ~Habitating, data = Habitating_data_micro.ps.css.df, permutations = 9999)
Habitating_data_micro.ps.css.adonis #


## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Habitating.permanova <- pairwise.adonis(Habitating_data_micro.ps.css.dist, Habitating_data_micro.ps.css.df$Habitating, perm = 9999, p.adjust.m = "BH")
Habitating.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Habitating.disper <- betadisper(Habitating_data_micro.ps.css.dist, Habitating_data_micro.ps.css.df$Habitating)
Habitating.permdisp <- permutest(Habitating.disper, permutations = 9999, pairwise = T)
Habitating.permdisp # looks like a few are significant



###
####
# Only subset "Co.Habitating" #####
####
### 

Co.Habitating_data_micro.ps.css <- subset_samples(data_micro.ps.css, Co.Habitating!="")

Co.Habitating_data_micro.ps.css <- prune_taxa(taxa_sums(Co.Habitating_data_micro.ps.css) > 0, Co.Habitating_data_micro.ps.css)


Co.Habitating_data_micro.ps.css.df <- as(sample_data(Co.Habitating_data_micro.ps.css),"data.frame")
Co.Habitating_data_micro.ps.css.dist <- vegdist(t(otu_table(Co.Habitating_data_micro.ps.css)), method = "bray")
Co.Habitating_data_micro.ps.css.ord <- vegan::metaMDS(comm = t(otu_table(Co.Habitating_data_micro.ps.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plt_Co.Habitating <- plot_ordination(Co.Habitating_data_micro.ps.css, Co.Habitating_data_micro.ps.css.ord, color = "Co.Habitating", shape = "Co.Habitating") +
  theme_bw() +
  labs(title ="Co.Habitating") +
  stat_ellipse(aes(fill= Co.Habitating), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =22, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 24),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_Co.Habitating 

# permanova test of beta diversity, by Co.Habitating
Co.Habitating_data_micro.ps.css.adonis <- adonis2(Co.Habitating_data_micro.ps.css.dist ~Co.Habitating, data = Co.Habitating_data_micro.ps.css.df, permutations = 9999)
Co.Habitating_data_micro.ps.css.adonis #


## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Co.Habitating.permanova <- pairwise.adonis(Co.Habitating_data_micro.ps.css.dist, Co.Habitating_data_micro.ps.css.df$Co.Habitating, perm = 9999, p.adjust.m = "BH")
Co.Habitating.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Co.Habitating.disper <- betadisper(Co.Habitating_data_micro.ps.css.dist, Co.Habitating_data_micro.ps.css.df$Co.Habitating)
Co.Habitating.permdisp <- permutest(Co.Habitating.disper, permutations = 9999, pairwise = T)
Co.Habitating.permdisp # looks like a few are significant


#
##
# Calculating average pairwise distances for each group ####
##
#

## Indices of samples from shared households
shared_indices <- which(Habitating_data_micro.ps.css.df$Habitating == "Shared")
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])


## Indices of samples from Shared_1
shared_1_indices <- which(Habitating_data_micro.ps.css.df$Co.Habitating == "Shared_1")
# Extract pairwise distances for shared households
shared_1_indices <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_1_indices, shared_1_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_1_indices[lower.tri(shared_1_indices)])


## Indices of samples from Shared_2
shared_2_indices <- which(Habitating_data_micro.ps.css.df$Co.Habitating == "Shared_2")
# Extract pairwise distances for shared households
shared_2_indices <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_2_indices, shared_2_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_2_indices[lower.tri(shared_2_indices)])

## Nonshared ###
# Indices of samples from nonshared households
nonshared_indices <- which(Habitating_data_micro.ps.css.df$Habitating == "Control")
# Extract pairwise distances for nonshared households
nonshared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[nonshared_indices, nonshared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(nonshared_distances[lower.tri(nonshared_distances)])




