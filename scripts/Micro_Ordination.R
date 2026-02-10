# Beta diversity
#############################################################################################
##############################         BETA DIVERSITY         #########################@#####
#############################################################################################
#############################################################################################
library(ggsci) # if not using this, make sure to switch out the palatte in the figures
library(metagMisc)
#############################################################################################
####################################   CSS TRANSFORM    ###################################
cb_colors <- c("#0072B2","#882255","#332288", "#117733","#D55E00","#E69F00", "#56B4E9", "#009E73", "#F0E442","#CC79A7")
trialcolors <- c("#e69f00", "#117733", "#FF0000", "#0000CC")


data_micro.ps <- prune_taxa(taxa_sums(data_micro.ps) > 0, data_micro.ps)
any(taxa_sums(data_micro.ps)==0) # QUADRUPLE CHECKING - nope good.

data_micro.ps.css <- phyloseq_transform_css(data_micro.ps, log = F)
data_micro.ps.css.df <- as(sample_data(data_micro.ps.css), "data.frame")

data_micro.ps.css <- transform_sample_counts(data_micro.ps.css, function(x) {x/sum(x)}*100)

# ordinate it based on Bray-Curtis
# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data.dist <- vegdist(decostand(t(otu_table(data_micro.ps.css)), "hell"), "euclidean")
data_micro.ps.css.ord <- vegan::metaMDS(comm = t(otu_table(data_micro.ps.css)), try = 50, trymax = 500, distance = "bray", autotransform = F)


# plot the ordination with 95% confidence ellipses coloured by groups
plt_ord_by_sampletype <- plot_ordination(data_micro.ps.css, data_micro.ps.css.ord, color = "sample_type") +
  theme_bw() +
  labs(title = "ASVs by Sample Type") +
  stat_ellipse(aes(color = sample_type), linetype = 1, size = 0.75) +
  scale_fill_manual(values = cb_colors)+
  scale_color_manual(values = cb_colors)+
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"),
    strip.background = element_rect(fill = "none", size = 1.0),
    strip.text = element_text(size = 24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(colour = "black", size = 0.7),
    plot.title = element_text(size = 18),
    panel.border = element_rect(colour = "black", size = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
plt_ord_by_sampletype

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
data_micro.ps.css.dist <- vegdist(t(otu_table(data_micro.ps.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
all.groups.permanova <- pairwise.adonis(data_micro.ps.css.dist, data_micro.ps.css.df$sample_type, perm = 9999, p.adjust.m = "BH")
all.groups.permanova 

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
all.groups.disper <- betadisper(data_micro.ps.css.dist, data_micro.ps.css$sample_type)
all.groups.permdisp <- permutest(all.groups.disper, permutations = 9999, pairwise = T)
all.groups.permdisp # looks like a few are significant




###
##### Combined order - Fecal samples #####
###
##
#

Fecal.data_micro.ps.css <- subset_samples(data_micro.ps.css, sample_type =="Fecal")
Fecal.data_micro.ps.css <-  prune_taxa(taxa_sums(Fecal.data_micro.ps.css) > 0, Fecal.data_micro.ps.css)

#### create d.f. for beta-diversity metadata
Fecal.data_micro.ps.css.df <- as(sample_data(Fecal.data_micro.ps.css),"data.frame")


# Change ordre of "Combined_order"
sample_data(Fecal.data_micro.ps.css)$Combined_order <- factor(
  sample_data(Fecal.data_micro.ps.css)$Combined_order,
  levels = c("Start_CONV", "Midpoint_CONV", "End_CONV", 
             "Start_RWA", "Midpoint_RWA", "End_RWA", 
             "Final_washout")
)


# ordinate it based on Bray-Curtis
Fecal.data_micro.ps.css.ord <- vegan::metaMDS(comm = t(otu_table(Fecal.data_micro.ps.css)), try = 50, trymax = 500, distance = "bray", autotransform = F)


# plot the ordination with 95% confidence ellipses coloured by groups
plot_ordination(Fecal.data_micro.ps.css, Fecal.data_micro.ps.css.ord, type = "samples", color = "Combined_order") +
  theme_bw() +
  labs(title = "ASVs by Time")+
  scale_color_manual(values = cb_colors)+
  scale_fill_manual(values = cb_colors)+
  stat_ellipse(aes(color = Combined_order),  # Remove fill aesthetic
               level = 0.95, 
               geom = "path",  # Use "path" for empty ellipses
               lty = 1, size = 0.75) +
    theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"),
    panel.border = element_rect(colour = "black", size = 0.75),
    legend.key.size = unit(2, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 18, vjust = 1.75),
    axis.title.x = element_text(size = 18, vjust = -1.5)
  )

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
Fecal.data_micro.ps.css.dist <- vegdist(t(otu_table(Fecal.data_micro.ps.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Fecal.groups.permanova <- adonis2(Fecal.data_micro.ps.css.dist ~ Combined_order, data = Fecal.data_micro.ps.css.df, permutations = 9999)
Fecal.groups.permanova

# Run pairwise PERMANOVA
Fecal.groups.pairwise <- pairwise.adonis2(Fecal.data_micro.ps.css.dist ~ Combined_order, data = Fecal.data_micro.ps.css.df, permutations = 9999, p.adjust.m = "BH") 
Fecal.groups.pairwise



#Fecal.groups.permanova <- pairwise.adonis(Fecal.data_micro.ps.css.dist, data_micro.ps.css.df$Combined_order, perm = 9999, p.adjust.m = "BH")
#Fecal.groups.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Fecal.groups.disper <- betadisper(Fecal.data_micro.ps.css.css.dist, Fecal.data_micro.ps.css.css.df$Combined_order)
Fecal.groups.permdisp <- permutest(Fecal.groups.disper, permutations = 9999, pairwise = T)
Fecal.groups.permdisp # looks like a few are significant

