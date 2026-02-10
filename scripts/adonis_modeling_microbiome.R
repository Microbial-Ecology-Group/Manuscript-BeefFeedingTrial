
Fecal_data_micro.ps <- subset_samples(data_micro.ps, sample_type=="Fecal" & Combined_order != "Final_Washout")

Fecal_data_micro.ps <- prune_taxa(taxa_sums(Fecal_data_micro.ps) > 0, Fecal_data_micro.ps)


Fecal_data_micro.ps.df <- as(sample_data(Fecal_data_micro.ps),"data.frame")
Fecal_data_micro.ps.dist <- vegdist(t(otu_table(Fecal_data_micro.ps)), method = "bray")
Fecal_data_micro.ps.ord <- vegan::metaMDS(comm = t(otu_table(Fecal_data_micro.ps)), distance = "bray", try = 10, trymax = 100, autotransform = F)
# plt_Combined_order <- plot_ordination(Fecal_data_micro.ps, Fecal_data_micro.ps.ord, color = "Combined_order", shape = "Combined_order") +
#   theme_bw() +
#   labs(title ="Combined_order") +
#   stat_ellipse(aes(fill= Combined_order), geom="polygon", alpha = 0.25) +
#   theme(legend.position = "right",
#         plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
#         strip.background = element_rect(fill= "grey91", size = 1.0),
#         strip.text = element_text(size =22, colour = "black"),
#         axis.text = element_text(size = 14, colour = "black"),
#         axis.title = element_text(size = 24),
#         axis.ticks = element_line(colour = "black", linewidth = 0.7),
#         plot.title = element_text(size = 30),
#         panel.border = element_rect(colour = "black", linewidth = 1.0),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank())
# plt_Combined_order 

# permanova test of beta diversity, by Combined_order
Fecal_data_micro.ps.adonis <- adonis2(Fecal_data_micro.ps.dist ~Combined_order, data = Fecal_data_micro.ps.df, permutations = 9999)
Fecal_data_micro.ps.adonis #


#adonis2(Fecal_data_micro.ps.dist ~ (trt * timepoint): treatment_order, strata=Fecal_data_micro.ps.df$participant_id , data = Fecal_data_micro.ps.df, permutations = 9999)

#We tested whether or not there was a difference in rsponse depending on the order in which the subjects received the treatments (3 way interaction - not significant.)

# Full model, global test
adonis2(
  Fecal_data_micro.ps.dist ~ trt * timepoint * treatment_order,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = NULL
)
# Full model, by margin
adonis2(
  Fecal_data_micro.ps.dist ~ trt * timepoint * treatment_order,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = "margin"
)

# no significance in three way interaction
# now specify each term individually without the three way interaction
# Global test
adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint + treatment_order + trt:timepoint + trt:treatment_order + timepoint:treatment_order,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = NULL
)
# By margin
adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint + treatment_order + trt:timepoint + trt:treatment_order + timepoint:treatment_order,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = "margin"
)

# No significant terms in model, reduce further to only include fixed effects and interaction between trt and timepoint
adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint + treatment_order + trt:timepoint,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = "margin"
)


# No significance, only keep fixed effects
# Now we can remove the two way interactions
adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint + treatment_order,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = "margin"
)

adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint + treatment_order,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = "margin"
)

pairwise.adonis(Fecal_data_micro.ps.dist, Fecal_data_micro.ps.df$timepoint, perm = 9999, p.adjust.m = "BH")

pairwise.adonis(Fecal_data_micro.ps.dist, Fecal_data_micro.ps.df$Combined_order, perm = 9999, p.adjust.m = "BH")


bd_order <- betadisper(Fecal_data_micro.ps.dist,
                       group = Fecal_data_micro.ps.df$timepoint,
                       )

permutest(bd_order, permutations = 9999)
TukeyHSD(bd_order)
plot(bd_order)


adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint + treatment_order,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = NULL
)

# no significance by treatment_order, just use two fixed effects
adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint ,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = NULL
)

adonis2(
  Fecal_data_micro.ps.dist ~ trt + timepoint ,
  data   = Fecal_data_micro.ps.df,
  strata = Fecal_data_micro.ps.df$participant_id,
  permutations = 9999,
  by = "margin"
)



