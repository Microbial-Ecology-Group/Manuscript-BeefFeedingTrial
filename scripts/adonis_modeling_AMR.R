
Fecal_AMR_data.ps <- subset_samples(AMR_data.ps, Sample_Type=="Fecal" & col_time != "Washout")

Fecal_AMR_data.ps <- prune_taxa(taxa_sums(Fecal_AMR_data.ps) > 0, Fecal_AMR_data.ps)


Fecal_AMR_data.ps.df <- as(sample_data(Fecal_AMR_data.ps),"data.frame")
Fecal_AMR_data.ps.dist <- vegdist(t(otu_table(Fecal_AMR_data.ps)), method = "bray")
Fecal_AMR_data.ps.ord <- vegan::metaMDS(comm = t(otu_table(Fecal_AMR_data.ps)), distance = "bray", try = 10, trymax = 100, autotransform = F)
# plt_Combined_order <- plot_ordination(Fecal_AMR_data.ps, Fecal_AMR_data.ps.ord, color = "Combined_order", shape = "Combined_order") +
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
Fecal_AMR_data.ps.adonis <- adonis2(Fecal_AMR_data.ps.dist ~Combined_order, data = Fecal_AMR_data.ps.df, permutations = 9999)
Fecal_AMR_data.ps.adonis #


#adonis2(Fecal_AMR_data.ps.dist ~ (trt_id * col_time): treatment_order, strata=Fecal_AMR_data.ps.df$Participant_ID , data = Fecal_AMR_data.ps.df, permutations = 9999)

#We tested whether or not there was a difference in rsponse depending on the order in which the subjects received the treatments (3 way interaction - not significant.)

# Full model, global test
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id * col_time * treatment_order,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = NULL
)
# Full model, by margin
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id * col_time * treatment_order,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = "margin"
)

# no significance in three way interaction
# now specify each term individually without the three way interaction
# Global test
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time + treatment_order + trt_id:col_time + trt_id:treatment_order + col_time:treatment_order,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = NULL
)
# By margin
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time + treatment_order + trt_id:col_time + trt_id:treatment_order + col_time:treatment_order,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = "margin"
)

# trt_id:treatment_order was nearly significant, reduce further to only include fixed effects and interaction between trt_id and col_time
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time + treatment_order + trt_id:col_time + trt_id:treatment_order,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = "margin"
)

# no Significance, but "trt_id:treatment_order" was close to significant (p = 0.07)
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time + treatment_order + treatment_order:trt_id,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = "margin"
)
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time + treatment_order + treatment_order:col_time,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = "margin"
)


# 
# Now we can remove the two way interactions
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time + treatment_order,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = NULL
)

adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time + treatment_order,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = "margin"
)

#  try to use just two fixed effects
adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time ,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = NULL
)

adonis2(
  Fecal_AMR_data.ps.dist ~ trt_id + col_time ,
  data   = Fecal_AMR_data.ps.df,
  strata = Fecal_AMR_data.ps.df$Participant_ID,
  permutations = 9999,
  by = "margin"
)



