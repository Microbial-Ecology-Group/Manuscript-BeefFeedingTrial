library(UpSetR)
#library(MicrobiotaProcess) #NEVER LOAD THIS



# Make new phyloseq object
beta_data_noSNP <- data_noSNP
beta_data_noSNP <- prune_taxa(taxa_sums(beta_data_noSNP) > 0, beta_data_noSNP)

# Fix random seed BEFORE normalization (for reproducibility)
set.seed(1234)

# Use CSS normalization (no need for transform_sample_counts)
#beta_data_noSNP.css <- phyloseq_transform_css(beta_data_noSNP, log = FALSE)

# Agglomerate by gene group
beta_ra_group <- tax_glom(beta_data_noSNP, taxrank = "group")
beta_ra_group.css <- phyloseq_transform_css(beta_ra_group, log = FALSE)

# Generate the upset matrix BEFORE any further scaling
upsetda_edit <- MicrobiotaProcess::get_upset(beta_ra_group, factorNames = "Combined_order")

# Continue your code from here
NEWupsetda_edit <- upsetda_edit
NEWupsetda_edit[] <- lapply(NEWupsetda_edit, function(x) {
  if (is.factor(x) || is.logical(x)) as.numeric(x) else x
})

# Rename columns to pretty publication labels
colnames(NEWupsetda_edit) <- c(
  "End No Claim Fecal",
  "End RWA Fecal",
  "Midpoint No Claim Fecal",
  "Midpoint RWA Fecal",
  "No Claim Beef",
  "RWA Beef",
  "Start No Claim Fecal",
  "Start RWA Fecal"
)

# Define plotting order
m_order <- c(
  "End RWA Fecal", "Midpoint RWA Fecal", "Start RWA Fecal",
  "End No Claim Fecal", "Midpoint No Claim Fecal", "Start No Claim Fecal",
  "RWA Beef", "No Claim Beef"
)

# Define color scheme
c_colors <- c(
  "No Claim Beef"           = "#FFD700",
  "RWA Beef"                = "#117744",
  "Start No Claim Fecal"    = "#660000",
  "Midpoint No Claim Fecal" = "#CC0066",
  "End No Claim Fecal"      = "#ff9999",
  "Start RWA Fecal"         = "#000999",
  "Midpoint RWA Fecal"      = "#006699",
  "End RWA Fecal"           = "#66CCFF"
)

# Plot (text scaled for paper)
AMR_Upset <- upset(
  NEWupsetda_edit,
  sets = m_order,
  sets.bar.color = rev(unname(c_colors)),
  text.scale = 3,
  order.by = "freq",
  empty.intersections = "on",
  keep.order = TRUE,
  mainbar.y.label = "Intersection Size (ARG Group)",
  sets.x.label = "Set Size (Samples)"
)
AMR_Upset

png("Figures/UpSet_AMR.png", width = 28, height = 20, units = "in", res = 300, bg = "white")
upset(
  NEWupsetda_edit,
  sets = m_order,
  sets.bar.color = rev(unname(c_colors)),
  text.scale = 4,
  order.by = "freq",
  empty.intersections = "on",
  keep.order = TRUE,
  mainbar.y.label = "Intersection Size (ARG Group)",
  sets.x.label = "Set Size (Samples)"
)
dev.off()










# Combine with ANCOM figure
final_panel <- plot_grid(
  ggdraw() + draw_image("AMR_UpSet_ARGs.png", scale = 1.0),
  ggdraw() + draw_image("Ancom_Fig.png", scale = 1.0),
  ncol = 1,
  rel_heights = c(0.8, 1.1),
  labels = c("A", "B"),
  label_size = 20,
  label_fontface = "bold",
  align = "v"
)
final_panel


ggsave("Panel_UpSet_ANCOMBC.png", final_panel,
       width = 15, height = 16, dpi = 300, bg = "white")

