theme_set(theme_bw(base_size = 30))

# ---- Subset and Prepare Data ----

amr_ps_filt <- subset_samples(AMR_data.ps, Combined_order != "Final_Washout")
tet_ps <- subset_taxa(amr_ps_filt, class == "Tetracyclines")

tet_mech_ps <- tax_glom(tet_ps, taxrank = "mechanism")
tet_group_ps <- tax_glom(tet_ps, taxrank = "group")

tet_mech_melt <- psmelt(tet_mech_ps)
tet_group_melt <- psmelt(tet_group_ps)

tet_mech_melt <- subset(tet_mech_melt, !is.na(mechanism) & mechanism != "" &
                          !is.na(Abundance) & Abundance > 0 &
                          !is.na(Trial) & Trial != "")
tet_group_melt <- subset(tet_group_melt, !is.na(group) & group != "" &
                           !is.na(Abundance) & Abundance > 0 &
                           !is.na(Trial) & Trial != "")

tet_mech_melt$mechanism <- factor(tet_mech_melt$mechanism)

sample_order <- c("No Claim Beef", "RWA Beef", "No Claim Fecal", "RWA Fecal")
tet_mech_melt$Trial <- factor(tet_mech_melt$Trial, levels = sample_order)
tet_group_melt$Trial <- factor(tet_group_melt$Trial, levels = sample_order)

mm_colors <- c(
  "#E41A1C", "#0066CC", "#4DAF4A", "#984EA3", "#FF7F00", "#33FF99",
  "#A65628", "#F781BF", "#cc0066", "#66C2A5", "#006600", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#3399ff", "#1B9E77",
  "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#ffff00",
  "#660066", "#FF33CC", "#33FFCC", "#CC33FF", "#FF3366", "#33CCFF",
  "#FF9933", "#993399", "#FF6666", "#9999FF", "#66FF33", "#FF33FF"
)

# ---- Figure A: Mechanism ----
p_mech <- ggplot(tet_mech_melt, aes(x = Trial, y = Abundance, fill = mechanism)) +
  geom_bar(stat = "summary", fun = "mean", color = "black") +
  labs(x = "", y = "Abundance", title = "Tetracycline Resistance by Mechanism") +
  scale_fill_manual(
    values = mm_colors,
    labels = gsub("(?i)tetracycline[ _-]*", "", levels(tet_mech_melt$mechanism), perl = TRUE)
  ) +
  guides(
    fill = guide_legend(
      nrow = 4, byrow = TRUE,
      keywidth = 2.3, keyheight = 2.3
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(size = 42, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 32),
    legend.box.margin = margin(t = 18, b = 18),
    axis.title.y = element_text(size = 40, vjust = 3),
    axis.text.x = element_text(size = 30, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 30, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---- Figure B: Gene Group ----
p_group <- ggplot(tet_group_melt, aes(x = Trial, y = Abundance, fill = group)) +
  geom_bar(stat = "summary", fun = "mean", color = "black") +
  labs(x = "", y = "Abundance", title = "Tetracycline Resistance by Gene Group") +
  scale_fill_manual(values = mm_colors) +
  guides(
    fill = guide_legend(
      nrow = 5, byrow = TRUE,
      keywidth = 2.3, keyheight = 2.3
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(size = 42, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 32),
    legend.box.margin = margin(t = 18, b = 18),
    axis.title.y = element_text(size = 40, vjust = 3),
    axis.text.x = element_text(size = 30, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 30, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---- Final Panel Layout ----
tet_panel <- plot_grid(
  ggdraw() + draw_label("Tetracycline Resistance", size = 46, fontface = "bold", hjust = 0.5),
  plot_grid(
    p_mech,
    NULL,
    p_group,
    labels = c("A", "", "B"),
    label_size = 36,
    label_fontface = "bold",
    ncol = 3,
    rel_widths = c(1, 0.08, 1)
  ),
  ncol = 1,
  rel_heights = c(0.15, 1.5)
)

tet_panel <- ggdraw() +
  draw_plot(tet_panel, x = 0.04, y = 0, width = 0.92, height = 1)


# ---- Save ----
ggsave("Figures/Tetracycline_Resistance_Mechanism_Group_Panel.png", tet_panel, width = 34, height = 28, dpi = 300, bg = "white")


