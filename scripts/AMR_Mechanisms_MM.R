#Source all AMR data


#Identify top four abundant classes
ggplot(joined_ra_class_ps_class_melt, aes(x = reorder(class, -Abundance), y = Abundance, fill = class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "class", y = "Relative Abundance", title = "Top 4 Most Abundant Classes") +
  theme_minimal() +
  scale_x_discrete(limits = unique(joined_ra_class_ps_class_melt$class)[1:6])
#Tetracyclines, betalactams, DB resistance, MLS
#According to the figure, MLS has more abundance than DB, but if you change it to top 3, it counts
#Db instead of mls. 



unique_classes <- unique(joined_ra_class_ps_class_melt$class)
print(unique_classes)

# Get unique classes
unique_classes <- unique(joined_ra_class_ps_class_melt$class)

# Create an empty dataframe to store the results
class_abundance <- data.frame(class = character(), total_abundance = numeric(), percentage = numeric(), stringsAsFactors = FALSE)

# Loop through each unique class
for (class in unique_classes) {
  # Subset the data for the current class
  class_data <- joined_ra_class_ps_class_melt[joined_ra_class_ps_class_melt$class == class, ]
  # Calculate the total abundance for the current class
  total_abundance <- sum(class_data$Abundance)
  # Calculate the percentage of abundance for the current class
  percentage <- (total_abundance / sum(joined_ra_class_ps_class_melt$Abundance)) * 100
  # Append the results to the dataframe
  class_abundance <- rbind(class_abundance, data.frame(class = class, total_abundance = total_abundance, percentage = percentage))
}

# Print the result
print(class_abundance)











#Create a color palette
my_colors <- c("#990033", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF", "#FFFF00", "#800080", "#008000", "#000080",
               "#FFA500", "#663333", "#008080", "#FF0000", "#808000", "#FFC0CB", "#99FF99", "#808080", "#C0C0C0",
               "#FF6347", "#993300", "#00CED1", "#6666FF", "#DC143C", "#00CC66", "#4B0082", "#FF8C00", "#9932CC",
               "#8B008B", "#ADFF2F", "#CC9966", "#000099","#CC0099","#66FF66")
num_colors <- 30
my_palette <- colorRampPalette(my_colors)(30)
#Changed the first hexadecimal number from red to maroon-ish

#Tetracycline Resistance by Group
tet_class_melted.css <- subset_taxa(beta_data_noSNP.css, class == "Tetracyclines") %>%
  psmelt() 

p1 <- ggplot(tet_class_melted.css, aes(x = Sample_Type, y = Abundance, fill = group)) +
  theme_bw() +
  labs(x = "", y = "Abundance") +
  scale_fill_manual(values = my_colors) +
  geom_bar(stat = "summary", colour = "black") +  
  #ggtitle("Tetracycline Resistance") +
  theme(legend.position = "right",
        plot.margin = unit(c(0.25, 4, 0.25, 1), "cm"), 
        panel.border = element_rect(colour = "black", size = 1),  # Normal thickness
        axis.ticks = element_line(size = 0.8, colour = "black"),  # Normal thickness
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 18, vjust = 2.5),
        axis.text.x = element_text(size = 18, colour = "black", vjust = 1, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 0.01)  
p1









p1_1 <- ggplot(tet_class_melted.css, aes(x = Sample_Type, y = Abundance, fill = mechanism)) +
  theme_bw() +
  #labs(y = "Tetracycline Resistance Mechanisms") +
  geom_bar(stat = "summary", colour = "black") +
  #ggtitle("Resistance Mechanisms Across Antibiotic Classes") +
  scale_fill_manual(values = my_colors) +  
  theme(legend.position = "right",
        plot.margin = unit(c(0.25, 4, 0.25, 1), "cm"),  
        panel.border = element_rect(colour = "black", size = 1),  
        axis.ticks = element_line(size = 0.8, colour = "black"),  
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 18, vjust = 2.5),
        axis.text.x = element_text(size = 18, colour = "black", vjust = 1, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 0.1)  
p1_1


#Plot final grid
tet_final_plot <- plot_grid(
  ggdraw() + draw_label("Tetracycline Resistance", size = 18, fontface = "plain"),
  plot_grid(p1, p1_1, ncol = 2, align = "hv", rel_widths = c(1, 1)),  # Aligns horizontally & vertically
  ncol = 1,
  rel_heights = c(0.15, 1)  # Slightly smaller title space
)
tet_final_plot







#MLS Resistance by Group
MLS_class_melted.css <- subset_taxa(beta_data_noSNP.css, class == "MLS") %>%
  psmelt() 
p2 <- ggplot(MLS_class_melted.css, aes(x= sample_type, y= Abundance, fill = group)) +
  theme_bw() +labs(y= "MLS Resistance Groups") +
  #facet_grid(~sample_type) + 
  geom_bar(stat = "summary", colour = "black")+
  scale_fill_manual(values = my_colors) +
  coord_fixed(ratio = 0.04)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p2

#MLS Resistance by Mechanism
MLS_class_melted.css <- subset_taxa(beta_data_noSNP.css, class == "MLS") %>%
  psmelt() 
p2_2 <- ggplot(MLS_class_melted.css, aes(x= sample_type, y= Abundance, fill = mechanism)) +
  theme_bw() +
  labs(y= "MLS Resistance Mechanisms") +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = my_colors) +
  theme(
    plot.margin = unit(c(0.25, 7.25, 0.25, 7.25), "cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 20),
    axis.title.y = element_text(size = 20, vjust = 2.5),
    axis.text.x = element_text(size = 10, colour = "black", vjust = 1, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 0.1)
p2_2

#Betalactam Resistance by Group
beta_class_melted.css <- subset_taxa(beta_data_noSNP.css, class == "betalactams") %>%
  psmelt() 
p3 <- ggplot(beta_class_melted.css, aes(x= sample_type, y= Abundance, fill = group)) +
  theme_bw() +labs(y= "Betalactam Resistance Genes") +
  #facet_grid(~sample_type) + 
  geom_bar(stat = "summary", colour = "black")+
  scale_fill_manual(values = my_colors) +
  coord_fixed(ratio = 0.06)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p3

#Betalactam Resistance by Mechanism
beta_class_melted.css <- subset_taxa(beta_data_noSNP.css, class == "betalactams") %>%
  psmelt() 
p3_3 <- ggplot(beta_class_melted.css, aes(x= sample_type, y= Abundance, fill = mechanism)) +
  theme_bw() +labs(y= "Betalactam Resistance Mechanisms")+
  #facet_grid(~sample_type) + 
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = my_palette) +
  #geom_text(aes(label = sample_type), vjust = 1, color = "black", size = 3) +
  theme(legend.position = "right",
    plot.margin = unit(c(0.25,7.25,0.25,7.25),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 20),
    axis.title.y = element_text(size = 20, vjust = 2.5),
    axis.text.x = element_text(size = 10, colour = "black", vjust = 1, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank (),
    panel.grid.major.x = element_blank())+
    coord_fixed(ratio = 0.7)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p3_3

#Drug and biocide Resistance by Group
#db_class_melted.css <- subset_taxa(beta_data_noSNP.css, class == "Drug and biocide resistance") %>%
psmelt() 

#ggplot(db_class_melted.css, aes(x= sample_type, y= Abundance, fill = group)) +
theme_bw() +labs(y= "Drug and biocide Resistance Mechanisms") +
  #facet_grid(~sample_type) + 
  geom_bar(stat = "summary", colour = "black")


#Drug and biocide Resistance by Mechanism
#db_class_melted.css <- subset_taxa(beta_data_noSNP.css, class == "Drug and biocide resistance") %>%
psmelt() 
#ggplot(db_class_melted.css, aes(x= sample_type, y= Abundance, fill = mechanism)) +
theme_bw() +labs(y= "Drug and Biocide Resistance Mechanisms")+
  #facet_grid(~sample_type) + 
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = my_palette) +
  #geom_text(aes(label = sample_type), vjust = 1, color = "black", size = 3) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.25,7.25,0.25,7.25),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 20),
    axis.title.y = element_text(size = 20, vjust = 2.5),
    axis.text.x = element_text(size = 10, colour = "black", vjust = 1, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank (),
    panel.grid.major.x = element_blank())




#Now use cowplot to make a panel containing all of the figures
library(cowplot)

#First, label plots how you want them to be laid out in the panel
#ReLabel Figures from above
#I am not including DB because it has so many groups and looks awful
#Make sure your first figure has labels appropriate for the panel



#Below panel works, but is messy. Probably be best to make this panel in BioRender

#Plot final grid
final_plot <- plot_grid(
  plot_grid(p1, p1_1, labels = c("a", "b"), ncol = 2),
  plot_grid(p2, p2_2, labels = c("c", "d"), ncol = 2),
  plot_grid(p3, p3_3, labels = c("e","f"), ncol = 2),
  nrow = 3,
  align = "hv")
final_plot

