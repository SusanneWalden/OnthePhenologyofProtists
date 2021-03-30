Redundancy Analysis
================

The Redundancy Analysis (RDA) is a direct gradient analysis technique which summarises linear relationships between components of response variables that are "redundant" with (i.e. "explained" by) a set of explanatory variables. To test the influence of explanatory variables i.e. microhabitats and tree species - on response variables (our OTU abundances) the following code for a RDA was used.

Load data
---------

``` r
rm(list = ls())

library(ggplot2)
library(vegan)
library(RColorBrewer)
library(funrar)
library(ggpubr)

OTU_Table = as.data.frame(read.csv2("../00_Data/05_Cercozoa_Seasonal_OTU_Table_min-freq-7633_transposed_withMetadata.csv",header = T,fileEncoding="UTF-8-BOM"))

SampleMetadata = OTU_Table[,1:17]
rownames(OTU_Table) = SampleMetadata$SampleID
rownames(SampleMetadata) = SampleMetadata$SampleID

species = OTU_Table[,18:ncol(OTU_Table)]
species_mat = as.matrix(species)
species_mat = make_relative(species_mat)
```

Data Preparation
----------------

``` r
# group Sample Metadata by sampling period

SampleMetadata_S18 = SampleMetadata[OTU_Table$SamplingDescription == "Spring 2018",]
SampleMetadata_S19 = SampleMetadata[OTU_Table$SamplingDescription == "Spring 2019",]
SampleMetadata_A17 = SampleMetadata[OTU_Table$SamplingDescription == "Autumn 2017",]
SampleMetadata_A18 = SampleMetadata[OTU_Table$SamplingDescription == "Autumn 2018",]

# group OTUS by sampling period

species_mat_S18 = species_mat[OTU_Table$SamplingDescription == "Spring 2018",]
species_mat_S19 = species_mat[OTU_Table$SamplingDescription == "Spring 2019",]
species_mat_A17 = species_mat[OTU_Table$SamplingDescription == "Autumn 2017",]
species_mat_A18 = species_mat[OTU_Table$SamplingDescription == "Autumn 2018",]
```

Calculate RDA
-------------

Next we normalize our species matrix with a hellinger transformation and perform the redundancy analysis for every sampling period:

``` r
#Spring 2018
species_mat_he = decostand(species_mat_S18, "hellinger")# choose the subtable depending on sampling period

Microhabitat = SampleMetadata_S18$Microhabitat
TreeSpecies = SampleMetadata_S18$TreeSpecies

species_mat_rda = rda(species_mat_he~Microhabitat+TreeSpecies)

summary_species_mat_rda = summary(species_mat_rda, scaling = 2)#scaling = 2 represents a linear correlation between the vectors.

#Spring 2019
#species_mat_he = decostand(species_mat_S19, "hellinger")
#Microhabitat = SampleMetadata_S19$Microhabitat
#TreeSpecies = SampleMetadata_S19$TreeSpecies
#species_mat_rda = rda(species_mat_he~Microhabitat+TreeSpecies)
#summary_species_mat_rda_S19 = summary(species_mat_rda, scaling = 2)

#Autumn 2017
#species_mat_he = decostand(species_mat_A17, "hellinger")
#Microhabitat = SampleMetadata_A17$Microhabitat
#TreeSpecies = SampleMetadata_A17$TreeSpecies
#species_mat_rda = rda(species_mat_he~Microhabitat+TreeSpecies)
#summary_species_mat_rda_A17 = summary(species_mat_rda, scaling = 2)

#Autumn 2018
#species_mat_he = decostand(species_mat_A18, "hellinger")
#Microhabitat = SampleMetadata_A18$Microhabitat
#TreeSpecies = SampleMetadata_A18$TreeSpecies
#species_mat_rda = rda(species_mat_he~Microhabitat+TreeSpecies)
#summary_species_mat_rda_A18 = summary(species_mat_rda, scaling = 2)
```

Create Dataframes
-----------------

From the RDA-summary object, we need to extract three dataframes:

1.  The positions for the microhabitat vectors
2.  The positions fot the OTUs
3.  The positions for the tree species vectors

``` r
# For the microhabitats (rows 1-9)
df1  <- data.frame(summary_species_mat_rda$centroids[1:9,1:2])
rownames(df1) = gsub("Microhabitat", "", rownames(df1))
df1$Microhabitat = rownames(df1)

# For the species (OTUs)
df2  <- data.frame(summary_species_mat_rda$species[,1:2])

# For the tree species (rows 10-12)
df3 = data.frame(summary_species_mat_rda$centroids[10:12,1:2])
rownames(df3) = gsub("TreeSpecies", "", rownames(df3))
df3$TreeSpecies = rownames(df3)
```

Choose you colors wisely by creating a colour vector
----------------------------------------------------

``` r
ColorVectorTreeSpecies = vector()
ColorVector2 = vector()

# Again, these colors correspond to the YlGn Palette
for (a in rownames(df1)){
  if(a == "Arboreal Soil"){ColorVector2 = c(ColorVector2, "#FFE4B5")}
  else if (a == "Bark"){ColorVector2 = c(ColorVector2, "#FF944D")}
  else if (a == "Deadwood"){ColorVector2 = c(ColorVector2, "#BF6692")}
  else if (a == "Fresh Leaves"){ColorVector2 = c(ColorVector2, "#138F6A")}
  else if (a == "Hypnum"){ColorVector2 = c(ColorVector2, "#A8ACFF")}
  else if (a == "Lichen"){ColorVector2 = c(ColorVector2, "#8B8989")}
  else if (a == "Orthotrichum"){ColorVector2 = c(ColorVector2, "#4363d8")}
  else if (a == "Leaf Litter"){ColorVector2 = c(ColorVector2, "#9A6324")}
  else if (a == "Soil"){ColorVector2 = c(ColorVector2, "#800000")}
}


for (a in rownames(df3)){
  if(a == "Tilia cordata"){ColorVectorTreeSpecies = c(ColorVectorTreeSpecies, "#6b4e71")}
  else if (a == "Quercus robur"){ColorVectorTreeSpecies = c(ColorVectorTreeSpecies, "#53687e")}
  else if (a == "Fraxinus excelsior"){ColorVectorTreeSpecies = c(ColorVectorTreeSpecies, "#c2b2b4")}
}
```

Plot with ggplot2
-----------------

``` r
# plot
S18 = ggplot(df2, aes(x=RDA1, y=RDA2)) + 
  geom_point(aes(fill = "OTUs"), size=3, color = "darkslategrey") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  geom_segment(data=df1, aes(x=0, xend=RDA1, y=0, yend=RDA2, 
                             linetype = "Microhabitat", color = Microhabitat), 
               arrow=arrow(length=unit(0.01,"npc")), 
               size = 1.5, alpha = 0.9) +
  scale_color_manual(values = unique(ColorVector2), 
                     breaks = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) +
  geom_label(data=df1, 
             aes(x=RDA1,y=RDA2,label=rownames(df1),
                 hjust=0.5*(1-sign(RDA1)),
                 vjust=0.5*(1-sign(RDA2))), 
             fill=ColorVector2, color = "black", size=4) + 
  geom_segment(data=df3, aes(x=0, xend=RDA1, y=0, yend=RDA2, 
                             linetype = "Tree Species", color = TreeSpecies), 
               arrow=arrow(length=unit(0.01,"npc")), 
               size = 1.5, alpha = 0.75, 
               color = ColorVectorTreeSpecies) +
  geom_label(data=df3, 
             aes(x=RDA1,y=RDA2,label=rownames(df3),
                 hjust=0.5*(1-sign(RDA1)),
                 vjust=0.5*(1-sign(RDA2))), 
             fill=ColorVectorTreeSpecies, color = "black", size=4, fontface="italic") + 
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_continuous(limits = c(-1, 1)) +
  labs(title = "Spring 2018", 
       x = paste0("RDA1\nVariance explained: ",
                  round(as.numeric(summary_species_mat_rda$cont$importance["Proportion Explained", "RDA1"]), digits = 4) * 100, 
                  "%"), 
       y = paste0("RDA2\nVariance explained: ",
                  round(as.numeric(summary_species_mat_rda$cont$importance["Proportion Explained", "RDA2"]), digits = 4) * 100, 
                  "%")) +
  theme_minimal() + 
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16, face = "bold"), 
        legend.position = "right",
        legend.direction = "vertical", 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5), 
        legend.title.align = 0.5,
        axis.title = element_text(size =16, face = "bold"), 
        axis.text = element_text(size = 14)) +
  scale_linetype_manual(name = "Explanatory Variable", 
                        values = c("Microhabitat"=1, "Tree Species"=2)) + 
  scale_fill_manual(name = "Response Variable", 
                    values = "darkslategrey") +
  guides(color = F, 
         linetype=guide_legend(keywidth = 3, keyheight = 1, order = 1), 
         fill = guide_legend(order = 2))
S18
```

![](RDA_files/figure-markdown_github/Plot%20RDA-1.png)

Combine all plots
=================

After you repeated the above mentioned code for every sampling period you can combine all plots with the *ggarrange* function from the ggpubr package.

``` r
# combine all plots

#combi = ggarrange( AH17, S18, A18, S19,
#                   labels = c("A", "B","C","D"), 
#                   ncol = 2, nrow = 2, font.label = list(size = 16, color = "black"), 
#                   common.legend = T, legend = "right", 
#                   align = "v")+
#        theme(plot.margin = margin(1,1,1,1, "cm")) 

#combi


# save plot 
#ggsave("Rarefaction_300_336x252.jpeg", plot = combi, 
#       device = "jpeg", dpi = 300, width = 336, height = 252, 
#       units = "mm")
#ggsave("Rarefaction_600_336x252.pdf", plot = combi, 
#       device = "pdf", dpi = 600, width = 336, height = 252, 
#       units = "mm")
#ggsave("Rarefaction_600_336x252.tif", plot = combi, 
#       device = "tiff", dpi = 600, width = 336, height = 252, 
#       units = "mm")
#ggsave("Rarefaction_600_336x252.png", plot = combi, 
#       device = "png", dpi = 600, width = 336, height = 252, 
#       units = "mm")
```
