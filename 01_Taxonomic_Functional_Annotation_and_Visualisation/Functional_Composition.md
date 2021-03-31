Functional Diversity
================

Functional diversity is a biodiversity measure based on functional traits of the species present in a community. Functional traits are those that define species in terms of their ecological roles - how they interact with the environment and with other species **[(Diaz and Cabido, 2001)](https://www.sciencedirect.com/science/article/pii/S0169534701022832?casa_token=t3ikAIDQjocAAAAA:sMQ2sSbZxMTJCQ70YjehZur0m9zd0MkHijOa8jax6uAE0OtAo4iLlAAxAr6GW_Q5ZqZiDNHyD35C)**.In this section we want to investigate differences in functional diversity, i.e. feeding mode of cercozaon genera detected from forest soil to the canopy.

Load data
---------

``` r
rm(list= ls())


library(ggplot2)
library(reshape2)
library(phenotypicForest)
library(magrittr)

TAX= read.csv2("../00_Data/AllSeason_Cerc_En_Taxonomy_R.csv",header = T)
OTU_Table = as.data.frame(read.csv2("../00_Data/05_Cercozoa_Seasonal_OTU_Table_min-freq-7633_transposed_withMetadata.csv",header = T))

SampleMetadata = OTU_Table[,1:17]
Species = OTU_Table[,18:ncol(OTU_Table)]
```

Prepare Data
------------

``` r
# group Sample Metadata by sampling period

SampleMetadata_S18 = SampleMetadata[OTU_Table$SamplingDescription == "Spring 2018",]
SampleMetadata_S19 = SampleMetadata[OTU_Table$SamplingDescription == "Spring 2019",]
SampleMetadata_A17 = SampleMetadata[OTU_Table$SamplingDescription == "Autumn 2017",]
SampleMetadata_A18 = SampleMetadata[OTU_Table$SamplingDescription == "Autumn 2018",]

# group OTUS by sampling period

species_mat_S18 = Species[OTU_Table$SamplingDescription == "Spring 2018",]
species_mat_S19 = Species[OTU_Table$SamplingDescription == "Spring 2019",]
species_mat_A17 = Species[OTU_Table$SamplingDescription == "Autumn 2017",]
species_mat_A18 = Species[OTU_Table$SamplingDescription == "Autumn 2018",]


Abundances = colSums(species_mat_S18)# choose subtable
TAX = cbind(TAX, Abundances)

TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)
TAX$Nutrition = as.factor(TAX$Nutrition)
```

Aggregate tables
----------------

``` r
NutritionTable = species_mat_S18# choose subtable
colnames(NutritionTable) = TAX$Nutrition

# First aggregate by Habitat
HabitatAggregatedNutritionTable = aggregate(NutritionTable, by = list(SampleMetadata_S18$Microhabitat), FUN = sum)# choose subtable
rownames(HabitatAggregatedNutritionTable) = HabitatAggregatedNutritionTable$Group.1
HabitatAggregatedNutritionTable = HabitatAggregatedNutritionTable[,-1]


# if you want to check the number of OTUs instead of the amount of reads,this is the point to convert it into a presence/absence matrix like this:
#HabitatAggregatedNutritionTable[,] = ifelse(HabitatAggregatedNutritionTable[,] > 0, 1, 0)

#Then aggregate the (transposed) table by Nutrition
HabitatAggregatedNutritionTable = aggregate(t(HabitatAggregatedNutritionTable), 
                                            by = list(TAX$Nutrition), 
                                            FUN = sum)
rownames(HabitatAggregatedNutritionTable) = HabitatAggregatedNutritionTable$Group.1
HabitatAggregatedNutritionTable = HabitatAggregatedNutritionTable[,-1]
HabitatAggregatedNutritionTable = as.data.frame(HabitatAggregatedNutritionTable)

# melt the table so ggplot can interpret it as a histogram
data = HabitatAggregatedNutritionTable
data$Nutrition = rownames(data)
data2 = melt(data, id.vars = "Nutrition")

# add the Stratum
data3 = data2
data3$Stratum = ifelse(data2$variable %in% c("Leaf Litter", "Soil"), "Ground", "Canopy")

# for the plotting we need to change the column names
colnames(data3) = c("score", "item", "value", "family")
data3$item <- factor(data3$item, levels=c("Soil","Leaf Litter","Lichen","Orthotrichum","Hypnum", "Arboreal Soil","Bark","Deadwood","Fresh Leaves"))
```

ggplot
------

``` r
#plot 

S18 =ggplot(data3, aes(fill=score, y=value, x=item)) + 
  scale_fill_manual(values = c("#4E9258","#7D053D", "#7E587E", 
                               "coral2","#A9A9A9","#348781","#d3a540"), 
                    limits = c("autotroph","bacterivore", "eukaryvore", 
                               "not plant parasite","omnivore","plant parasite", "undetermined")) +
  geom_bar(position="fill", stat="identity") + 
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  labs(fill = "Nutrition", y= "Relative Abundances", x = "Microhabitat",title = "Spring 2018")+
  theme_minimal()+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=14), 
        axis.title=element_text(size=16,face="bold"))+
  geom_label(x = 6, y = 0.7, fill = "darkolivegreen4", color = "white", 
             label = "Canopy", size = 8) +
  geom_label(x = 1.5, y = 0.7, fill = "burlywood4", color = "white", 
             label = "Ground", size = 8)+
  geom_vline(xintercept = 2.5, size = 0.6, linetype=1.5, color = "#6A785F")+
  coord_flip()

S18
```

![](Functional_Composition_files/figure-markdown_github/Plot%20RDA-1.png)

Combine all plots
-----------------

After you repeated the above mentioned code for every sampling period you can combine all plots with the *ggarrange* function from the ggpubr package.

``` r
#combi = ggarrange( H17, S18, A18, S19,
#                   labels = c("A", "B","C","D"), 
#                   ncol = 2, nrow = 2, font.label = list(size = 16, color = "black"), 
#                   common.legend = T, legend = "right", 
#                   align = "v")+
#          theme(plot.margin = margin(1,1,1,1, "cm")) 
#combi


#save
#ggsave("NutritionCombined_300_336x252.jpeg", plot = combi, 
#      device = "jpeg", dpi = 300, width = 336, height = 252, 
#     units = "mm")
#ggsave("NutritionCombined_600_336x252.pdf", plot = combi, 
#       device = "pdf", dpi = 600, width = 336, height = 252, 
#       units = "mm")
#ggsave("NutritionCombined_600_336x252.tif", plot = combi, 
#       device = "tiff", dpi = 600, width = 336, height = 252, 
#       units = "mm")
#ggsave("NutritionCombined_600_336x252.png", plot = combi, 
#       device = "png", dpi = 600, width = 336, height = 252, 
#       units = "mm")
```
