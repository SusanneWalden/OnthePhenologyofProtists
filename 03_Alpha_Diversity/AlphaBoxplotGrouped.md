Plot alpha diversity indices grouped by microhabitat in a boxplot
================

To compare different alpha diversity indices, we calculated species richness, Simpson index, Shannon index as well as eveness for each microhabitat per sampling period

Load data
---------

``` r
rm(list = ls())

library(vegan)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsignif)


#setwd("03_Alpha_Diversity")


OTU_Table = as.data.frame(read.csv2("../00_Data/05_Cercozoa_Seasonal_OTU_Table_min-freq-7633_transposed_withMetadata.csv",header = T))
SampleMetadata = OTU_Table[,1:17]
OTU_Table = OTU_Table[,18:ncol(OTU_Table)]
```

Calculate Alpha Diversity Indices
---------------------------------

To calculate the diversity indices, simply load the table and calculate the species richness per microhabitat and sampling period, then convert it into a dataframe and add all other diversity `indices` - in this case: The Simpson index, Shannon Index and eveness. Finally add the metadata and melt the dataframe.

``` r
richness = specnumber(OTU_Table)
df = as.data.frame(richness)

df$simpson = diversity(OTU_Table, index = "simpson")
df$shannon = diversity(OTU_Table, index = "shannon")
df$eveness = df$shannon/log(df$richness)

df$Microhabitat = SampleMetadata$Microhabitat
df$Season = SampleMetadata$Season


#df$TreeSpecies = SampleMetadata$TreeSpecies
#df$SamplingDescription = SampleMetadata$SamplingDescription
#df$MicrohabitatPerSampling = SampleMetadata$MicrohabitatPerSampling
#df$SamplingDescription = factor(df$SamplingDescription, levels = c('Autumn 2017','Spring 2018','Autumn 2018', 'Spring 2019' ))


df_melted = melt(df)
```

    ## Using Microhabitat, Season as id variables

Plot the Figure
---------------

Now we want to visualise the alpha diversity indices in a combined boxplot and calculate the differences:

``` r
df_melted_sort <- df_melted
df_melted_sort$Microhabitat <- factor(df_melted_sort$Microhabitat,      # Reordering group factor levels
                         levels = c("Fresh Leaves","Deadwood","Arboreal Soil","Bark","Hypnum","Lichen", 
                                     "Orthotrichum","Leaf Litter","Soil"))

g = ggplot(df_melted_sort, aes(x = Season, y = value, fill = Season)) + 
    stat_boxplot(geom = "errorbar", width = 0.1, show.legend = F) +
    geom_boxplot(show.legend = F) + 
    scale_fill_manual(values = alpha(c("#014636","#660033"), 0.8), 
                    limits = c("Spring","Autumn")) + 
    scale_x_discrete(limits = c("Spring","Autumn")) + 
    theme_minimal() + 
    labs(y = "Alpha Diversity", x = "Season")+ 
    theme(axis.text=element_text(size=14, face = "bold"), 
        axis.title=element_text(size=16, face = "bold"), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle=45, hjust = 1)) +
    facet_grid(variable ~ Microhabitat, scales = "free", switch = "y")+
    theme(plot.margin = margin(1,1,1,1, "cm")) +
    stat_compare_means(comparisons = list(c("Spring","Autumn")), label = "p.signif", 
                     size = 3.5, method = "wilcox.test", exact=FALSE, vjust = 0.1, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                       symbols = c("***", "**", "*","NS")))

g
```

![](AlphaBoxplotGrouped_files/figure-markdown_github/CercozoaAlphaBoxPlot-1.png)

``` r
#write.csv(df_melted_sort, file = "AlphaDiv_Values.csv")



#ggsave("AlphaDivGroup_300_336x252.jpeg", plot = g, 
#       device = "jpeg", dpi = 300, width = 336, height = 252, 
#       units = "mm")
#ggsave("AlphaDivGroup_600_336x252.pdf", plot = g, 
#       device = "pdf", dpi = 600, width = 336, height = 252, 
#       units = "mm")
#ggsave("AlphaDivGroup_600_336x252.tif", plot = g, 
#       device = "tiff", dpi = 600, width = 336, height = 252, 
#       units = "mm")
#ggsave("AlphaDivGroup_600_336x252.png", plot = g, 
#      device = "png", dpi = 600, width = 336, height = 252, 
#      units = "mm")
```
