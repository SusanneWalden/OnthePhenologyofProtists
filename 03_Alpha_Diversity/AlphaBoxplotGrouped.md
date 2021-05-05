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

![](AphaBoxplotGrouped_files/figure-markdown_github/CercozoaAlphaBoxPlot-1.png)

``` r
#write.csv(df_melted_sort, file = "AlphaDiv_Values.csv")



ggsave("AlphaDivGroup_300_336x252.jpeg", plot = g, 
       device = "jpeg", dpi = 300, width = 336, height = 252, 
       units = "mm")
```

    ## Warning in wilcox.test.default(c(267, 398, 236, 435, 277, 321, 358, 336, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(317, 334, 356, 356, 457, 323, 391, 355, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(318, 455, 417, 410, 380, 274, 308, 345, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(402, 404, 392, 407, 444, 367, 443, 369, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(433, 461, 457, 387, 301, 385, 464, 452, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(372, 381, 486, 394, 354, 361, 468, 386, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(364, 360, 478, 390, 404, 376, 425, 425, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(488, 313, 369, 450, 448, 450, 454, 374, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(505, 507, 511, 487, 520, 476, 515, 497, :
    ## cannot compute exact p-value with ties

``` r
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
