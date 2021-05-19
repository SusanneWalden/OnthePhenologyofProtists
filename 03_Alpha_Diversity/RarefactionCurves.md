Rarefaction Curves
================

A sufficient sequencing depth is essential for microbial community ecologists. We need to be sure that our sequencing approach (nearly) covers the whole diversity in a sampled habitat, to avoid unexplained variance as a result of shallow sequencing depth. Thus, calculation of rarefaction curves is an important tool check the overall sequencing success and the suitability of statistical tests (e.g. beta diversity).

Here, we use the **[iNEXT Package](http://chao.stat.nthu.edu.tw/wordpress/software_download/inext-online/)** by Anne Chao. It performs interpolation (rarefaction) as well as extrapolation. First we load the packages and data:

Load Data
---------

``` r
rm(list = ls())

library(iNEXT)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)


OTU_Table = as.data.frame(read.csv2("../00_Data/05_Cercozoa_Seasonal_OTU_Table_min-freq-7633_transposed_withMetadata.csv",header = T))

SampleMetadata = OTU_Table[,1:17]
species = OTU_Table[,18:ncol(OTU_Table)]
species_mat = as.matrix(species)
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


# Spring 2018 - aggregate by microhabitat
species_mat_aggregated_S18 = aggregate(species_mat_S18,
                                       by = list(SampleMetadata_S18$Microhabitat),
                                       FUN = "sum")
rownames(species_mat_aggregated_S18) = levels(list(SampleMetadata_S18$Microhabitat)[[1]])
species_mat_aggregated_S18 = species_mat_aggregated_S18[,-1]


# Spring 2019 - aggregate by microhabitat
#species_mat_aggregated_S19 = aggregate(species_mat_S19, 
#                                       by = list(SampleMetadata_S19$Microhabitat), 
#                                      FUN = "sum")
#rownames(species_mat_aggregated_S19) = levels(list(SampleMetadata_S19$Microhabitat)[[1]])
#species_mat_aggregated_S19 = species_mat_aggregated_S19[,-1]


# Autumn 2017 - aggregate by microhabitat
#species_mat_aggregated_A17 = aggregate(species_mat_A17, 
#                                       by = list(SampleMetadata_A17$Microhabitat), 
#                                       FUN = "sum")
#rownames(species_mat_aggregated_A17) = levels(list(SampleMetadata_A17$Microhabitat)[[1]])
#species_mat_aggregated_A17 = species_mat_aggregated_A17[,-1]

# Autumn 2018 - aggregate by microhabitat
#species_mat_aggregated_A18 = aggregate(species_mat_A18, 
#                                       by = list(SampleMetadata_A18$Microhabitat), 
#                                       FUN = "sum")
#rownames(species_mat_aggregated_A18) = levels(list(SampleMetadata_A18$Microhabitat)[[1]])
#species_mat_aggregated_A18 = species_mat_aggregated_A18[,-1]
```

Calculate rarefaction curves for every sampling period with the *iNEXT* function
--------------------------------------------------------------------------------

``` r
# Spring 2018
Cercozoa_out_S18 = iNEXT(t(species_mat_aggregated_S18), q = 0, #q = 0 means species richness
    datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
#save(Cercozoa_out_S18, file = "Cercozoa_out_S18.RData")
#load("Cercozoa_out_S18.RData")

# Spring 2019
#Cercozoa_out_S19 = iNEXT(t(species_mat_aggregated_S19), q = 0, #q = 0 means species richness
                         #datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
#save(Cercozoa_out_S19, file = "Cercozoa_out_S19.RData")
#load("Cercozoa_out_S19.RData")

# Autumn 2017
#Cercozoa_out_A17 = iNEXT(t(species_mat_aggregated_A17), q = 0, #q = 0 means species richness
                        # datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
#save(Cercozoa_out_A17, file = "Cercozoa_out_A17.RData")
#load("Cercozoa_out_H17.RData")

# Autumn 2018
#Cercozoa_out_A18 = iNEXT(t(species_mat_aggregated_A18), q = 0, #q = 0 means species richness
                         #datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
#save(Cercozoa_out_A18, file = "Cercozoa_out_A18.RData")
#load("Cercozoa_out_A18.RData")
```

Plot with ggplot2
-----------------

``` r
# plot
S18 = ggiNEXT(Cercozoa_out_S18, type = 1, color.var = "site") +
  scale_fill_manual(name = "Microhabitat", 
                    values = c("#006d2c","#fd8d3c","#fed976","#d4b9da","#045a8d","#0570b0","#a6bddb","#8e8878","#524640"), 
                    limits = c("Fresh Leaves","Bark","Arboreal Soil", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) + 
  scale_color_manual(name = "Microhabitat", 
                     values = c("#006d2c","#fd8d3c","#fed976","#d4b9da","#045a8d","#0570b0","#a6bddb","#8e8878","#524640"), 
                     limits = c("Fresh Leaves","Bark","Arboreal Soil", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) + 
  scale_shape_manual(name = "Microhabitat", 
                     values = rep(16, 9), 
                     limits = c("Fresh Leaves","Bark","Arboreal Soil", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) + 
  theme_minimal() +
  labs(title = "Spring 2018", 
       x = "Number of sequences", y = "Number of OTUs") +
  scale_x_continuous(labels = scales::comma) +
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16, face = "bold"), 
        legend.position = "right",
        legend.direction = "vertical", 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5), 
        legend.title.align = 0.5,
        axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 14)) +
  guides(fill = guide_legend(title="Microhabitat"), 
         color = guide_legend(title="Microhabitat"), 
         shape = guide_legend(title="Microhabitat"))

# we rearrange the layers to make them more visible
S18$layers = c(S18$layers[[3]], S18$layers[[2]], S18$layers[[1]])

S18
```

![](RarefactionCurves_files/figure-markdown_github/Plot%20rarefaction-1.png)

Combine all plots
=================

After you repeated the above mentioned code for every sampling period you can combine all plots with the *ggarrange* function from the ggpubr package.
