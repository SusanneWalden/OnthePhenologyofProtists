rm(list = ls())

library(iNEXT)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)




#AllSeason
setwd("C:/Users/Susi/Desktop/season_server/Functional_Traits/F19/")
OTU_Table = as.data.frame(read.csv2("AllSeason_F19_Cerc_En_Final_OTU_R_7633.csv",header = T))

SampleMetadata = OTU_Table[,1:5]
rownames(OTU_Table) = SampleMetadata$SampleID
rownames(SampleMetadata) = SampleMetadata$SampleID
species = OTU_Table[,6:ncol(OTU_Table)]
species_mat = as.matrix(species)

#We want to test the rarefaction for the microhabitats, so now we aggregate the OTU table
species_mat_aggregated = aggregate(species_mat, by = list(SampleMetadata$Microhabitat), FUN = "sum")
rownames(species_mat_aggregated) = levels(list(SampleMetadata$Microhabitat)[[1]])
species_mat_aggregated = species_mat_aggregated[,-1]





#H17
Cercozoa_out_A = iNEXT(t(species_mat_aggregated), q = 0, #q = 0 means species richness
    datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
save(Cercozoa_out_A, file = "Cercozoa_out_A.RData")
#load("Cercozoa_out_A.RData")
#F18
Cercozoa_out_B = iNEXT(t(species_mat_aggregated), q = 0, #q = 0 means species richness
                     datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
save(Cercozoa_out_B, file = "Cercozoa_out_B.RData")
#load("Cercozoa_out_B.RData")
#H18
Cercozoa_out_C = iNEXT(t(species_mat_aggregated), q = 0, #q = 0 means species richness
                       datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
save(Cercozoa_out_C, file = "Cercozoa_out_C.RData")
#load("Cercozoa_out_C.RData")
#F19
Cercozoa_out_D = iNEXT(t(species_mat_aggregated), q = 0, #q = 0 means species richness
                       datatype = "abundance", nboot = 99, conf = 0.97, knots = 250, endpoint = NULL)
save(Cercozoa_out_D, file = "Cercozoa_out_D.RData")
#load("Cercozoa_out_D.RData")

#Cercozoa_out_A
#Cercozoa_out_B
#Cercozoa_out_C
#Cercozoa_out_D



################### Plot



g = ggiNEXT(Cercozoa_out_D, type = 1, color.var = "site") +
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
  labs(title = "Spring 2019", 
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

# In the default plot, brighter colors are hard to see.
# So here we rearrange the layers to make them more visible
g$layers = c(g$layers[[3]], g$layers[[2]], g$layers[[1]])

g

A=g
B=g
C=g
D=g


############plot together 
combi = ggarrange( A,B,C,D,
                   labels = c("A", "B","C","D"), 
                   ncol = 2, nrow = 2, font.label = list(size = 16, color = "black"), 
                   common.legend = T, legend = "right", 
                   align = "v")+
  theme(plot.margin = margin(1,1,1,1, "cm")) 

combi

#############save
ggsave("Rarefaction_300_336x252.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 336, height = 252, 
       units = "mm")
ggsave("Rarefaction_600_336x252.pdf", plot = combi, 
       device = "pdf", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("Rarefaction_600_336x252.tif", plot = combi, 
       device = "tiff", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("Rarefaction_600_336x252.png", plot = combi, 
       device = "png", dpi = 600, width = 336, height = 252, 
       units = "mm")


