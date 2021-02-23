

##########################################################################################################################################
##########################################################################################################################################
######################################################### Seasonnal Effect ###############################################################
##########################################################################################################################################
##########################################################################################################################################


rm(list = ls())

library(ggplot2)
library(vegan)
library(RColorBrewer)
library(reshape2)#for melting
library(agricolae)

#Seasonal Susi Cerco + Endo

setwd("C:/Users/Susi/Desktop/season_server/")


OTU_Table = as.data.frame(read.csv2("AllSeason_Cerc_En_Final_OTU_R_7633_edit_plottogether_extended.csv",header = T))

SampleMetadata = as.data.frame(OTU_Table[,1:17])
OTU_Table = OTU_Table[,18:ncol(OTU_Table)]



############ Calculate Alpha Diversity Indices
richness = specnumber(OTU_Table)
df = as.data.frame(richness)


df$simpson = diversity(OTU_Table, index = "simpson")
df$shannon = diversity(OTU_Table, index = "shannon")
df$eveness = df$shannon/log(df$richness)


rownames(df) = SampleMetadata$ï..SampleID # watch out: sometimes R adds "ï.." to the titel of the first column, probably due to .csv

############################## extend table with metadata

df$Microhabitat = SampleMetadata$Microhabitat
df$TreeSpecies = SampleMetadata$TreeSpecies
df$Season = SampleMetadata$Season
df$SamplingDescription = SampleMetadata$SamplingDescription
df$MicrohabitatPerSampling = SampleMetadata$MicrohabitatPerSampling
df$SamplingDescription = factor(df$SamplingDescription, levels = c('Autumn 2017','Spring 2018','Autumn 2018', 'Spring 2019' ))


######## melt table so you end up with only one row of alpha div values

df_melted = melt(df)



################################

g = ggplot(df_melted, aes(x = Microhabitat, y = value, fill = Microhabitat)) + 
  stat_boxplot(geom = "errorbar", width = 0.1, show.legend = F) +
  geom_boxplot(show.legend = F) + 
  scale_fill_manual(values = c("#138F6A","#FF944D","#FFE4B5","#BF6692","#A8ACFF","#4363d8","#8B8989","#9A6324","#800000"), 
                    limits = c("Fresh Leaves","Bark","Arboreal Soil", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) + 
  scale_x_discrete(limits = c("Fresh Leaves","Deadwood","Lichen","Bark","Orthotrichum","Hypnum","Arboreal Soil","Leaf Litter", "Soil")) + 
  #scale_fill_manual(values = c("#138F6A","#FF944D","#FFE4B5","#BF6692","#A8ACFF","#4363d8","#8B8989","#9A6324","#800000"), 
   #                 limits = c("Fresh Leaves","Bark","Arboreal Soil", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) + 
  #scale_x_discrete(limits = c("Fresh Leaves","Bark","Arboreal Soil", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) +  
  #scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  labs(y = "Alpha Diversity", 
       x = "Microhabitat")+ 
      
  theme(axis.text=element_text(size=14, face = "bold"), 
        axis.title=element_text(size=16, face = "bold"), 
        #plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 16, hjust = 0.5), 
        #legend.text = element_text(size = 10), 
       # legend.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle=45, hjust = 1))+
 facet_grid(variable ~ SamplingDescription, scales = "free", switch = "y")+
 theme(plot.margin = margin(1,1,1,1, "cm")) 

g
#write.csv(df_melted, file = "AlphaDiv_Values.csv")

###### export 
#write.csv(df_melted, file = "Melted_AlphaDiv.csv")
#write.csv(df, file = "Unmelted_AlphaDiv.csv")

#save

ggsave("AlphaDiv_300_336x252.jpeg", plot = g, 
       device = "jpeg", dpi = 300, width = 336, height = 252, 
       units = "mm")
ggsave("AlphaDiv_600_336x252.pdf", plot = g, 
       device = "pdf", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("AlphaDiv_600_336x252.tif", plot = g, 
       device = "tiff", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("AlphaDiv_600_336x252.png", plot = g, 
       device = "png", dpi = 600, width = 336, height = 252, 
       units = "mm")