
rm(list = ls())

library(ggplot2)
library(vegan)
library(RColorBrewer)
library(reshape2)#for melting
library(agricolae)
library(cowplot)

#Seasonal Susi Cerco + Endo

setwd("C:/Users/Susi/Desktop/season_server/Submission/Supplementary/")


OTU_Table = as.data.frame(read.csv2("Boxplot_Parasites.csv",header = T))

names(OTU_Table)[names(OTU_Table) == "ï..Microhabitat"] <- "Microhabitat"




#TukeyLetters = HSD.test(aov(OTU_Table$Reads ~ OTU_Table$Microhabitat),
                       # "OTU_Table$Microhabitat")$groups
#TukeyLetters$Microhabitat = rownames(TukeyLetters)


#SampleMetadata = as.data.frame(OTU_Table[,1])
#Species = OTU_Table[,2:ncol(OTU_Table)]

# check col names
colnames(OTU_Table)
#colnames(Species)
#colnames(SampleMetadata)
#names(SampleMetadata)[names(SampleMetadata) == "OTU_Table[, 1]"] <- "Microhabitat"
OTU.274 <- OTU_Table[ which(OTU_Table$OTUs == 'OTU 274 Polymyxa betae'),]
OTU.230 <- OTU_Table[ which(OTU_Table$OTUs == 'OTU 230 Spongospora nasturtii'),]
OTU.566 <- OTU_Table[ which(OTU_Table$OTUs == 'OTU 566 unassigned Leptophryidae'),]



#HSD by microhabitat
#OTU274
TukeyLetters_274 = HSD.test(aov(OTU.274$Reads ~ OTU.274$Microhabitat),
                        "OTU.274$Microhabitat")$group
TukeyLetters_274$Microhabitat = rownames(TukeyLetters_274)
#OTU230
TukeyLetters_230 = HSD.test(aov(OTU.230$Reads ~ OTU.230$Microhabitat),
                            "OTU.230$Microhabitat")$group
TukeyLetters_230$Microhabitat = rownames(TukeyLetters_230)
#OTU566
TukeyLetters_566 = HSD.test(aov(OTU.566$Reads ~ OTU.566$Microhabitat),
                            "OTU.566$Microhabitat")$group
TukeyLetters_566$Microhabitat = rownames(TukeyLetters_566)

#########################################
#OTU.274
a <- ggplot(OTU.274, aes(x=Microhabitat, y=Reads, fill="OTU 274 Polymyxa betae")) + 
  geom_boxplot()+
  stat_boxplot(geom = "errorbar", width = 0.74, show.legend = F)+
  geom_text(data = TukeyLetters_274, aes(label = groups), y=75, size=5) +
  theme_classic()+
  theme(plot.margin = margin(1,1,1,1, "cm"))+
  theme(axis.text=element_text(size=14, face = "bold"), 
        axis.title=element_text(size=16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle=45, hjust = 1),
        legend.text = element_text(size = 15.5, face = "italic"), 
        legend.title = element_blank(), 
        #legend.position = "top",
        legend.direction = "horizontal",
        legend.position=c(0.4, 1))+
  scale_fill_manual(values = c("#014636"), 
                    limits = c("OTU 274 Polymyxa betae")) +
  scale_x_discrete(limits = c("Fresh Leaves","Bark","Detritus", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) 

a


#########################################
#OTU.230
b <- ggplot(OTU.230, aes(x=Microhabitat, y=Reads, fill="OTU 230 Spongospora nasturtii")) + 
  geom_boxplot()+
  stat_boxplot(geom = "errorbar", width = 0.74, show.legend = F)+
  geom_text(data = TukeyLetters_230, aes(label = groups), y=90, size=5) +
  theme_classic()+
  theme(plot.margin = margin(1,1,1,1, "cm"))+
  theme(axis.text=element_text(size=14, face = "bold"), 
        axis.title=element_text(size=16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle=45, hjust = 1),
        legend.text = element_text(size = 15.5, face = "italic"), 
        legend.title = element_blank(), 
        #legend.position = "top",
        legend.direction = "horizontal",
        legend.position=c(0.4, 1))+
  scale_fill_manual(values = c("#660033"), 
                    limits = c("OTU 230 Spongospora nasturtii")) +
  scale_x_discrete(limits = c("Fresh Leaves","Bark","Detritus", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) 

b
 

#########################################


#OTU.566
c <- ggplot(OTU.566, aes(x=Microhabitat, y=Reads, fill="OTU 566 unassigned Leptophryidae")) + 
  geom_boxplot()+
  expand_limits(y=c(NA, 33))+
  stat_boxplot(geom = "errorbar", width = 0.74, show.legend = F)+
  geom_text(data = TukeyLetters_566, aes(label = groups), y=31.5, size=5) +
  theme_classic()+
  theme(plot.margin = margin(1,1,1,1, "cm"))+
  theme(axis.text=element_text(size=14, face = "bold"), 
        axis.title=element_text(size=16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle=45, hjust = 1),
        legend.text = element_text(size = 15.5), 
        legend.title = element_blank(), 
        #legend.position = "top",
        legend.direction = "horizontal",
        legend.position=c(0.4, 1))+
  scale_fill_manual(values = c("#c4782b"), 
                    limits = c("OTU 566 unassigned Leptophryidae")) +
  scale_x_discrete(limits = c("Fresh Leaves","Bark","Detritus", "Deadwood","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil")) 

c 


#
library(ggpubr)

combi = ggarrange( a,b,c,
                   labels = c("A", "B","C"), 
                   ncol = 3, nrow = 1, font.label = list(size = 16, color = "black"),
                   common.legend = F, legend = "bottom", align = "h")+
  theme(plot.margin = margin(1,2,1,1, "cm")) 

combi

######## save plot


ggsave("Boxplot_300_336x168.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 336, height = 168, 
       units = "mm")
ggsave("Boxplot_600_336x168.pdf", plot = combi, 
       device = "pdf", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Boxplot_600_336x168.tif", plot = combi, 
       device = "tiff", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Boxplot_600_336x168.png", plot = combi, 
       device = "png", dpi = 600, width = 336, height = 168, 
       units = "mm")

