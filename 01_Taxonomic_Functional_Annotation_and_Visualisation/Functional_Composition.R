

########################################################################################

########################################## Nutrition ####################################

#########################################################################################


rm(list= ls())


library(ggplot2)
library(reshape2)
library(phenotypicForest)
library(magrittr)


setwd("C:/Users/Susi/Desktop/season_server/Functional_Traits/F19/")

TAX =read.csv2("AllSeason_Cerc_En_Taxonomy_R_nicerNames.csv", header = T)    #same taxonom< file for seasonal comparison, data had been clustered together 
OTU_Table = as.data.frame(read.csv2("AllSeason_F19_Cerc_En_Final_OTU_R_7633.csv",header = T))


SampleMetadata = OTU_Table[,1:5]
OTU_Table = OTU_Table[,6:ncol(OTU_Table)]

# Calculate the total abundance of each OTU by summing the columns
Abundances = colSums(OTU_Table)
TAX = cbind(TAX, Abundances)


######################################################### Nutrition ##################################################

colnames(TAX) = c("OTU_Number", "Order", "Family", "Genus", "Species", "ReferenceID", "PercentID", "Nutrition", "Morphology", "Abundance")
TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)
TAX$Nutrition = as.factor(TAX$Nutrition)





## Aggregate Table
#Now we need to aggregate the Table. We combine the forms of Nutrition 
#and the Microhabitats, so we have a 9x4 dataframe containing the number of reads per Nutrition per Habitat

#```{r Aggregate Table}
NutritionTable = OTU_Table
colnames(NutritionTable) = TAX$Nutrition

# First aggregate by Habitat
HabitatAggregatedNutritionTable = aggregate(NutritionTable, by = list(SampleMetadata$Microhabitat), FUN = sum)
rownames(HabitatAggregatedNutritionTable) = HabitatAggregatedNutritionTable$Group.1
HabitatAggregatedNutritionTable = HabitatAggregatedNutritionTable[,-1]


# if you want to check the number of OTUs instead of the amount of reads:
# this is the point to convert it into a presence/absence matrix like this:
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

#put in the right Nutrition so that it has the right Nutrition as the legend
#data3$item <- factor(data3$item, levels=c("Fresh Leaves","Deadwood","Bark", "Arboreal Soil","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil"))
data3$item <- factor(data3$item, levels=c("Soil","Leaf Litter","Lichen","Orthotrichum","Hypnum", "Arboreal Soil","Bark","Deadwood","Fresh Leaves"))
#data3$score <- factor(data3$score, levels=c("bacterivore", "eukaryvore","omnivore","autotroph","plant parasite","not plant parasite", "undetermined"))


#plot stacked barplot with nutrition

o=ggplot(data3, aes(fill=score, y=value, x=item)) + 
  scale_fill_manual(values = c("#4E9258","#7D053D", "#7E587E", 
                               "coral2","#A9A9A9","#348781","#d3a540"), 
                    limits = c("autotroph","bacterivore", "eukaryvore", 
                               "not plant parasite","omnivore","plant parasite", "undetermined")) +
  geom_bar(position="fill", stat="identity") + 
  #facet_grid(cols = vars(Stratum), scales = "free", space = "free") 
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  labs(fill = "Nutrition", y= "Relative Abundances", x = "Microhabitat",title = "Spring 2019")+
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
  geom_vline(xintercept = 2.5, size = 0.6, linetype=1.5, color = "#6A785F")


#geom_rect(data=NULL,aes(xmin=0,xmax=2.5,ymin=-Inf,ymax=Inf),
#fill="lightgrey")+
#geom_rect(data=NULL,aes(xmin=2.5,xmax=9.5,ymin=-Inf,ymax=Inf),
#fill="lightgreen", alpha=0.9)  
#A= o + annotate("rect", xmin=2.5, xmax=9.5, ymin=0, ymax=Inf, alpha=0, fill = "darkolivegreen4",colour="darkolivegreen4")

#A = A + annotate("rect", xmin=0, xmax=2.5, ymin=0, ymax=Inf, alpha=0, fill="burlywood4",colour="burlywood4")
o

D = o + coord_flip()
D 


# you have to repeat it with every season before you can combine all graphs

###################### print data for percentages

write.csv(data, file = "Spring19_Nutrition_Values_cor.csv")


####plot together


combi = ggarrange( A,B,C,D,
                   labels = c("A", "B","C","D"), 
                   ncol = 2, nrow = 2, font.label = list(size = 16, color = "black"), 
                   common.legend = T, legend = "right", 
                   align = "v")+
          theme(plot.margin = margin(1,1,1,1, "cm")) 
combi


##############################
ggsave("NutritionCombined_300_336x252.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 336, height = 252, 
       units = "mm")
ggsave("NutritionCombined_600_336x252.pdf", plot = combi, 
       device = "pdf", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("NutritionCombined_600_336x252.tif", plot = combi, 
       device = "tiff", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("NutritionCombined_600_336x252.png", plot = combi, 
       device = "png", dpi = 600, width = 336, height = 252, 
       units = "mm")


