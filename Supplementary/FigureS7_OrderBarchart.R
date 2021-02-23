

################################################### ORDER BARCHART#########################################


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



######################################## Order #################################################


colnames(TAX) = c("OTU_Number", "Order", "Family", "Genus", "Species", "ReferenceID", "PercentID", "Nutrition", "Morphology", "Abundance")
TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)
TAX$Order = as.factor(TAX$Order)

#```{r Aggregate Table}
OrderTable = OTU_Table
colnames(OrderTable) = TAX$Order

# First aggregate by Habitat
HabitatAggregatedOrderTable = aggregate(OrderTable, by = list(SampleMetadata$Microhabitat), FUN = sum)
rownames(HabitatAggregatedOrderTable) = HabitatAggregatedOrderTable$Group.1
HabitatAggregatedOrderTable = HabitatAggregatedOrderTable[,-1]


# if you want to check the number of OTUs instead of the amount of reads:
#HabitatAggregatedOrderTable[,] = ifelse(HabitatAggregatedOrderTable[,] > 0, 1, 0)

#Then aggregate the (transposed) table by Order
HabitatAggregatedOrderTable = aggregate(t(HabitatAggregatedOrderTable), 
                                        by = list(TAX$Order), 
                                        FUN = sum)
rownames(HabitatAggregatedOrderTable) = HabitatAggregatedOrderTable$Group.1
HabitatAggregatedOrderTable = HabitatAggregatedOrderTable[,-1]
HabitatAggregatedOrderTable = as.data.frame(HabitatAggregatedOrderTable)

# melt the table so ggplot can interpret it as a histogram
data = HabitatAggregatedOrderTable
data$Order = rownames(data)
data2 = melt(data, id.vars = "Order")

# add the Stratum
data3 = data2
data3$Stratum = ifelse(data2$variable %in% c("Leaf Litter", "Soil"), "Ground", "Canopy")

# for the plotting we need to change the column names
colnames(data3) = c("score", "item", "value", "family")

###############plot as  barchat


#data3$item <- factor(data3$item, levels=c("Fresh Leaves","Deadwood","Bark", "Arboreal Soil","Hypnum","Orthotrichum","Lichen","Leaf Litter", "Soil"))
data3$item <- factor(data3$item, levels=c("Soil","Leaf Litter","Lichen","Orthotrichum","Hypnum", "Arboreal Soil","Bark","Deadwood","Fresh Leaves"))



####### plot
o=ggplot(data3, aes(fill=score, y=value, x=item)) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("#008B8B","#ADD8E6","#EC6B25","#ffd27f","#228B22","#90EE90","#6f0000","#f6beb2","#730073","#cc99cc","#4c3706",
                               "#998b6a","#b25162","#ffabba","#656565","#c2c2c2","#c7b823","#e9df81","#2ca898","#9be8df","#73315d","#e0a9cd",
                               "#ffff4c","#ffffb2","#0a2d22"), 
                    limits = c("Cercomonadida","Cryomonadida","Cryptofilida", "Desmothoracida","Ebriida","Euglyphida",
                               "Glissomonadida","Limnofilida","Marimonadida","Novel-clade-10","Novel-clade-12","Novel-clade-9",
                               "Pansomonadida","Plasmodiophorida","Spongomonadida","Tectofilosida","Thaumatomonadida",
                               "Tremulida","unassigned Cercozoa","unassigned Filosa","unassigned Granofilosea",
                               "unassigned Imbricatea","unassigned Thecofilosea","Vampyrellida","Ventricleftida"),
                    guide = guide_legend(ncol = 1))+
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  labs(fill = "Order", y= "Relative Abundances", x = "Microhabitat",title = "Spring 2019")+
  theme_minimal()+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16, face = "bold"),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.text=element_text(size=14), 
        axis.title=element_text(size=16,face="bold"))+
  geom_label(x = 6.5, y = 0.75, fill = "darkolivegreen4", color = "white", 
             label = "Canopy", size = 8 ) +
  geom_label(x = 1.5, y = 0.75, fill = "burlywood4", color = "white", 
             label = "Ground", size = 8)+
  geom_vline(xintercept = 2.50, size = 0.6, linetype=1.5, color = "#6A785F")+
  coord_flip()
o

#geom_rect(data=NULL,aes(xmin=0,xmax=2.5,ymin=-Inf,ymax=Inf),
#fill="lightgrey")+
#geom_rect(data=NULL,aes(xmin=2.5,xmax=9.5,ymin=-Inf,ymax=Inf),
#fill="lightgreen", alpha=0.9)  
#A= o + annotate("rect", xmin=2.5, xmax=9.5, ymin=0, ymax=Inf, alpha=0, fill = "darkolivegreen4",colour="darkolivegreen4")

#A = A + annotate("rect", xmin=0, xmax=2.5, ymin=0, ymax=Inf, alpha=0, fill="burlywood4",colour="burlywood4")

A = o
B = o
C = o
D = o

###################### print data for percentages

#write.csv(data, file = "Spring19_Order_Values.csv")


####plot together


combi = ggarrange( A,B,C,D,
                   labels = c("A", "B","C","D"), 
                   ncol = 2, nrow = 2, font.label = list(size = 16, color = "black"), 
                   common.legend = T, legend = "right", 
                   align = "v")+
  theme(plot.margin = margin(1,1,1,1, "cm")) 
combi


##############################
ggsave("OrderCombined_300_336x252.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 336, height = 252, 
       units = "mm")
ggsave("OrderCombined_600_336x252.pdf", plot = combi, 
       device = "pdf", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("OrderCombined_600_336x252.tif", plot = combi, 
       device = "tiff", dpi = 600, width = 336, height = 252, 
       units = "mm")
ggsave("OrderCombined_600_336x252.png", plot = combi, 
       device = "png", dpi = 600, width = 336, height = 252, 
       units = "mm")


