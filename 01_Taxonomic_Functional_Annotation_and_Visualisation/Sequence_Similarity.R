rm(list = ls())

library(ggplot2)
library(plyr)




# LEIPZIG CERC-EN-ALLSEASON


setwd("C:/Users/Susi/Desktop/season_server/")

TAX =read.csv2("AllSeason_Cerc_En_Taxonomy_R.csv", header = T)    #same taxonom< file for seasonal comparison, data had been clustered together 
OTU_Table = as.data.frame(read.csv2("AllSeason_Cerc_En_Final_OTU_R_7633.csv",header = T))


SampleMetadata = OTU_Table[,1:7]####aufpassen bei AllSeason_Cerc_En_Final_OTU_R existiert eine zusätzliche zeile für die jahreszeit
OTU_Table = OTU_Table[,8:ncol(OTU_Table)]
Abundances = colSums(OTU_Table)
TAX = cbind(TAX, Abundances)
TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)

######################################################################################

# This function rounds the value to the provided base
mround <- function(x,base){
  base*round(x/base)
} 


# Here I divided by 100 to obtain values between 0 and 1
TAX$PercentIDrounded = mround(TAX$PercentID, 1) / 100    #######  !!!! darauf achten, dass percent ID durch komma und nicht durch punkt getrennt sind

# Aggregate Abundances with the same percent identity
AggregatedTAXpercent = ddply(TAX, "PercentID", numcolwise(sum))#counts abundances of seq. with same percent ID
AggregatedTAXpercent$PercentID = AggregatedTAXpercent$PercentID / 100

# Expand the table. 
# This means that if e.g. 1000 sequences have the percent identity of 70.5 we print 70.5 1000 times

AggregatedTAXpercent_expanded = AggregatedTAXpercent[rep(row.names(AggregatedTAXpercent), AggregatedTAXpercent$Abundance), c(1,3)]



#####################################################################################################
###########################  READ SIMILARITIES ######################################################
#####################################################################################################


g = ggplot(AggregatedTAXpercent_expanded[AggregatedTAXpercent_expanded > 0, ], aes(x = PercentID)) +
  geom_histogram(aes(y = ..count..), 
                 alpha = 0.7, color = "darkgrey", fill="#094d35",
                 bins = 25)+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  labs(x = "Similarity to reference database", 
       y = "Number of reads") +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=16, face = "bold"))+ 
        #plot.title = element_text(size = 20, face = "bold", hjust = 0.7), 
        #legend.text = element_text(size = 12), 
        #legend.title = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.3)))

g


#####################################################################################################
###########################  OTU SIMILARITIES #######################################################
#####################################################################################################

# LEIPZIG CERC-EN-ALLSEASON
#rm(list = ls())

#setwd("C:/Users/Susi/Desktop/season_server/")

#TAX =read.csv2("AllSeason_Cerc_En_Taxonomy_R.csv", header = T)    #same taxonom< file for seasonal comparison, data had been clustered together 
#OTU_Table = as.data.frame(read.csv2("AllSeason_Cerc_En_Final_OTU_R_7633.csv",header = T))

#SampleMetadata = OTU_Table[,1:7]####aufpassen bei AllSeason_Cerc_En_Final_OTU_R existiert eine zusätzliche zeile für die jahreszeit
#OTU_Table = OTU_Table[,8:ncol(OTU_Table)]
#Abundances = colSums(OTU_Table)
#TAX = cbind(TAX, Abundances)
#TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)


################## subset 

L_TAX_subset = subset(TAX, select = "PercentID")
L_TAX_subset$Biome = "Temperate Forest"
L_TAX_subset$PercentID = L_TAX_subset$PercentID / 100




######## only Leipzig
g_OTU = ggplot(L_TAX_subset[L_TAX_subset$PercentID > 0, ], aes(x = PercentID)) +
  geom_histogram(aes(y = ..count..), 
                 alpha = 0.7, 
                 color = "darkgrey",fill="#094d35", 
                 bins = 25) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  labs(x = "Similarity to reference database", 
       y = "Number of OTUs") +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=16, face = "bold"))+ 
  #plot.title = element_text(size = 20, face = "bold", hjust = 0.7), 
  #legend.text = element_text(size = 12), 
  #legend.title = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.3)))


g_OTU 


###################### plot both together
library(ggpubr)



combi = ggarrange( g,g_OTU,
                   labels = c("A", "B"), 
                   ncol = 2, nrow = 1, font.label = list(size = 16, color = "black"),
                   common.legend = T, legend = "right", align = "h")+
        theme(plot.margin = margin(1,1,1,1, "cm")) 

combi

# save plot 
 
ggsave("Similarity_300_336x168.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 336, height = 168, 
       units = "mm")
ggsave("Similarity_600_336x168.pdf", plot = combi, 
       device = "pdf", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Similarity_600_336x168.tif", plot = combi, 
       device = "tiff", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Similarity_600_336x168.png", plot = combi, 
       device = "png", dpi = 600, width = 336, height = 168, 
       units = "mm")





