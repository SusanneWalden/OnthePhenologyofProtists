
#################################################################################################
################################# Differential Abbundances Spring vs Autumn #####################

rm(list = ls())


library(DESeq2)# for this you need the package "XML" which ist not availabla at cran, type:install.packages("XML", type = "binary")
library(apeglm)
library(magrittr)
library(ggplot2)
library(ggpubr)
#library(egg)
library(gridExtra)# for theme black
library(plyr)
library(scatterpie)

setwd("C:/Users/Susi/Desktop/season_server/Diff_Abb/")


#OTU_Table = as.data.frame(read.csv2("AllSeason_Cerc_En_Final_OTU_R_7633.csv",header = T, stringsAsFactors=TRUE))
OTU_Table = as.data.frame(read.csv2("AllSeason_Cerc_En_Final_OTU_R_7633_Tax_DiffAbb_OTU.csv",header = T, stringsAsFactors=TRUE))#season

species = OTU_Table[,8:ncol(OTU_Table)]
species_mat = as.matrix(species)
SampleMetadata = as.matrix(OTU_Table[,1:7])

#SampleMetadata=as.factor(SampleMetadata)
#####################

#OTU table into a DESeq2 object
ddsSeason = DESeqDataSetFromMatrix(countData = t(species_mat),
                                    colData = SampleMetadata,
                                    design = ~ Season)

#SampleMetadata=as.factor(SampleMetadata)

#Do the differential abundance analysis
ddsSeason = DESeq(ddsSeason, 
                   minReplicatesForReplace = Inf, #This means that Outliers should not be excluded, we want the maximum amount of information
                   sfType="poscounts",
                   quiet = T)

#Convert the analysis into a clearly readable results table
resSeason = results(ddsSeason, pAdjustMethod = "BH")

#Shrink the LFC (Log Fold Change)
#This shrinks the lfc into  a normal distribution, which is very handy for visualisations
resSeasonLFC = lfcShrink(ddsSeason, 
                          type = "apeglm", 
                          res = resSeason, 
                          coef = "Season_Spring_vs_Autumn", 
                          quiet = T)


######################
DESeq_Season_dependend_Subset = subset(resSeasonLFC, 
                                        resSeasonLFC$padj < 0.01 #highly significant pvalues less than 0.01
                                        & abs(resSeasonLFC$log2FoldChange) >= 1) #OTUs highly differential abundant, with an absolute LFC greater than 1

#Convert into a data_bframe and sort by LFC
data_b = data.frame(DESeq_Season_dependend_Subset[order(DESeq_Season_dependend_Subset$log2FoldChange),]$log2FoldChange)
colnames(data_b) = "log2FoldChange"

#For the visualisation, the values are inverted so ground-associated OTUs appear at the bottom and "Canopy-OTUs" at the top
data_b$log2FoldChange = data_b$log2FoldChange * -1



############

#Create empty matrices, one for the Ground Season and one for the Canopy
SpringMatrix = matrix(numeric(0), 0,0)
AutumnMatrix = matrix(numeric(0), 0,0)

#Then we fill the matrices with the number of reads per OTU and Season
for (a in as.factor(rownames(data_b))){ #The rownames are the names of the OTU
  a = get("a")
  get("OTU_Table")
  SpringMatrix[[a]] = sum(OTU_Table[[a]][OTU_Table$Season == "Spring"]) #All reads associated with Bark, Deadwood etc are summed
  AutumnMatrix[[a]] = sum(OTU_Table[[a]][OTU_Table$Season == "Autumn"]) #All reads associated with Soil and Leaf Litter are summed
}

#Convert Matrices to data_bframes and add column names
AutumnMatrix = as.data.frame(AutumnMatrix, header = F)
SpringMatrix = as.data.frame(SpringMatrix, header = F)
colnames(AutumnMatrix) = "Autumn"
colnames(SpringMatrix) = "Spring"

#Add the new columns to our data_b
data_b = cbind(data_b, AutumnMatrix, SpringMatrix)

#We don't need the baseMean or pvalues etc., so we subset the data_b and keep only the LFC (which will be the y-coordinates in the plot) and the Autumn and Spring Columns
data_b = subset(data_b, select = c("log2FoldChange", "Autumn", "Spring"))

#We need to add an extra column for the x-coordinates in the plot, this is just a sequence from 1 to the number of OTUs
data_b$xPosition = seq(1:nrow(data_b))


########## PLOT


P = data_b %>% #pipe the data_b into ggplot
  arrange(data_b$log2FoldChange, decreasing = T) %>% #Not sure if this is actually needed, because we already sorted the data_b. I keep it nevertheless
  mutate(x = factor(rownames(data_b), rownames(data_b))) %>% #Define the OTUs as x-coordinates
  ggplot(aes(x = x, y = data_b$log2FoldChange)) +  #x are the OTUs, y the LFC values
  geom_segment(aes(x=x, xend=x, #geom_segment is similar to a barplot, but you can control the width and other aesthetics
                   y=0, yend = data_b$log2FoldChange), 
               size=1.5, alpha=0.9, #size is the width of the bars, alpha the transparency
               #The color is a bit tricky. I want to plot smaller LFCs in lighter green/brown and higher LFC in darker colors
               #So I generate a sequential palette first from dark green to light green
               #Then I specify how many colors I need, this is the total number of OTUs with a positive LFC (sum data_b$lfc > 0)
               color = c(colorRampPalette(c("660033","#99004c"))(sum(data_b$log2FoldChange > 0)), #this one was purple
                         #For the Autumn OTUs, I generate a palette from light to dark brown as the first Autumn OTU has the "lowest" LFC
                         #Then I again specify how many brown colors I want, namely the number of OTUs with a negative LFC
                         colorRampPalette(c("#02818a","#014636"))(sum(data_b$log2FoldChange < 0))), 
               position = position_dodge(width = 0.5)) + #Can't really remember what this does
  theme_minimal(base_size = 5) + #The theme for the plot, less is more ;)
  theme(legend.position = "none", #We don't need a legend
        panel.border = element_blank(), #No panel borders
        axis.text.y = element_text(size = 14, vjust = 0.5, hjust = 0.5), #make the axis tick labels (OTU names) bigger and rotate them 90 degrees
        axis.text.x = element_text(size =12, angle = 90, vjust = 0.5, hjust = 0.5), #make the axis tick labels (OTU names) bigger and rotate them 90 degrees
        axis.title = element_text(size=16, face = "bold", hjust = 0.5)) + #Title should be in the middle (hjust 0.5) and bold (actual text not yet specified here)
  xlab("OTU") + #Axis title for x axis
  ylab("log2FoldChange") + #Axis title for y axis
  #labs(title = "Differentially Abundant OTUs") + #labs targets all labels, here I specify the text of the title
  #With geom_label we can specify our own labels, in this case "Spring" and "Autumn"
  theme(axis.title.y = element_text(size = rel(1), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1), angle = 00))+
  geom_label(x = 20, y = 6, label = "Autumn", label.padding = unit(0.8, "lines"), fill = "#660033", color = "white", size = 10) +
  geom_label(x = 65, y = -6, label = "Spring", label.padding = unit(0.8, "lines"), fill = "#014636", color = "white", size = 10)+ 
#Scatterpie is a mixture of piecharts and scatterplots
#I want to show the number of reads per Season per OTU in a piechart to identify OTUs which are differentialy abundant but not exclusively present in one Season
geom_scatterpie(data = data_b, #read the data_b
                #IMPORTANT: very contraintuitive, but dont put the columns in "quotationmarks" and also dont call the columns with e.g. data_b$xPosition, this will give errors
                aes(x = xPosition, y = log2FoldChange, r = 0.8), #define x and y coordinates and radius for each pie
                cols = c("Autumn", "Spring"), #the columns which will be used for the pie charts, here the number of reads per Ground and Canopy
                color = "darkgrey") + #color will only color the borders of the pies, the actual fill color is specified below
  scale_fill_manual(values=c("#660033", "#014636"))+  #fill the piecharts with green and brown
  #geom_point(size = 1.4, fill = "white", color = "white") #Here I plot white circles to make the piecharts look like donuts, not necessary, but I think it looks better
  geom_hline(yintercept = 0, size = 0.5, linetype=1)+ 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

P

#PF=P + coord_flip()

#PF


############### save & export

ggsave("Differential_Abbundance_OTU_300_336x168.jpeg", plot = P, 
       device = "jpeg", dpi = 300, width = 336, height =168, 
       units = "mm")
ggsave("Differential_Abbundance_OTU_600_336x168.pdf", plot = P,  
       device = "pdf", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Differential_Abbundance_OTU_600_336x168.tif", plot = P,  
       device = "tiff", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Differential_Abbundance_OTU_600_336x168.png", plot = P,  
       device = "png", dpi = 600, width = 336, height = 168, 
       units = "mm")
#################


####plot together NMDS + diff abb. (if you want to plot diffenrent grahics from different scripts together make sure to use specific variable/table names, otherwise you´ll overwrite them

#both = ggarrange(P,g, nrow=2,labels = c("A", "B"))
#both

#ggsave("Diff_NMDS_combined_5_diffcol.png", plot = both, 
#       device = "png", dpi = 300, width = 25, height = 25, #frontiers = length max 1 site, 18 mm width, 300 dpi
#       units = "cm")

#ggexport(L,P, filename = "test.pdf",
#         nrow = 2, ncol = 1)




############### wer kommt wo vor?
library(vegan)

specnumber(data_b$Spring)
specnumber(data_b$Autumn)

#wer ist diff abundant?
#Autumn
specnumber(data_b$log2FoldChange > 0)
#Spring
specnumber(data_b$log2FoldChange < 0)


#write.csv(data_b, file = "Differential_Abbundaces_LFC__Spring_Fall_withposcounts.csv")



