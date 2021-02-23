rm(list = ls())
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(funrar)






#Leizig cerco + Endomyxa 7633
setwd("C:/Users/Susi/Desktop/season_server/Functional_Traits/F19/")
OTU_Table = as.data.frame(read.csv2("AllSeason_F19_Cerc_En_Final_OTU_R_7633.csv",header = T))


SampleMetadata = OTU_Table[,1:5]
rownames(OTU_Table) = SampleMetadata$SampleID
rownames(SampleMetadata) = SampleMetadata$SampleID
species = OTU_Table[,6:ncol(OTU_Table)]
species_mat = as.matrix(species)
species_mat = make_relative(species_mat)

Microhabitat = SampleMetadata$Microhabitat
TreeSpecies = SampleMetadata$TreeSpecies




#Next we normalise our species matrix with a hellinger transformation and perform 
#the redundancy analysis. We want to test both the microhabitats as well as the tree species as an explanatory variable:



species_mat_he = decostand(species_mat, "hellinger")
species_mat_rda = rda(species_mat_he~Microhabitat+TreeSpecies)
summary_species_mat_rda = summary(species_mat_rda, scaling = 2)

# For the microhabitats (rows 1-9)
df1  <- data.frame(summary_species_mat_rda$centroids[1:9,1:2])
rownames(df1) = gsub("Microhabitat", "", rownames(df1))
df1$Microhabitat = rownames(df1)

# For the species (OTUs)
df2  <- data.frame(summary_species_mat_rda$species[,1:2])

# For the tree species (rows 10-12)
df3 = data.frame(summary_species_mat_rda$centroids[10:12,1:2])
rownames(df3) = gsub("TreeSpecies", "", rownames(df3))
df3$TreeSpecies = rownames(df3)




ColorVectorTreeSpecies = vector()
ColorVector2 = vector()

# Again, these colors correspond to the YlGn Palette
for (a in rownames(df1)){
  if(a == "Arboreal Soil"){ColorVector2 = c(ColorVector2, "#FFE4B5")}
  else if (a == "Bark"){ColorVector2 = c(ColorVector2, "#FF944D")}
  else if (a == "Deadwood"){ColorVector2 = c(ColorVector2, "#BF6692")}
  else if (a == "Fresh Leaves"){ColorVector2 = c(ColorVector2, "#138F6A")}
  else if (a == "Hypnum"){ColorVector2 = c(ColorVector2, "#A8ACFF")}
  else if (a == "Lichen"){ColorVector2 = c(ColorVector2, "#8B8989")}
  else if (a == "Orthotrichum"){ColorVector2 = c(ColorVector2, "#4363d8")}
  else if (a == "Leaf Litter"){ColorVector2 = c(ColorVector2, "#9A6324")}
  else if (a == "Soil"){ColorVector2 = c(ColorVector2, "#800000")}
}


for (a in rownames(df3)){
  if(a == "T.cordata"){ColorVectorTreeSpecies = c(ColorVectorTreeSpecies, "#6b4e71")}
  else if (a == "Q.robur"){ColorVectorTreeSpecies = c(ColorVectorTreeSpecies, "#53687e")}
  else if (a == "F.excelsior"){ColorVectorTreeSpecies = c(ColorVectorTreeSpecies, "#c2b2b4")}
}

g = ggplot(df2, aes(x=RDA1, y=RDA2)) + 
  geom_point(aes(fill = "OTUs"), size=3, color = "darkslategrey") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  geom_segment(data=df1, aes(x=0, xend=RDA1, y=0, yend=RDA2, 
                             linetype = "Microhabitat", color = Microhabitat), 
               arrow=arrow(length=unit(0.01,"npc")), 
               size = 1.5, alpha = 0.9) +
  scale_color_manual(values = unique(ColorVector2), 
                     breaks = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) +
  geom_label(data=df1, 
             aes(x=RDA1,y=RDA2,label=rownames(df1),
                 #hjust=c(0,0,1,1,0,0,1,1,0), 
                 hjust=0.5*(1-sign(RDA1)),
                 vjust=0.5*(1-sign(RDA2))), 
             fill=ColorVector2, color = "black", size=4) + 
  geom_segment(data=df3, aes(x=0, xend=RDA1, y=0, yend=RDA2, 
                             linetype = "Tree Species", color = TreeSpecies), 
               arrow=arrow(length=unit(0.01,"npc")), 
               size = 1.5, alpha = 0.75, 
               color = ColorVectorTreeSpecies) +
  geom_label(data=df3, 
             aes(x=RDA1,y=RDA2,label=rownames(df3),
                 hjust=0.5*(1-sign(RDA1)),
                 vjust=0.5*(1-sign(RDA2))), 
             fill=ColorVectorTreeSpecies, color = "black", size=4) + 
  scale_x_continuous(limits = c(-0.8, 0.8)) + # wenn daten außerhalb der skala liegen können sie nicht mehr abgebildet werden
  scale_y_continuous(limits = c(-0.8, 0.8)) +
  labs(title = "Spring 2019", 
       x = paste0("RDA1\nVariance explained: ",
                  round(as.numeric(summary_species_mat_rda$cont$importance["Proportion Explained", "RDA1"]), digits = 4) * 100, 
                  "%"), 
       y = paste0("RDA2\nVariance explained: ",
                  round(as.numeric(summary_species_mat_rda$cont$importance["Proportion Explained", "RDA2"]), digits = 4) * 100, 
                  "%")) +
  theme_minimal() + 
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16, face = "bold"), 
        legend.position = "right",
        legend.direction = "vertical", 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5), 
        legend.title.align = 0.5,
        axis.title = element_text(size =16, face = "bold"), 
        axis.text = element_text(size = 14)) +
  scale_linetype_manual(name = "Explanatory Variable", 
                        values = c("Microhabitat"=1, "Tree Species"=2)) + 
  scale_fill_manual(name = "Response Variable", 
                    values = "darkslategrey") +
  guides(color = F, 
         linetype=guide_legend(keywidth = 3, keyheight = 1, order = 1), 
         fill = guide_legend(order = 2))
g
A=g
B=g
C=g
D=g


# plot
combi = ggarrange( A,B,C,D,
                   labels = c("A", "B","C","D"), 
                   ncol = 2, nrow = 2, font.label = list(size = 16, color = "black"), 
                   common.legend = T, legend = "right", 
                   align = "v")+
  theme(plot.margin = margin(1,1,1,1, "cm")) 
combi


##############################
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