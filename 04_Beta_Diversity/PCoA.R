rm(list = ls())

library(vegan)
library(ggplot2)
library(plyr)
library(funrar)# for make relative function



######################################## All microhabitats(290) 7633  ######################

rm(list = ls()) 
setwd("C:/Users/Susi/Desktop/season_server/")
OTU_Table = as.data.frame(read.csv2("AllSeason_Cerc_En_Final_OTU_R_7633.csv",header = T))


species = OTU_Table[,8:ncol(OTU_Table)]
species_mat = as.matrix(species)
species_mat = make_relative(species_mat)
#species_mat = log(species_mat +1)
SampleMetadata = as.data.frame(OTU_Table[,1:7])

dim(species_mat)#  290 samples 783   OTUs as column, sites as rows
OTUCer<-species_mat


############################################### Only leaves ###################################

#rm(list = ls())
setwd("C:/Users/Susi/Desktop/season_server/PCA/")

DatabaseCer<-read.table(file="AllSeason_OTU_PCA_Leaves_7633.txt", header=T)
dim(DatabaseCer)   # 783 25
# Select the OTUs and leave the taxonomy and function - the first
#591 columns and transpose

OTUCer<-t(DatabaseCer[ ,c(1:25)])
OTUCer_mat = as.matrix(OTUCer)
OTUCer_mat = make_relative(OTUCer_mat)
#OTUCer_mat = log(OTUCer_mat +1)

dim(OTUCer)    #  25 783   OTUs as column, sites as rows

##################################################################################


###functions to visualize groups of samples in ordination with the function chull:
  # Computes the subset of points which lie on the convex hull of the
  #set of points specified, using a two-column matrix
find_hull.nmds <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
find_hull.pcoa <- function(df) df[chull(df$V1, df$V2), ]

## normalize by sample
OTUCer.tss <- decostand(OTUCer, "total")
OTUCer.tss <- OTUCer.tss*100 ### all samples combine to 100 observations
rowSums(OTUCer.tss)  # OTUs as columns and sites as rows, 591   2101

# Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis (Gower, 1966).
# Needs as input a distance matrix (full), a k value (max number of dimensions, default=2)
# to know the variance explained by each axe, add eig = TRUE the output will contain item GOF (= "goodness of fit")
# which is as close to the R2 as you can get in PCoA.
pcoa <- cmdscale(vegdist(OTUCer.tss, "bray"), eig=TRUE, k=2) 
#pcoa <- cmdscale(vegdist(OTUCer.tss, "bray"), eig=TRUE, k=2)   
#$GOF 0.4303537 0.4616277 # here you have to add k = maximum number of dimensions to represent the data
# explained variation by axis 1 and axis 2 of PCoA
explainedvar1<-round(pcoa$eig[1]/sum(pcoa$eig),2)*100
explainedvar2<-round(pcoa$eig[2]/sum(pcoa$eig),2)*100

explainedvar1 #  all microhabitats: the first x-axis explains 25% of variation, the second y-axis 11% of variation
explainedvar2 # only leaves: the first x-axis explains 21% of variation, the second y-axis 17% of variation 


sum_eig<-sum(explainedvar1,explainedvar2)
sum_eig # all microhabitats = 36%, Leaves 38% = half of the variation is explained by the first two axes

# There are two criteria to assess whether the first few PCoA axes capture a disproportionally large amount of variation:
# Kaise-Guttman criterion: the first few axes should be larger than the average of all eigenvalues
# Broken Stick Model: compare first few axes to expectations of Broken Stick Model

# plot eigenvalues
plot(pcoa$eig, xlab = "PCoA", ylab="Eigenvalue", las=1, cex.lab =
       1.5, pch=16, xlim =c(0,20))
abline (h= mean(pcoa$eig), lty=2, col= "blue")
b_stick<-bstick(8,sum(pcoa$eig))
lines(1:8, b_stick, lty=4, lwd=2, col= "red")
legend("topright", legend=c("Avg Eigenvalue", "Broken Stick"),
       lty=c(2,4), bty="n", col=c("blue", "red"))

# Result: eigenvalues of the first 7 axes are larger than the average of all eigenvalues but smaller than the expectation of the Broken Stick Model
# still ok i guess!
############################################################

names(pcoa)  # "points" = Sample_name, then (1),(2)=coordinates "eig"=vector   "x" NULL??     "ac"??  "GOF"  OK
pcoa<-as.data.frame(cbind(pcoa$points, pcoa$eig))   # keep only the first 3 columns as data frame
dim(pcoa)  # 290 oder 25  3

# create 4 dataframes for each season 
####### all microhabitats all samples
#pcoaH17<-as.data.frame(pcoa[c(244:324),])
#pcoaF18<-as.data.frame(pcoa[c(1:81),])
#pcoaH18<-as.data.frame(pcoa[c(163:243),])
#pcoaF19<-as.data.frame(pcoa[c(82:162),])

####### all microhabitats 7633 only 290 samples
pcoaH17<-as.data.frame(pcoa[c(211:290),])
pcoaF18<-as.data.frame(pcoa[c(1:65),])
pcoaH18<-as.data.frame(pcoa[c(142:210),])
pcoaF19<-as.data.frame(pcoa[c(66:141),])

###leaves 
#pcoaH17<-as.data.frame(pcoa[c(28:36),])
#pcoaF18<-as.data.frame(pcoa[c(1:9),])
#pcoaH18<-as.data.frame(pcoa[c(19:27),])
#pcoaF19<-as.data.frame(pcoa[c(10:18),])

###leaves 7633
pcoaH17<-as.data.frame(pcoa[c(11:19),])
pcoaF18<-as.data.frame(pcoa[c(1:5),])
pcoaH18<-as.data.frame(pcoa[c(20:25),])
pcoaF19<-as.data.frame(pcoa[c(6:10),])

# create a column of names called "season" so we can select them
#later and add it to the dataframe

# all microhabitats
pcoaH17$Season=rep("Autumn 2017",nrow(pcoaH17))
pcoaF18$Season=rep("Spring 2018",nrow(pcoaF18))
pcoaH18$Season=rep("Autumn 2018",nrow(pcoaH18))
pcoaF19$Season=rep("Spring 2019",nrow(pcoaF19))


pcoabyseason<-rbind(pcoaH17, pcoaF18, pcoaH18, pcoaF19)

# calculate the points for the PCoA by applying the function to visualize data
hulls.pcoaSe <- ddply(pcoabyseason,"Season", find_hull.pcoa)
hulls.pcoaSe
#V1           V2            V3      season
#1   0.26133702 -0.085069505 -9.983158e-02 Autumn 2017
#2   0.15555069 -0.182483670 -7.917436e-02 Autumn 2017
#3   0.12067476 -0.214267439 -6.784525e-02 Autumn 2017



# graph by season
#PCOA_Se <- ggplot()+
#  geom_point(data=pcoabyseason, aes(x=V1,y=V2, col=Season,
#                                    shape=Season), size=4)+
#  scale_shape_manual(values=c(15,15,16,16,17,17))+
#  
#  scale_color_manual(values=c("cadetblue","#2F4F4F","yellowgreen","olivedrab"))+
#  xlab('PCoA1')+
#  ylab('PCoA2')+
#  ggtitle("All Microhabitats - PCoA by Season")
#PCOA_Se

############ Final Colors

# graph by season
PCOA_Se <- ggplot()+
  geom_point(data=pcoabyseason, aes(x=V1,y=V2, col=Season,
                                    shape=Season), size=3)+
  geom_polygon(data=hulls.pcoaSe,aes(x=V1, y=V2, group = Season, fill = Season),alpha=0.2)+
  scale_shape_manual(values=c(16,16,17,17))+
  theme_minimal()+
  scale_fill_manual(values = c("#094d35","#660033","#B5A642","#1c9099"), 
                    limits = c("Autumn 2017","Autumn 2018","Spring 2018","Spring 2019")) +  
  scale_color_manual(values=c("#094d35","#660033","#B5A642","#1c9099"))+
  xlab('PCoA1')+
  ylab('PCoA2')+
  #ggtitle("Only Leaves - PCoA by Season")+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=16, face = "bold"), 
        legend.text = element_text(size = 16), 
        #legend.title = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),#delete legend title
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))+
  labs(y= "PCoA2 - 17%", x = "PCoA1 - 21%")#change axis titles
        
PCOA_Se

A<-PCOA_Se

C<-PCOA_Se # Only Leaves

###################### plot both together



combi = ggarrange( A,C,
                   labels = c("A", "B"), 
                   ncol = 2, nrow = 1, font.label = list(size = 16, color = "black"),
                   common.legend = T, legend = "bottom", align = "h")+
          theme(plot.margin = margin(1,1,1,1, "cm")) 
 
combi



###save plot

ggsave("PCA_300_336x168.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 336, height = 168, 
       units = "mm")
ggsave("PCA_600_336x168.pdf", plot = combi, 
       device = "pdf", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("PCA_600_336x168.tif", plot = combi, 
       device = "tiff", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("PCA_600_336x168.png", plot = combi, 
       device = "png", dpi = 600, width = 336, height = 168, 
       units = "mm")

