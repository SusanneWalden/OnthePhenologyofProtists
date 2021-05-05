Variation Partitioning Analysis
================

In the following section we want to calculate the proportion of variation in species composition explained by each variable.

Load data
---------

``` r
rm(list = ls()) 

library(vegan)
library(plyr)
library(funrar) 
library(ggplot2)
library(ggforce)
library(ggpubr)


OTU_Table = as.data.frame(read.csv2("../00_Data/05_Cercozoa_Seasonal_OTU_Table_min-freq-7633_transposed_withMetadata.csv",header = T))

#All
species = OTU_Table[,18:ncol(OTU_Table)]
species.mat = as.matrix(species)
species.mat = make_relative(species.mat)
SampleMetadata = as.data.frame(OTU_Table[,1:17])

#We need the Hellinger-transformed OTUs table with OTUs as columns and a table of metadata: Microhabitat, season and treespecies.
OTUCer.hel = decostand(species.mat, "hellinger")


#Ground
OTU_Ground <- OTU_Table[ which(OTU_Table$Stratum=='Ground'), ]

species.Ground = OTU_Ground[,18:ncol(OTU_Ground)]
species.mat.Ground = as.matrix(species.Ground)
species.mat.Ground = make_relative(species.mat.Ground)
SampleMetadata.Ground = as.data.frame(OTU_Ground[,1:17])

OTUCer.hel.Ground = decostand(species.mat.Ground, "hellinger")

#Canopy
OTU_Canopy <- OTU_Table[ which(OTU_Table$Stratum=='Canopy'), ]

species.Canopy = OTU_Canopy[,18:ncol(OTU_Canopy)]
species.mat.Canopy = as.matrix(species.Canopy)
species.mat.Canopy = make_relative(species.mat.Canopy)
SampleMetadata.Canopy = as.data.frame(OTU_Canopy[,1:17])

OTUCer.hel.Canopy = decostand(species.mat.Canopy, "hellinger")
```

Calulate Variation Partitioning
-------------------------------

First we plot the variation partitioning of our variables with the `varprat` function within the *vegan* package. Variation partitioning (using RDA, CCA or db-RDA) among up to four matrices of environmental variables. First argument (Y) is dependent variable, usually the matrix of species composition (the function calculates RDA, or, if chisquare = TRUE, CCA), but could be also only a single variable (in that case it calculates linear regression) or distance matrix (applying db-RDA using the function capscale). Further arguments (up to four) are (groups of) explanatory variables. The function uses either formula interface (with ~, see examples) or matrices.

``` r
#All
Cervarp.all <- varpart(OTUCer.hel, ~ Season, ~ Microhabitat, ~TreeSpecies, data = SampleMetadata)
#Ground
Cervarp.Ground <- varpart(OTUCer.hel.Ground, ~ Season, ~ Microhabitat, ~TreeSpecies, data = SampleMetadata.Ground)
#Canopy
Cervarp.Canopy <- varpart(OTUCer.hel.Canopy, ~ Season, ~ Microhabitat, ~TreeSpecies, data = SampleMetadata.Canopy)
```

Plot Variation Partitioning
---------------------------

``` r
#All
plot(Cervarp.all, digits = 1, Xnames = c('Season', 'Microhabitat', 'TreeSpecies'), bg = c('navy', 'tomato', 'green'))
```

![](Variation_Partitioning_files/figure-markdown_github/Plot%20Varprat-1.png)

``` r
#Ground
plot(Cervarp.Ground, digits =1, Xnames = c('Season', 'Microhabitat', 'TreeSpecies'), bg = c('navy', 'tomato', 'green'))
```

![](Variation_Partitioning_files/figure-markdown_github/Plot%20Varprat-2.png)

``` r
#Canopy
plot(Cervarp.Canopy, digits = 1, Xnames = c('Season', 'Microhabitat', 'TreeSpecies'), bg = c('navy', 'tomato', 'green'))
```

![](Variation_Partitioning_files/figure-markdown_github/Plot%20Varprat-3.png)

Plot with ggplot
----------------

Unfortunately it is not so easy to modify the plot afterwards. Consequently, we repeat the plotting with a dummy dataframe in *ggplot* to modify coloration and save settings.

``` r
#create dummy dataframe
df.venn.all <- data.frame(x = c(3, 1, 2),y = c(1, 1,2.8),labels = c('Microhabitat\n0.32', 'Season\n0.01',"Tree species\n0.03"))
#Ground 
df.venn.Ground<- data.frame(x = c(3, 1, 2),y = c(1, 1,2.8),labels = c('Microhabitat\n0.31', 'Season\n0.05',"Tree species\n0.02"))
#Canopy
df.venn.Canopy<- data.frame(x = c(3, 1, 2),y = c(1, 1,2.8),labels = c('Microhabitat\n0.18', 'Season\n0.02',"Tree species\n0.05"))

#Plot
#All
Varpart <- ggplot(df.venn.all, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
                  geom_circle(alpha = .6, size = 1, colour = "#3f3a3a",show.legend = FALSE) +
                  coord_fixed()+
                  annotate("text", x = df.venn.all$x , y = df.venn.all$y,label=df.venn.all$labels ,size = 6) +
                  scale_fill_manual(values = c("#660033","#B5A642","#094d35")) +
                  theme_void()
                  #annotate("text", x = 2 , y =1,label="A and B" ,size = 4) +
                  #annotate("text", x = 1.35 , y =2,label="B and C" ,size = 4) +
                  #annotate("text", x = 2.7 , y =2,label="A and C" ,size = 4) +
                  #annotate("text", x = 2 , y =1.6,label="A and B and C" ,size = 2)+theme_void()

Varpart
```

![](Variation_Partitioning_files/figure-markdown_github/Ggplot%20Varprat-1.png)

``` r
#Ground
A <- ggplot(df.venn.Ground, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
                  geom_circle(alpha = .6, size = 1, colour = "#3f3a3a",show.legend = FALSE) +
                  coord_fixed()+
                  annotate("text", x = df.venn.Ground$x , y = df.venn.Ground$y,label=df.venn.Ground$labels ,size = 6) +
                  scale_fill_manual(values = c("#82662a","#bea773","#daba73")) +
                  theme_void()
                  #annotate("text", x = 2 , y =1,label="A and B" ,size = 4) +
                  #annotate("text", x = 1.35 , y =2,label="B and C" ,size = 4) +
                  #annotate("text", x = 2.7 , y =2,label="A and C" ,size = 4) +
                  #annotate("text", x = 2 , y =1.6,label="A and B and C" ,size = 2)+theme_void()

A
```

![](Variation_Partitioning_files/figure-markdown_github/Ggplot%20Varprat-2.png)

``` r
#Canopy

B <- ggplot(df.venn.Canopy, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
                  geom_circle(alpha = .6, size = 1, colour = "#3f3a3a",show.legend = FALSE) +
                  coord_fixed()+
                  annotate("text", x = df.venn.Canopy$x , y = df.venn.Canopy$y,label=df.venn.Canopy$labels ,size = 6) +
                  scale_fill_manual(values = c("#42614a","#61a072","#8fd7a2")) +
                  theme_void()
                  #annotate("text", x = 2 , y =1,label="A and B" ,size = 4) +
                  #annotate("text", x = 1.35 , y =2,label="B and C" ,size = 4) +
                  #annotate("text", x = 2.7 , y =2,label="A and C" ,size = 4) +
                  #annotate("text", x = 2 , y =1.6,label="A and B and C" ,size = 2)+theme_void()

B
```

![](Variation_Partitioning_files/figure-markdown_github/Ggplot%20Varprat-3.png) Combine plots:

``` r
combi = ggarrange(A, B,
                  labels = c("A", "B"), 
                  ncol = 2, nrow = 1, font.label = list(size = 16, color = "black"),
                  common.legend = T, legend = "bottom", align = "h")+
        theme(plot.margin = margin(1,1,1,1, "cm")) 

combi
```

![](Variation_Partitioning_files/figure-markdown_github/combine%20Plots-1.png)

``` r
#save plot
#ggsave("VarPar_300_336x168.jpeg", plot = combi, 
#       device = "jpeg", dpi = 300, width = 336, height = 168, 
#       units = "mm")
#ggsave("VarPar_600_336x168.pdf", plot = combi, 
#       device = "pdf", dpi = 600, width = 336, height = 168, 
#       units = "mm")
#ggsave("VarPar_600_336x168.tif", plot = combi, 
#       device = "tiff", dpi = 600, width = 336, height = 168, 
#       units = "mm")
#ggsave("VarPar_600_336x168.png", plot = combi, 
#       device = "png", dpi = 600, width = 336, height = 168, 
#       units = "mm")
```

Variation Partitioning Analysis - calculation of adjusted R2 values
-------------------------------------------------------------------

The interpretation the results should be based on adjusted R2, although raw R2 is also reported (for CCA, adjusted R2 is calculated by permutation method of [Peres-Neto et al. 2006](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/0012-9658%282006%2987%5B2614%3AVPOSDM%5D2.0.CO%3B2?casa_token=tRcDjcXcbNwAAAAA%3AugcpGMzywSAYRFWifRSZJ-0fLtStFiHFuqDP1LnvoSLTdkJS_EtBhbnMCIjd0fl7nmpkonuTVIraoDafeQ) and may slightly vary between re-analyses of the same data; the argument permutations specifies the number of permutations used to calculate adjusted R2 in CCA).

``` r
rda.all<-rda(OTUCer.hel~Season+Microhabitat+TreeSpecies, data=SampleMetadata)
rda.all
R2all <- RsquareAdj(rda.all)$adj.r.squared
R2all #0.3552025


rda.MS<-rda(OTUCer.hel~Microhabitat+Season, data=SampleMetadata)
rda.MS 
R2MS <- RsquareAdj(rda.MS)$adj.r.squared
R2MS   #  0.3267351


rda.TS<-rda(OTUCer.hel~TreeSpecies+Season, data=SampleMetadata)
rda.TS
R2TS <- RsquareAdj(rda.TS)$adj.r.squared
R2TS   #  0.03582976


rda.MT<-rda(OTUCer.hel~Microhabitat+TreeSpecies, data=SampleMetadata)
rda.MT
R2MT <- RsquareAdj(rda.MT)$adj.r.squared
R2MT   #  0.341761

#Getting values for individual fractions, S(Season), M (Microhabitat), T(Tree species), substract: 
R2S<-(R2MS+R2TS)-R2all
R2S   #   0.007362318
R2M<-R2MS-R2S
R2M   #  0.3193728
R2T<-R2TS-R2S
R2T   # 0.02846744

# conditional effect of Microhabitat is highest  
```

Now, when we know both simple and conditional effect of each variables, we may want to know whether these variances are significant, and hence worth of interpreting. Results from varpart contain the column testable with logical values indicating whether given fraction is testable or not. To test each of them, we will need the models defined above, and the function anova, which (if applied on single object resulting from rda or cca method, returns Monte Carlo permutation test of the predictor effect). For this, we need to first define also partial ordination models with one variable as exlanatory and the other as covariable.

``` r
anova(rda.MS <- rda (OTUCer.hel ~ Microhabitat + Condition (Season), data = SampleMetadata))
#Model: rda(formula = OTUCer.hel ~ Microhabitat + Condition(Season), data = SampleMetadata)
#          Df Variance      F Pr(>F)    
#Model      8  0.13447 17.819  0.001 ***
#Residual 280  0.26413   

#anova(model <- capscale(OTUCer.hel~Season + Condition(Microhabitat), data = SampleMetadata))
#zeigt das selbe ergebnis wie mit der rda funktion
           
anova(rda.SM <- rda (OTUCer.hel ~ Season + Condition (Microhabitat), data = SampleMetadata))      
#Model: rda(formula = OTUCer.hel ~ Season + Condition(Microhabitat), data = SampleMetadata)
#          Df Variance      F Pr(>F)    
#Model      1 0.006217 6.5906  0.001 ***
#Residual 280 0.264131 


anova(rda.ST <- rda (OTUCer.hel ~ Season + Condition (TreeSpecies), data = SampleMetadata))
#Model: rda(formula = OTUCer.hel ~ Season + Condition(TreeSpecies), data = SampleMetadata)
#          Df Variance      F Pr(>F)    
#Model      1  0.00629 4.6534  0.001 ***
#Residual 286  0.38636                  


anova(rda.TS <- rda (OTUCer.hel ~ TreeSpecies + Condition (Season), data = SampleMetadata))
#Model: rda(formula = OTUCer.hel ~ TreeSpecies + Condition(Season), data = SampleMetadata)
#          Df Variance      F Pr(>F)    
#Model      2  0.01224 4.5298  0.001 ***
#Residual 286  0.38636                  


anova(rda.TM <- rda (OTUCer.hel ~ TreeSpecies + Condition (Microhabitat), data = SampleMetadata))
#Model: rda(formula = OTUCer.hel ~ TreeSpecies + Condition(Microhabitat), data = SampleMetadata)
#          Df Variance      F Pr(>F)    
#Model      2 0.013034 7.0663  0.001 ***
#Residual 279 0.257314  

anova(rda.MT <- rda (OTUCer.hel ~ Microhabitat + Condition (TreeSpecies), data = SampleMetadata))
#Model: rda(formula = OTUCer.hel ~ Microhabitat + Condition(TreeSpecies), data = SampleMetadata)
#          Df Variance      F Pr(>F)    
#Model      8  0.13534 18.343  0.001 ***
#Residual 279  0.25731  
       
anova(rda.all)
#Model: rda(formula = OTUCer.hel ~ Season + Microhabitat + TreeSpecies, data = SampleMetadata)
#          Df Variance      F Pr(>F)    
#Model     11  0.15377 15.473  0.001 ***
#Residual 278  0.25116       
```
