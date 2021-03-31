Sequence Similarity
================

One very first step is to check the sequence similarity of our OTUs to the used reference database (PR2 database). Basically we plot the percentage of identity from our taxonomic annotation against the number of sequences, as well as OTUs in a histogram.

**[The PR2 sequence database](https://pr2-database.org/)** was initiated in 2010 in the frame of the BioMarks project from work that had developed in the previous ten years in the Plankton Group of the Station Biologique of Roscoff. Its aim is to provide a reference database of carefully annotated 18S rRNA sequences using eight unique taxonomic fields (from kingdom to species).

Load Data
---------

For this script we need two files: the taxonomy file, giving the percent identity of each detected OTU with the best hit within the used reference database, while the OTU table will be used to calculate the total number sequence abundances:

``` r
rm(list = ls())

library(ggplot2)
library(plyr)
library(ggpubr)


TAX= read.csv2("../00_Data/AllSeason_Cerc_En_Taxonomy_R.csv",header = T)
OTU_Table = as.data.frame(read.csv2("../00_Data/05_Cercozoa_Seasonal_OTU_Table_min-freq-7633_transposed_withMetadata.csv",header = T))

SampleMetadata = OTU_Table[,1:17]
OTU_Table = OTU_Table[,18:ncol(OTU_Table)]


Abundances = colSums(OTU_Table)
TAX = cbind(TAX, Abundances)

TAX$OTU_ID = paste0("OTU", TAX$OTU_Number, "_", TAX$Species)
```

In the next step we aggregate all reads depending on the percentage of their sequence similarity and expand the table to make it readable for *geom\_histogram*.

``` r
# Aggregate Abundances with the same percent identity
AggregatedTAXpercent = ddply(TAX, "PercentID", numcolwise(sum))
AggregatedTAXpercent$PercentID = AggregatedTAXpercent$PercentID / 100

# Expand the table. 
AggregatedTAXpercent_expanded = AggregatedTAXpercent[rep(row.names(AggregatedTAXpercent), AggregatedTAXpercent$Abundance), c(1,2)]
```

Read similarities to reference database
---------------------------------------

``` r
g_reads = ggplot(AggregatedTAXpercent_expanded[AggregatedTAXpercent_expanded > 0, ], aes(x = PercentID)) +
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
  guides(fill = guide_legend(override.aes = list(alpha = 0.3)))

g_reads
```

![](Sequence_Similarity_files/figure-markdown_github/Plot%20reads-1.png)

OTU similarities to reference database
--------------------------------------

Now, we want to take a look at the similarities of detected OTUs without taking read abundances into account:

``` r
OTU_TAX = subset(TAX, select = "PercentID")
OTU_TAX$Abundances = TAX$Abundances
OTU_TAX$PercentID = OTU_TAX $PercentID / 100



g_OTU = ggplot(OTU_TAX[OTU_TAX $PercentID > 0, ], aes(x = PercentID)) +
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
  guides(fill = guide_legend(override.aes = list(alpha = 0.3)))


g_OTU 
```

![](Sequence_Similarity_files/figure-markdown_github/Plot%20OTUs-1.png)

Next, we combine both plots:

``` r
combi = ggarrange(g_reads,g_OTU,
                  labels = c("A", "B"), 
                  ncol = 2, nrow = 1, font.label = list(size = 16, color = "black"),
                  common.legend = T, legend = "right", align = "h")+
        theme(plot.margin = margin(1,1,1,1, "cm")) 

combi
```

![](Sequence_Similarity_files/figure-markdown_github/combine%20Plots-1.png)

``` r
# save plot 
 
#ggsave("Similarity_300_336x168.jpeg", plot = combi, 
#       device = "jpeg", dpi = 300, width = 336, height = 168, 
#       units = "mm")
#ggsave("Similarity_600_336x168.pdf", plot = combi, 
#       device = "pdf", dpi = 600, width = 336, height = 168, 
#       units = "mm")
#ggsave("Similarity_600_336x168.tif", plot = combi, 
#       device = "tiff", dpi = 600, width = 336, height = 168, 
#       units = "mm")
#ggsave("Similarity_600_336x168.png", plot = combi, 
#       device = "png", dpi = 600, width = 336, height = 168, 
#       units = "mm")
```
