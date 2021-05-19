Shared OTUs
================

To compare species richness of different sampling periods, we calculated shared OTUs.

Load data
---------

``` r
rm(list = ls())

library(vegan)
library(plyr)
library(ggplot2)
library(ggpubr)
library(UpSetR)
library(ggupset)
library(tidyverse)
library(cowplot)


OTU_Table = as.data.frame(read.csv2("../00_Data/05_Cercozoa_Seasonal_OTU_Table_min-freq-7633_transposed_withMetadata.csv",header = T))



SampleMetadata = OTU_Table[,1:18]
Period = SampleMetadata$SamplingDescription
OTU_Table = OTU_Table[,18:ncol(OTU_Table)]

#Aggregate by Sampling Period
Aggregated_Period = ddply(OTU_Table, "Period", numcolwise(sum))
rownames(Aggregated_Period) = Aggregated_Period$Period
Aggregated_Period = Aggregated_Period[,-1]

#make incidence based
Aggregated_Period[,] = ifelse(Aggregated_Period[,] > 0, 1, 0)
t_Aggregated_Period = as.data.frame(t(Aggregated_Period))
```

UpSetR
------

Now we want to visualise shared OTUs between each sampling period using the UpSetR package.

``` r
U = upset(t_Aggregated_Period, 
          nsets = 4, # Sets are the sampling periods in this case
          order.by = c("degree","freq"), # choose between frequency or degree
          decreasing = T, group.by = "sets", 
          nintersects = 6, # number of combinations to display
          empty.intersections = TRUE, 
          mainbar.y.label = "Shared OTUs", 
          sets.x.label = "Number of OTUs per Year", 
          text.scale = 1.8, set_size.show = F)

U
```

![](SharedOTUs_files/figure-markdown_github/UpsetSharedOTUs-1.png)

Data Preparation for ggupset
----------------------------

As an input, we need a tibble (another version of a dataframe). To generate a tibble from a presence/absence matrix, we need to convert it with using the tidyverse package.

``` r
#sort data
Aggregated_Period_TF = ifelse(Aggregated_Period[,] > 0, 
                                    TRUE, FALSE)

tidy_Periods = as_tibble(Aggregated_Period_TF, 
                               rownames = "Periods") %>% 
  gather(OTU, Presence, -Periods) %>% 
  filter(Presence) %>% 
  select(- Presence)

tidy_Dataframe = tidy_Periods %>% 
  group_by(OTU) %>% 
  dplyr::summarise(shared_OTUs = list(Periods))
```

ggupset
-------

Now we want to visualise shared OTUs between each sampling period using the ggupset package.

``` r
a = ggplot(tidy_Dataframe, aes(x = shared_OTUs)) +
  # plot the shard OTUs
  geom_bar() +
  geom_label(aes(label=..count..), 
             stat="count", 
             position=position_stack(),
             size=5) +
  # convert to combination matrix
  scale_x_upset() +
  # display only 8 intersections
  axis_combmatrix(xlim = c(0, 8.35)) +
  theme_combmatrix(combmatrix.label.make_space = F, 
                   combmatrix.label.text = element_text(size=16, 
                                                        face = "bold"),
                   axis.text.y = element_text(size = 16), 
                   axis.title=element_text(size=16, face = "bold"), 
                   panel.background = element_blank(), 
                   panel.grid.major.y = element_line(colour = "black"), 
                   panel.grid.minor.y = element_line(color = "black"), 
                   panel.grid.major.x = element_blank(), 
                   panel.grid.minor.x = element_blank(), 
                   plot.margin = unit(c(0.5,0.5,0.5,1.8), "cm"))+ 
  labs(x = NULL, y = "Number of shared OTUs")

a
```

![](SharedOTUs_files/figure-markdown_github/ggupsetSharedOTUs-1.png)

Add a seperate Barchart to the Plot
-----------------------------------

Now we want to add the missing left bar chart of the number of OTUs.

``` r
tidy_Periods$Periods = factor(tidy_Periods$Periods, levels = rownames(as.data.frame(specnumber(Aggregated_Period) %>% sort())))

b = ggplot(tidy_Periods, aes(x = Periods)) + 
  geom_bar() + 
  coord_flip() +  
  scale_y_reverse() + 
  theme_minimal() + 
  labs(y = "Number of OTUs", x = NULL) + 
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size = 16), 
        axis.text.x = element_text(size = 15), 
        axis.ticks.x = element_line(),
        axis.line.x = element_line(),
        axis.line.y = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        plot.background = element_blank())
ab = ggdraw() + 
  draw_plot(a, x = 0.23, y = 0.06, height = 0.9, width = 0.8)+
  draw_plot(b, x = 0.03, y = 0, width = 0.2, height = 0.23)+
  theme(plot.margin = margin(1,1,1,1, "cm")) 
ab
```

![](SharedOTUs_files/figure-markdown_github/Add%20Barchart-1.png)
