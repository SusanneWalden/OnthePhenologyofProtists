
rm(list = ls())

library(vegan)
library(plyr)
library(ggplot2)
library(ggpubr)
library(UpSetR)
library(ggupset)
library(tidyverse)
library(cowplot)

setwd("C:/Users/Susi/Desktop/season_server/")
OTU_Table = as.data.frame(read.csv2("AllSeason_Cerc_En_Final_OTU_R_7633_forShared.csv",header = T))#jahres Zahlen geändert damit sie in der legende richtig abgebildet werden



SampleMetadata = OTU_Table[,1:7]# habe die metadaten noch um season+Jahr erweitert
Year = SampleMetadata$Year
OTU_Table = OTU_Table[,8:ncol(OTU_Table)]

####################

Aggregated_Year = ddply(OTU_Table, "Year", numcolwise(sum))
rownames(Aggregated_Year) = Aggregated_Year$Year
Aggregated_Year = Aggregated_Year[,-1]
# make incidence based
Aggregated_Year[,] = ifelse(Aggregated_Year[,] > 0, 1, 0)
t_Aggregated_Year = as.data.frame(t(Aggregated_Year))

# Plot

u = upset(t_Aggregated_Year, 
          nsets = 4, # Sets are the Years in this case
          order.by = c("degree","freq"), # choose between frequency or degree
          decreasing = T, group.by = "sets", 
          nintersects = 10, # number of combinations to display
          empty.intersections = TRUE, 
          mainbar.y.label = "Shared OTUs", 
          sets.x.label = "Number of OTUs per Year", 
          text.scale = 1.8, set_size.show = F)

u


################ Data preparation "Tidy"


Aggregated_Year_TF = ifelse(Aggregated_Year[,] > 0, 
                                    TRUE, FALSE)

tidy_Years = as_tibble(Aggregated_Year_TF, 
                               rownames = "Years") %>% 
  gather(OTU, Presence, -Years) %>% 
  filter(Presence) %>% 
  select(- Presence)

tidy_Dataframe = tidy_Years %>% 
  group_by(OTU) %>% 
  dplyr::summarise(shared_OTUs = list(Years))




############## plot


g = ggplot(tidy_Dataframe, aes(x = shared_OTUs)) +
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

g


###############

#### add barchart for total number of OTUs

tidy_Years$Years = factor(tidy_Years$Years, levels = rownames(as.data.frame(specnumber(Aggregated_Year) %>% sort())))

a = ggplot(tidy_Years, aes(x = Years)) + 
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
ga = ggdraw() + 
  draw_plot(g, x = 0.23, y = 0.06, height = 0.9, width = 0.8)+
  draw_plot(a, x = 0.04, y = 0, width = 0.2, height = 0.3)+### hier muss etwas mit der height gespielt werden deutschland 0.26, papua 0.3
  theme(plot.margin = margin(1,1,1,1, "cm")) 
ga


#write.csv(t_Aggregated_Year,file="shared_OTUs.csv")


########## save Plot


ggsave("Upset_300_336x168.jpeg", plot = ga, 
       device = "jpeg", dpi = 300, width = 336, height = 168, 
       units = "mm")
ggsave("Upset_600_336x168.pdf", plot = ga, 
       device = "pdf", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Upset_600_336x168.tif", plot = ga, 
       device = "tiff", dpi = 600, width = 336, height = 168, 
       units = "mm")
ggsave("Upset_600_336x168.png", plot = ga, 
       device = "png", dpi = 600, width = 336, height = 168, 
       units = "mm")





