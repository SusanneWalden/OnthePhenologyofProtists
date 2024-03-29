# On the Phenology of Protists: Recurrent patterns reveal seasonal variation of protistan (Rhizaria: Cercozoa, Endomyxa) communities in tree canopies 


Welcome to the **On the Phenology of Protists** repository!

This repository is a collection of several scripts and mini-tutorials guiding you through the methods of metabarcoding analyses which were performed in the paper by [Walden et al., 2021](https://doi.org/10.1093/femsec/fiab081). The data and methods are mainly based on the paper by Jauss & Walden et al., 2020 and the corresponding Github repository [FromForestSoilToCanopy](https://github.com/RJauss/FromForestSoilToCanopy) .

The raw data can be downloaded [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB37525), plots and figures were generated with the final OTU Table (not provided) and annotation files accessible in the folder [00_Data](00_Data/). 

## Table of Content
*Relevant scripts and tutorials can be found in the corresponding markdown files, which are* **linked and highlighted** *for each category.*

### 00 Data
This folder contains the final OTU table, taxonomic and functional annotation, sample metadata as well as  supplementary material of this study.


### 01 Sequencing Results and Taxonomic and Funtional Annotation
The visualisation of the taxonomy then includes **[plotting the sequence similarity to reference sequences](01_Taxonomic_Functional_Annotation_and_Visualisation/Sequence_Similarity.md)** and a diagram showing the **[total taxonomic composition per microhabitat and sampling period](01_Taxonomic_Functional_Annotation_and_Visualisation/Taxonomic_Composition.md)** and the **[total functional trait composition per microhabitat and sampling period](01_Taxonomic_Functional_Annotation_and_Visualisation/Functional_Composition.md)**

### 02 Seasonal variation and spatial structuring
In this section we were looking for **[differentially abundundant OTUs per sampling season](02_Seasonal_Variation/DifferentiallyAbundantOTUs.md)** per sampling season in order to find seasonal patterns within detected protistan OTUs. In addition we provide the script for a **[Variation Partitioning Analysis](02_Seasonal_Variation/Variation_Partitioning.md)** between the season, microhabitats and tree species.

### 03 Determination of Alpha Diversity
Here we deal with the methods of **[plotting rarefaction curves](03_Alpha_Diversity/RarefactionCurves.md)** for earch sampling period and of course how to **[plot different alpha diversity indices in a combined boxplot](03_Alpha_Diversity/AlphaBoxplotGrouped.md)**.

### 04 Exploring Beta Diversity
One of the most straightforward methods of visualising beta diversity on the base of Bray-Curtis-Distances is an NMDS plot, the script is provided **[here](04_Beta_Diversity/NMDS.md)**. But we also perform and plot a **[Principle Coordinate Analysis](04_Beta_Diversity/PCoA.md)** to investigate differences between the foliar communites of different sampling periods, aswell as a **[Redundancy Analysis](04_Beta_Diversity/RDA_NewOrder.md)**.

### 05 Shared OTUs between Sampling Periods
This section deals with the **[visualisation of shared OTUs](05_Shared_OTUs/SharedOTUs.md)** of all investigated sampling periods.