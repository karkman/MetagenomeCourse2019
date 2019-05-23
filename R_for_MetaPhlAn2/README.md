MetaPhlAn2
================
Tommi

MetaPhlAn2 analysis and visualization in R
==========================================

Set your working directory to where you have your data on your own computer (and install) and load needed libraries
-------------------------------------------------------------------------------------------------------------------

download data from <https://www.dropbox.com/s/agzc1fecgoj708z/infants_merged_table.txt?dl=0>

``` r
setwd("~/Dropbox/teaching/HY_metagenomic_course_2019/MetagenomeCourse2019/R_for_MetaPhlAn2/")
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("devtools")

library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
    ## ✔ tibble  1.4.2     ✔ dplyr   0.7.4
    ## ✔ tidyr   0.8.0     ✔ stringr 1.4.0
    ## ✔ readr   1.1.1     ✔ forcats 0.3.0

    ## Warning: package 'stringr' was built under R version 3.5.2

    ## ── Conflicts ─────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-1

``` r
library(devtools)
```

    ## 
    ## Attaching package: 'devtools'

    ## The following object is masked from 'package:permute':
    ## 
    ##     check

``` r
# Use this to install microbiomics if needed:
#install_github("tvatanen/microbiomics")
library(microbiomics)
```

Read in the species from MetaPhlan table
----------------------------------------

``` r
metaphlan_species <- read_metaphlan_table("infants_merged_table.txt", lvl = 7, normalize = T)

rownames(metaphlan_species) <- sapply(rownames(metaphlan_species), function(x) strsplit(x, ".", fixed = T)[[1]][1])

mds_obj <- metaMDS(metaphlan_species)
```

    ## Run 0 stress 0.1564851 
    ## Run 1 stress 0.1922584 
    ## Run 2 stress 0.1361479 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1521212  max resid 0.359429 
    ## Run 3 stress 0.1686252 
    ## Run 4 stress 0.1720604 
    ## Run 5 stress 0.1361479 
    ## ... Procrustes: rmse 7.598371e-06  max resid 1.567332e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.1879219 
    ## Run 7 stress 0.1686306 
    ## Run 8 stress 0.1643489 
    ## Run 9 stress 0.1361479 
    ## ... New best solution
    ## ... Procrustes: rmse 1.044081e-06  max resid 1.81197e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.1686373 
    ## Run 11 stress 0.1720604 
    ## Run 12 stress 0.1723213 
    ## Run 13 stress 0.1564851 
    ## Run 14 stress 0.1698597 
    ## Run 15 stress 0.1564851 
    ## Run 16 stress 0.1723214 
    ## Run 17 stress 0.1686305 
    ## Run 18 stress 0.151936 
    ## Run 19 stress 0.213171 
    ## Run 20 stress 0.1564851 
    ## *** Solution reached

``` r
data.frame(mds_obj$points) %>%
  rownames_to_column("sampleID") %>%
  ggplot(aes(x=MDS1, y=MDS2, color = sampleID)) + 
  geom_point() +
  coord_equal() +
  theme_bw()
```

![](README_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
metaphlan_species_long <-
  metaphlan_species %>%
  rownames_to_column("sampleID") %>%
  gather(taxon_name, relative_abundance, -sampleID) %>%
  separate(taxon_name, sep = "\\.", into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
  
species_stats <- 
  metaphlan_species_long %>%
  group_by(species) %>%
  summarize(mean_relative_abundance = mean(relative_abundance),
            median_relative_abundance = median(relative_abundance),
            max_relative_abundance = max(relative_abundance),
            prevalence = sum(relative_abundance > 0) / n())

head(species_stats %>% arrange(-mean_relative_abundance))
```

    ## # A tibble: 6 x 5
    ##   species   mean_relative_ab… median_relative… max_relative_ab… prevalence
    ##   <chr>                 <dbl>            <dbl>            <dbl>      <dbl>
    ## 1 s__Esche…            0.238          0.0887              0.654        0.7
    ## 2 s__Esche…            0.141          0.0616              0.367        0.8
    ## 3 s__Veill…            0.0997         0.000247            0.530        0.6
    ## 4 s__Klebs…            0.0917         0                   0.912        0.4
    ## 5 s__Bifid…            0.0779         0                   0.641        0.4
    ## 6 s__Actin…            0.0599         0                   0.599        0.1

Generate barplot of 10 most abundant species
--------------------------------------------

``` r
species_stats %>% 
  arrange(-mean_relative_abundance) %>%
  top_n(10) %>% 
  left_join(metaphlan_species_long) %>%
  ggplot(aes(y=relative_abundance, x=sampleID, fill = species)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Relative abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

    ## Selecting by prevalence

    ## Joining, by = "species"

![](README_files/figure-markdown_github/barplot-1.png)

GeneratebBarplot of 10 most abundant genera
===========================================

``` r
metaphlan_genera_long <- 
  metaphlan_species_long %>%
  group_by(genus,sampleID) %>%
  summarise(relative_abundance = sum(relative_abundance)) %>%
  ungroup()
  
metaphlan_genera_long %>%
  group_by(genus) %>%
  summarize(mean_relative_abundance = mean(relative_abundance)) %>%
  arrange(-mean_relative_abundance) %>%
  top_n(10) %>%
  left_join(metaphlan_genera_long) %>%
  ggplot(aes(y=relative_abundance, x=sampleID, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Relative abundance")
```

    ## Selecting by mean_relative_abundance

    ## Joining, by = "genus"

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)
