Humann2
================
Tommi

Humann2 analysis and visualization in R
=======================================

Set your working directory to where you have your data on your own computer (and install) and load needed libraries
-------------------------------------------------------------------------------------------------------------------

``` r
setwd("/Users/kparnane/Documents/MetagenomeCourse2019/R_for_Humann2/")
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("broom")

library(tidyverse)
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

    ## -- Attaching packages ---------------------------------- tidyverse 1.2.1 --

    ## v ggplot2 3.1.1     v purrr   0.3.2
    ## v tibble  2.1.1     v dplyr   0.8.1
    ## v tidyr   0.8.3     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.4.0

    ## -- Conflicts ------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-5

``` r
library(broom)
```

Read in the normalized pathways file
------------------------------------

``` r
humann2_pathways <- as.data.frame(read_tsv("pathways_norm.tsv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   `# Pathway` = col_character(),
    ##   `02078-BA_R1_trim_Abundance-CPM` = col_double(),
    ##   `07004-B_R1_trim_Abundance-CPM` = col_double(),
    ##   `07005-B_R1_trim_Abundance-CPM` = col_double(),
    ##   `08012-B_R1_trim_Abundance-CPM` = col_double(),
    ##   `08016-B_R1_trim_Abundance-CPM` = col_double(),
    ##   `09024-B_R1_trim_Abundance-CPM` = col_double(),
    ##   `09069-B_R1_trim_Abundance-CPM` = col_double(),
    ##   `10029-BA_R1_trim_Abundance-CPM` = col_double(),
    ##   `10030-B_R1_trim_Abundance-CPM` = col_double(),
    ##   `12042-B_R1_trim_Abundance-CPM` = col_double()
    ## )

Re-format the data frame from wide to long format
-------------------------------------------------

tidyverse is a really handy package for data wrangling and reformatting tables. This block is commented for what each line does. You can check the help page for each command by typing ? and the command's name for exaple ?rename for the next blocks of code

``` r
humann2_pathways_long <- 
  # %>% is used to pipe the output to the next command
  humann2_pathways %>% 
  #Rename  `# Pathway` as pathway
  rename(pathway = `# Pathway`) %>%
  #Gather cmp by pathway and sampleID
  gather(sampleID, cpm, -pathway) %>% 
  #Separate by sampleID and drop any extra values without warning
  separate(sampleID, "sampleID", sep = "_", extra = "drop") %>% 
  #Separate pathways from organisms using |
  separate(pathway, c("pathway", "organism"), sep = "\\|", fill = "right")
```

Continue processing the humann2 output
======================================

``` r
# generate pathway table with no organism stratifications
humann2_pathways_no_stratifications_long <-
  humann2_pathways_long %>%
  filter(is.na(organism)) %>%
  select(-organism) %>%
  filter(!(grepl("^UN", pathway))) 

# Compute pathway alpha divertities per sample
humann2_pathways_no_stratifications_long %>%
  group_by(sampleID) %>%
  summarise(shannons_div = vegan::diversity(cpm),
            num_pathways = sum(cpm>0)) 
```

    ## # A tibble: 10 x 3
    ##    sampleID shannons_div num_pathways
    ##    <chr>           <dbl>        <int>
    ##  1 02078-BA         5.35          273
    ##  2 07004-B          5.43          281
    ##  3 07005-B          4.98          160
    ##  4 08012-B          4.39          186
    ##  5 08016-B          5.29          274
    ##  6 09024-B          5.08          255
    ##  7 09069-B          4.92          185
    ##  8 10029-BA         5.45          289
    ##  9 10030-B          4.83          301
    ## 10 12042-B          5.39          297

``` r
# .. continue with any statistical comparisons etc.

# Work with organism level stratifications
humann2_pathway_stratifications_long <- 
  humann2_pathways_long %>%
  filter(!(is.na(organism))) %>%
  filter(!(grepl("^UN", pathway))) 

# number of organisms per pathway
humann2_organisms_per_pathway <- 
  humann2_pathway_stratifications_long %>%
  group_by(pathway) %>%
  summarise(num_organisms = length(unique(organism)))

# average contributional alpha diverity (Gini simpson diversity) per pathway
humann2_pathway_alpha_div <- 
  humann2_pathway_stratifications_long %>%
  filter(cpm > 0) %>%
  group_by(pathway, sampleID) %>%
  summarise(alpha_div = vegan::diversity(cpm, index = "simpson")) %>%
  group_by(pathway) %>%
  summarise(mean_alpha_div = mean(alpha_div),
            median_alpha_dv = median(alpha_div)) %>%
  arrange(-mean_alpha_div)

head(humann2_pathway_alpha_div)
```

    ## # A tibble: 6 x 3
    ##   pathway                                    mean_alpha_div median_alpha_dv
    ##   <chr>                                               <dbl>           <dbl>
    ## 1 PWY-7219: adenosine ribonucleotides de no~          0.510           0.527
    ## 2 PWY-7221: guanosine ribonucleotides de no~          0.465           0.475
    ## 3 PWY-7220: adenosine deoxyribonucleotides ~          0.418           0.480
    ## 4 PWY-7222: guanosine deoxyribonucleotides ~          0.418           0.480
    ## 5 PWY-6386: UDP-N-acetylmuramoyl-pentapepti~          0.381           0.480
    ## 6 PWY-6126: superpathway of adenosine nucle~          0.380           0.486

``` r
# add number of samples where pathway present
humann2_pathway_alpha_div_with_n_samples <- 
  humann2_pathways_no_stratifications_long %>%
  filter(cpm > 0) %>%
  group_by(pathway) %>%
  summarise(n_samples = n()) %>%
  left_join(humann2_pathway_alpha_div) %>%
  left_join(humann2_organisms_per_pathway)
```

    ## Joining, by = "pathway"
    ## Joining, by = "pathway"

``` r
humann2_pathway_alpha_div_with_n_samples %>%
  filter(n_samples == 10) %>%
  arrange(mean_alpha_div) %>%
  print(n = 20)
```

    ## # A tibble: 121 x 5
    ##    pathway           n_samples mean_alpha_div median_alpha_dv num_organisms
    ##    <chr>                 <int>          <dbl>           <dbl>         <int>
    ##  1 ANAEROFRUCAT-PWY~        10        0               0                   1
    ##  2 DENOVOPURINE2-PW~        10        0               0                   1
    ##  3 GLUCONEO-PWY: gl~        10        0               0                   2
    ##  4 HEXITOLDEGSUPER-~        10        0               0                   1
    ##  5 P4-PWY: superpat~        10        0               0                   1
    ##  6 PWY-5154: L-argi~        10        0               0                   1
    ##  7 PWY-5265: peptid~        10        0               0                   1
    ##  8 PWY-6270: isopre~        10        0               0                   2
    ##  9 PWY66-389: phyto~        10        0               0                   1
    ## 10 PYRIDOXSYN-PWY: ~        10        0               0                   2
    ## 11 THISYNARA-PWY: s~        10        0               0                   1
    ## 12 GLYCOLYSIS-E-D: ~        10        0.00106         0                   2
    ## 13 GLYCOLYSIS: glyc~        10        0.00122         0                   3
    ## 14 PWY-841: superpa~        10        0.00212         0                   3
    ## 15 PWY0-162: superp~        10        0.00428         0                   6
    ## 16 ASPASN-PWY: supe~        10        0.00786         0.00101             5
    ## 17 NONMEVIPP-PWY: m~        10        0.0114          0                   6
    ## 18 PWY-7388: octano~        10        0.0119          0                   7
    ## 19 PWY66-422: D-gal~        10        0.0142          0                   5
    ## 20 PWY-7234: inosin~        10        0.0191          0                   8
    ## # ... with 101 more rows

``` r
# PWY0-162 as an example of low diversity pathway

humann2_pathway_alpha_div_with_n_samples %>%
  filter(n_samples == 10) %>%
  arrange(-mean_alpha_div)
```

    ## # A tibble: 121 x 5
    ##    pathway           n_samples mean_alpha_div median_alpha_dv num_organisms
    ##    <chr>                 <int>          <dbl>           <dbl>         <int>
    ##  1 PWY-7219: adenos~        10          0.510           0.527            34
    ##  2 PWY-7221: guanos~        10          0.465           0.475            28
    ##  3 PWY-7220: adenos~        10          0.418           0.480            17
    ##  4 PWY-7222: guanos~        10          0.418           0.480            17
    ##  5 PWY-6386: UDP-N-~        10          0.381           0.480            19
    ##  6 PWY-6126: superp~        10          0.380           0.486            15
    ##  7 PWY-6122: 5-amin~        10          0.379           0.464            26
    ##  8 PWY-6277: superp~        10          0.379           0.464            26
    ##  9 PWY-7208: superp~        10          0.376           0.441            19
    ## 10 PWY-7228: superp~        10          0.373           0.424            21
    ## # ... with 111 more rows

``` r
# PWY-7219 as an example of high divertsity pathway
```
