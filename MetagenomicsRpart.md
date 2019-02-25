# Data analysis in R
## Antti Karkman, 2017


Open R console or Rstudio and run the commands there.

```
library(vegan)
````

read in the data
```
birds_metaxa <- read.table("PATH/TO/birds_metaxa.txt", sep="\t", header=TRUE, row.names=1)
birds_card <- read.table("PATH/TO//birds_CARD.csv", sep=";", header=TRUE, row.names=1)
```
remove bad entries from the CARD annotations 
```
bad_card <- c("gb|AFH35853.1|ARO:3001328|Escherichia",
"gb|NP_414996.1|ARO:3004043|Escherichia",
"gb|NP_418813.1|ARO:3004109|Escherichia",
"gb|CCE36834|ARO:3003784|Mycobacterium",
"gb|AAO47226.2|ARO:3003318|Streptomyces",
"gb|YP_003971446|ARO:3003730|Bifidobacteria",
"gb|BAD59497.1|ARO:3000501|Nocardia",
"gb|CAA67349.1|ARO:3003359|Streptomyces")

birds_card <- birds_card[!(rownames(birds_card) %in% bad_card),]
```

normalise the metaxa results with the trimmed read count
```
birds_metaxa_norm <- sweep(birds_metaxa, 2, c(4.5, 3.3, 2.9, 2.7, 0.1, 2.5), '/')
```
normalise the ARG results with 16S counts from Metaxa2
```
birds_card_norm <- sweep(birds_card, 2, colSums(birds_metaxa), '/')
```

```
par(mfrow=c(1,2))
barplot(colSums(birds_metaxa), main="Raw counts")
barplot(colSums(birds_metaxa_norm), main="Normalised counts")
```

Ordination plots
```
card_dist <- vegdist(t(birds_card_norm))
card_ord <- cmdscale(card_dist)
plot(card_ord, pch=21, cex=3, bg=c(rep("red", 3), rep("blue", 3)), main="CARD PCoA")

card_dist <- vegdist(t(birds_metaxa_norm))
card_ord <- cmdscale(card_dist)
plot(card_ord, pch=21, cex=3, bg=c(rep("red", 3), rep("blue", 3)), main="Metaxa2 PCoA")
```

Most abundant microbes and ARGs
```
tail(birds_metaxa_norm[order(rowSums(birds_metaxa_norm)),])
tail(birds_card_norm[order(rowSums(birds_card_norm)),])
```
