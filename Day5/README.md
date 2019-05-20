# HUMAnN2
Prepare environment and go to the humann2 output folder
```
module load bioconda/3
source activate humann2_env
cd $WRKDIR/Metagenomics2019/Humann2
```
See how humann2_join_tables works
```
humann2_join_tables -h
```

Combine gene family profiles
```
humann2_join_tables -i . \
  --file_name genefamilies.tsv \
  -o genefamilies.tsv
```

Look at gene family table
```
less -S genefamilies.tsv
```
 - value are RPKs
 - first line measures unmapped reads
 - each gene family is stratified to contributing organisms

Generate a table with no stratifications.  
(for any analysis where you want to look at gene family abundances but don't care contributing organisms)  

Note the difference in file size; for large data sets this makes a big difference
```
grep -v "|" genefamilies.tsv > genefamilies_nostratification.tsv
```

Add gene family names
```
humann2_rename_table -h
humann2_rename_table -i genefamilies_nostratification.tsv \
  -n uniref90 \
  -s \
  -o genefamilies_nostratification_names.tsv
```

Inspect the results
```
less -S genefamilies_nostratification_names.tsv
```

Generate new "grouping": MetaCyc enzymes.  
(examples without the stratifications work speed; similar commands work for files with stratifications)
```
humann2_regroup_table -h
humann2_regroup_table -i genefamilies_nostratification.tsv \
  -g uniref90_level4ec \
  -o level4ec.tsv
```
Add names
```
humann2_rename_table -i level4ec.tsv \
  -n ec \
  -s \
  -o level4ec_names.tsv
```
Normalize to "counts per million" (CPM)
```
humann2_renorm_table -h
humann2_renorm_table -i level4ec_names.tsv \
  -p \
  -o level4ec_names_norm.tsv

less -S level4ec_names_norm.tsv
```
SAMPLE_ID_pathabundance.tsv output files quantify MetaCyc pathway abundances  
Combine profiles in one table and normalize (this time with the organisms level stratifications)

```
humann2_join_tables -i . \
  --file_name pathabundance.tsv \
  -o pathways.tsv

humann2_renorm_table -i pathways.tsv \
  -p \
  -o pathways_norm.tsv
```

Analyse pathway table in R (in separate .R file)

.. visualize pathway identified in R analysis

```
humann2_barplot -h

humann2_barplot -i pathways_norm.tsv \
  -f PWY0-162 \
  -s usimilarity \
  -o PWY0-162.pdf \
  -a pseudolog \
```

Different normalization and more species for the second one
```
humann2_barplot -i pathways_norm.tsv \
  -f PWY-7219 \
  -s usimilarity \
  -o PWY-7219.pdf \
  -a normalize \
  -t 15
```
## Extra assignment if you have time:
Generate combined table of MetaCyc pathways and visualize pathways with high and low contributional diversities (one each).
