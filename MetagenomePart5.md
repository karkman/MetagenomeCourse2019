# Metagenome analysis of infant gut metagenomes - part 5

## MetaPhlAn2


```
#!/bin/bash -l
#SBATCH -J metaphlan2
#SBATCH -o metaphlan2_out_%A_%a.txt
#SBATCH -e metaphlan2_err_%A_%a.txt
#SBATCH -t 2:00:00
#SBATCH --mem=10000
#SBATCH --array=1-10
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH -p serial

module load biokit
cd /wrk/antkark/Metagenomics2019/Metaphlan2
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)
metaphlan2.py ../trimmed_data/$name"_R1_trimmed.fastq",../trimmed_data/$name"_R2_trimmed.fastq" \
              --input_type fastq --nproc  $SLURM_CPUS_PER_TASK \
              --mpa_pkl /appl/bio/metaphlan/db_v20/mpa_v20_m200.pkl \
              --bowtie2db /appl/bio/metaphlan/db_v20/mpa_v20_m200 \
              --bowtie2out $name".bowtie2.bz2" \
              -o $name"_metaphlan.txt"
```
Merge the metaphlan2 outputs

```
module load biokit
merge_metaphlan_tables.py *_metaphlan.txt > infants_merged_table.txt
module purge
```

```
source activate metaphlan_plot_env
grep -E "(s__)|(^ID)" infants_merged_table.txt | grep -v "t__" | sed 's/^.*s__//g' > infants_metaphlan_species.txt
hclust2.py -i infants_metaphlan_species.txt -o abundance_heatmap_species.png --ftop 25 \
            --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 \
            -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 \
            --max_slabel_len 100 --minv 0.1 --dpi 300
 ```

 ## GraPhlAn

 ```
export2graphlan.py --skip_rows 1,2 -i infants_merged_table.txt --tree infants_merged.tree.txt \
                    --annotation merged_abundance.annot.txt --most_abundant 100 \
                     --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 \
                     --external_annotations 7 --min_clade_size 1
```
## HUMAnN2

```
#!/bin/bash -l
#SBATCH -J humann2
#SBATCH -o humann2_out_%A_%a.txt
#SBATCH -e humann2_err_%A_%a.txt
#SBATCH -t 1:00:00
#SBATCH --mem=2000
#SBATCH --array=1-10
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH -p serial

#module load biokit
source activate humann2_env
cd /wrk/antkark/Metagenomics2019/Humann2
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)
cat ../trimmed_data/$name"_R1_trimmed.fastq" ../trimmed_data/$name"_R2_trimmed.fastq" > $name"_tmp.fastq"
humann2 --input $name"_tmp.fastq" --output $name"_humann.txt"
rm $name"_tmp.fastq"
```

## Optional

### StrainPhlAn