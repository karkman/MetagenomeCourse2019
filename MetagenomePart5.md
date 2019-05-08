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
cd Metaphlan2
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
merge_metaphlan_tables.py *_metaphlan.txt > infants_metaphlan.txt
 ```

## HUMAnN2

```
#!/bin/bash -l
#SBATCH -J humann2
#SBATCH -o humann2_out_%A_%a.txt
#SBATCH -e humann2_err_%A_%a.txt
#SBATCH -t 2:00:00
#SBATCH --mem=10000
#SBATCH --array=1-10
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH -p serial


#module load biokit
source activate humann2_env
cd Humann2
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)
humann2 --input ../trimmed_data/$name"_R1_trimmed.fastq",../trimmed_data/$name"_R2_trimmed.fastq" \
              --input_type fastq --nproc  $SLURM_CPUS_PER_TASK \
              --mpa_pkl /appl/bio/metaphlan/db_v20/mpa_v20_m200.pkl \
              --bowtie2db /appl/bio/metaphlan/db_v20/mpa_v20_m200 \
              --bowtie2out $name".bowtie2.bz2" \
              -o $name"_metaphlan.txt"
```

## Optional

### StrainPhlAn
