# Metagenome analysis of infant gut metagenomes - part 2

## Antibiotic resistance gene annotation - reads
Next step is the antibiotic resistance gene annotation.  

We will map all the reads against ResFinder database to annote the antibiotic resistance genes present in the reads.  
First download the Resfinder database.  
Go to the ResFinder database [BitbBucket repository](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) and clone it to your work directory (`$WRKDIR`).   

```
cd $WRKDIR
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
```

This copies the whole repository to your work directory in a folder called `resfinder_db`. Go to that folder and concatenate the different antibiotic resistance gene sub-classes to one file
```
cat *.fsa > ResFinder.fasta
```

The annotation of resistance genes will be done as an *array* job. You can learn more about array jobs in Taito from [here.](https://research.csc.fi/fi/taito-array-jobs)   

We use Bowtie2 for the mapping and Samtools for analysing the mapping results. Both can be found from Taito, so we only need to load the biokit.  
But first we need to index the database file for mapping. This is done using `bowtie2-build` command. Bowtie2 part of the biokit, load it with `module load biokit`
```
bowtie2-build ResFinder.fasta ResFinder
```
This makes few files with the prefix `ResFinder` that Bowtie2 understands and uses in the mapping step.  

Then make the array job script and submit it as previously. (This takes some hours + possible time in the queue).  
```
#!/bin/bash -l
#SBATCH -J ARG_mapping
#SBATCH -o ARG_out_%A_%a.txt
#SBATCH -e ARG_err_%A_%a.txt
#SBATCH -t 05:00:00
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -p serial
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=1000
#

module load biokit
cd $WRKDIR/Metagenomics2019/

# get the n:th sample name
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p sample_names.txt)

# map fastq read pairs to Resfinder database and create bam files
bowtie2 -x $WRKDIR/resfinder_db/ResFinder -1 trimmed_data/$name"_R1_trimmed.fastq" -2  trimmed_data/$name"_R2_trimmed.fastq" \
        -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
        --threads $SLURM_CPUS_PER_TASK | samtools view -Sb - > $name.bam

# sort bam file and count reads which map in pairs and read pairs for which only one maps as one
samtools view -h $name".bam" | awk '$7!="=" || ($7=="=" && and($2,0x40)) {print $0}' | samtools view -Su  - \
        | samtools sort -o $name"_sort.bam"

# index the bam file
samtools index $name"_sort.bam"

# Name of the sample
echo -e $name > $name"_counts"

# take the counts from column 3
samtools idxstats $name"_sort.bam" | grep -v "*" | cut -f3 >> $name"_counts"
```

After the mapping every sample has its own result file. So we need to combine the results into one table.

```
#
module load biokit

# name the gene column in the result matrix
echo -e "GENE" > gene_names

# take the gene name from column 1
samtools idxstats 07004-B_sort.bam | grep -v "*" | cut -f1 >> gene_names

# create the final ARG genematrix
paste  gene_names *_counts > ARG_genemat.txt
```

After this inspect the results in the ARG gene matrix.  

## Assembly quality statistics
Let's take a look at the assembly file from yesterday. From the log file at `$WRKDIR/Metagenomics2019/co-assembly` you can check how the assembly went and from the last rows you can see some summary statistics of the assembly. However, for more detailed analysis we ran [MetaQUAST](http://bioinf.spbau.ru/metaquast) together with the assembly. Copy folder called "assembly_QC" to your computer. You can view the results (`report.html`) in your favorite browser.
