# Day 1
*Antti Karkman, Katariina Pärnänen and Jenni Hultman*

# Metagenome analysis of infant gut metagenomes
Log into Taito, either with ssh (Mac/Linux) or PuTTy (Windows)

__All the scripts are to be run in `Metagenomics2019` folder!__

## Data download
First set up the course directory, make some folders and then download the data.  
Move to work directory at CSC and make a course directory there. Or a subfolder under a previous course folder.  
```
cd $WRKDIR
mkdir Metagenomics2019
cd Metagenomics2019
mkdir raw_data
mkdir scripts
mkdir trimmed_data
```

Download the metagenomic data (takes few minutes)  
```
cd raw_data
cp /wrk/parnanen/shared/COURSE_DATA/* .
```

Make a file containing the sample names. This is the second field separated by `_`.    
Use only the forward reads for this.  
 `ls *_R1.fastq.gz |awk -F "_" '{print $2}' > ../sample_names.txt`  

## QC and trimming
QC for the raw data (takes few min, depending on the allocation).  
Go to the folder that contains the raw data and make a folder called e.g. `FASTQC` for the QC reports.  
Then run the QC for the raw data and combine the results to one report using `multiqc`.  

Can be done on the interactive nodes using `sinteractive`.   
In that case do not allocate resources.  Just go to the right folder, activate `QC_env` and run `FastQC`. After it is finished you can just log out from the computing node. You don´t need to free any resources.

```
# allocate the computing resources and log in to the computing node.
salloc -n 1 --cpus-per-task=4 --mem=3000 --nodes=1 -t 00:30:00 -p serial
srun --pty $SHELL

# activate the QC environment
module load bioconda/3
source activate QC_env

# Run fastqc
fastqc ./*.fastq.gz -o FASTQC/ -t 4

# Then combine the reports with multiqc
multiqc ./ --interactive

# deactivate the virtual env
source deactivate

# log out from the computing node
exit

# and free the resources after the job is done
exit
```

Copy the resulting HTML file to your local machine with `scp` from the command line (Mac/Linux) or *WinSCP* on Windows. Have a look at the QC report with your favourite browser.  

After inspecting the output, it should be clear that we need to do some trimming.  
__What kind of trimming do you think should be done?__

You can check the adapter sequences from Illumina's website. (e.g. search with "*Illumina adapter sequences*"). Or from [here](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-10.pdf).

For trimming we'll make a bash script that runs `cutadapt` for each file using the `sample_names.txt` file.    
Go to your scripts folder and make a bash script for cutadapt with any text editor. Specify the adapter sequences that you want to trim after `-a` and `-A`. What is the difference with `-a` and `-A`? And what is specified with option `-O`? You can find the answers from Cutadapt [manual](http://cutadapt.readthedocs.io).

Save the file as `cutadapt.sh` in the scripts folder.  
__Make sure that the PATHS are correct in your own script__  

```
#!/bin/bash

while read i
do
        cutadapt  -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -O 10 --max-n 0 \
        -o ../trimmed_data/$i"_R1_trimmed.fastq" -p ../trimmed_data/$i"_R2_trimmed.fastq" \
        *$i*_R1.fastq.gz *$i*_R2.fastq.gz > ../trimmed_data/$i"_trim.log"
done < $1
```

Then we need a batch job file to submit the job to the SLURM system. More about CSC batch jobs here: https://research.csc.fi/taito-batch-jobs  
Make another file with text editor that runs the script above.

```
#!/bin/bash -l
#SBATCH -J cutadapt
#SBATCH -o cutadapt_out_%j.txt
#SBATCH -e cutadapt_err_%j.txt
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -p serial
#SBATCH --mem=100
#

module load biokit
cd $WRKDIR/Metagenomics2019/raw_data
bash ../scripts/cutadapt.sh ../sample_names.txt
```

After it is done, we can submit it to the SLURM system. Do it from the course main folder, so go one step back in your folders.  

`sbatch scripts/cut_batch.sh`  

You can check the status of your job with:  

`squeue -l -u $USER`  

After the job has finished, you can see how much resources it actually used and how many billing units were consumed. `JOBID` is the number after the batch job error and output files.  

`seff JOBID`  

Then let's check the results from the trimming.

Go to the folder containing the trimmed reads (`trimmed_data`) and make a new folder (`FASTQC`) for the QC files.  
Allocate some resources and then run FASTQC and MultiQC again.  

```
# Allocate resources
salloc -n 1 --cpus-per-task=4 --mem=3000 --nodes=1 -t 00:30:00 -p serial
srun --pty $SHELL

# activate the QC environment
module load bioconda/3
source activate QC_env

# run QC on the trimmed reads
fastqc ./*.fastq -o FASTQC/ -t 4
multiqc ./ --interactive

# deactivate the virtual env
source deactivate

# log out from the computing node
exit

# and free the resources after the job is done
exit
```

Copy it to your local machine as earlier and look how well the trimming went.  

## Assembly
We will assemble all 10 samples together (co-assembly) and use [Megahit assembler](https://github.com/voutcn/megahit) for the job. In addition, we will use MetaQuast to get some statistics about our assembly.  

Megahit is an ultra fast assembly tool for metagenomics data. It is installed to CSC and be loaded with following commands:
```
module purge
module load intel/16.0.0
module load megahit
```
module purge is needed to remove wrong Python versions you might have loaded earlier today.

Assembling metagenomic data can be very resource demanding and we need to do it as a batch job. As we want to do both individual and co-assemblies the R1 and R2 reads need to be merged into two files with `cat`

```
cd trimmed_data
cat *R1_trimmed.fastq > all_R1.fastq
cat *R2_trimmed.fastq > all_R2.fastq
```
Make a script called co_assembly.sh in a text editor
```
#!/bin/bash
#SBATCH -J megahit
#SBATCH -o megahit_out_%j.txt
#SBATCH -e megahit_err_%j.txt
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40000
#

module purge
module load intel/16.0.0
module load megahit

cd $WRKDIR/Metagenomics2019/

megahit -1 trimmed_data/all_R1.fastq -2 trimmed_data/all_R2.fastq \
         -o co-assembly -t $SLURM_CPUS_PER_TASK --min-contig-len 1000

# MetaQUAST assembly statistics
module purge
module load biokit
cd co-assembly
metaquast.py -t $SLURM_CPUS_PER_TASK --no-plots -o assembly_QC final.contigs.fa
```
Submit the batch job as previously

## HUMAnN2
HUMAnN2 might take some time to run.   
To be able to analyse the results on Thursday, we will run it already today.  

Make a batch job script for HUMAnN2.  
__NOTE! Use only the R1 reads AND the pre-downloaded databases, don't change the path to them__

```
#!/bin/bash -l
#SBATCH -J humann2
#SBATCH -o humann2_out_%A_%a.txt
#SBATCH -e humann2_err_%A_%a.txt
#SBATCH -t 20:00:00
#SBATCH --mem=20000
#SBATCH --array=1-10
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH -p serial

module load bioconda/3
source activate humann2_env
cd $WRKDIR/Metagenomics2019
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p sample_names.txt)

humann2 --input trimmed_data/$name"_R1_trimmed.fastq" --output Humann2 \
        --nucleotide-database /wrk/antkark/shared/chocophlan \
        --protein-database /wrk/antkark/shared/uniref \
        --metaphlan-options "--mpa_pkl /appl/bio/metaphlan/db_v20/mpa_v20_m200.pkl \
        --bowtie2db /appl/bio/metaphlan/db_v20/mpa_v20_m200"

```

## Antibiotic resistance gene annotation - reads
Next step is the antibiotic resistance gene annotation.  

We will map all the reads against ResFinder database to annote the antibiotic resistance genes present in the reads.  
First download the Resfinder database.  
Go to the ResFinder database [BitbBucket repository](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) and get the link to clone it to your work directory (`$WRKDIR`).   

```
cd $WRKDIR
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
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
bowtie2 -x $WRKDIR/resfinder_db/ResFinder -1 trimmed_data/$name"_R1_trimmed.fastq" \
        -2 trimmed_data/$name"_R2_trimmed.fastq" \
        -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
        --threads $SLURM_CPUS_PER_TASK | samtools view -Sb - > $name.bam

# sort bam file and count reads which map in pairs and read pairs for which only one maps as one
samtools view -h $name".bam" | awk '$7!="=" || ($7=="=" && and($2,0x40)) {print $0}' \
        | samtools view -Su  - \
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


# Optional

## (Fairly) Fast MinHash signatures with Sourmash

Go to the documentation of Sourmash and learn more about minhashes and the usage of Sourmash. https://sourmash.readthedocs.io/en/latest/

Before you start, make sure you're working at Taito-shell with command `sinteractive`.
```
sourmash compute *R1_trimmed.fastq -k 31 --scaled 10000
sourmash compare *.sig -o comparisons
sourmash plot comparisons
# annotate one
sourmash gather 09069-B_R1_trimmed.fastq.sig /wrk/antkark/shared/genbank-d2-k31.sbt.json -o OUTPUT_sour.txt
```

## Taxonomic profiling with Metaxa2

The microbial community profiling for the samples can alsp be done using a 16S/18S rRNA gene based classification software [Metaxa2](http://microbiology.se/software/metaxa2/).  
It identifies the 16S/18S rRNA genes from the short reads using HMM models and then annotates them using BLAST and a reference database.
We will run Metaxa2 as an array job in Taito. More about array jobs at CSC [here](https://research.csc.fi/taito-array-jobs).  
Make a folder for Metaxa2 results and direct the results to that folder in your array job script. (Takes ~6 h for the largest files)

```
#!/bin/bash -l
#SBATCH -J metaxa
#SBATCH -o metaxa_out_%A_%a.txt
#SBATCH -e metaxa_err_%A_%a.txt
#SBATCH -t 10:00:00
#SBATCH --mem=15000
#SBATCH --array=1-10
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH -p serial

cd $WRKDIR/Metagenomics2019/Metaxa2
# Metaxa uses HMMER3 and BLAST, so load the biokit first
module load biokit
# each job will get one sample from the sample names file stored to a variable $name
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)
# then the variable is used in running metaxa2
metaxa2 -1 ../trimmed_data/$name"_R1_trimmed.fastq" -2 ../trimmed_data/$name"_R2_trimmed.fastq" \
            -o $name --align none --graphical F --cpu $SLURM_CPUS_PER_TASK --plus
metaxa2_ttt -i $name".taxonomy.txt" -o $name
```

When all Metaxa2 array jobs are done, we can combine the results to an OTU table. Different levels correspond to different taxonomic levels.  
When using any 16S rRNA based software, be cautious with species (and beyond) level classifications. Especially when using short reads.  
We will look at genus level classification.
```
# Genus level taxonomy
cd Metaxa2
metaxa2_dc -o metaxa_genus.txt *level_6.txt
```
