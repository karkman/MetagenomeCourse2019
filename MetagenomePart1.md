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
cp /wrk/antkark/shared/Metagenomics2019data.tar.gz .
```

The md5 sum for the file is `2d86933b525034dcdcddf9592bc0e05c`. Check that the md5 sum of the file you downloaded matches with this.

```
md5sum Metagenomics2019data.tar.gz
```
And then unpack the zipped tar file with `tar`.  
The options are: `x` = extract, `z` = unzip, `v` = verbose, `f` = file

```
tar -xzvf Metagenomics2019data.tar.gz
```

Make a file containing the sample names. This is the second field separated by `_`.    
Use only the forward reads for this.  
 `ls *_R1_001.fastq |awk -F "_" '{print $2}' > ../sample_names.txt`  

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
fastqc ./*.fastq -o FASTQC/ -t 4

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
        *$i*_R1_001.fastq *$i*_R2_001.fastq > ../trimmed_data/$i"_trim.log"
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
cat *R1_trimmed.fastq > all_R1_trimmed.fastq
cat *R2_trimmed.fastq > all_R2_trimmed.fastq
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

megahit -1 trimmed_data/all_R1_trimmed.fastq -2 trimmed_data/all_R2_trimmed.fastq \
         -o co-assembly -t $SLURM_CPUS_PER_TASK --min-contig-len 1000

# MetaQUAST assembly statistics
module purge
module load biokit
cd co-assembly
metaquast.py -t $SLURM_CPUS_PER_TASK --no-plots -o assembly_QC final.contigs.fa
```
Submit the batch job as previously

Things to add:
- Humann2
- ARG annotation w/ mapping
-

**Optional**

## (Fairly) Fast MinHash signatures with Sourmash
```
sourmash compute *R1_trimmed.fastq -k 31 --scaled 1000
sourmash compare *.sig -o comparisons
sourmash plot comparisons
sourmash gather SIGNATURE.sig ../../shared/genbank-d2-k31.sbt.json -o OUTPUT_sour.txt
```

## Taxonomic profiling with Metaxa2

The microbial community profiling for the samples will be done using a 16S/18S rRNA gene based classification software [Metaxa2](http://microbiology.se/software/metaxa2/).  
It identifies the 16S/18S rRNA genes from the short reads using HMM models and then annotates them using BLAST and a reference database.
We will run Metaxa2 as an array job in Taito. More about array jobs at CSC [here](https://research.csc.fi/taito-array-jobs).  
Make a folder for Metaxa2 results and direct the results to that folder in your array job script. (Takes ~6 h for the largest files)

```
#!/bin/bash -l
#SBATCH -J metaxa
#SBATCH -o metaxa_out_%A_%a.txt
#SBATCH -e metaxa_err_%A_%a.txt
#SBATCH -t 10:00:00
#SBATCH --mem=10000
#SBATCH --array=1-6
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
