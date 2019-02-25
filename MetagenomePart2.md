*Jenni Hultman and Antti Karkman*

# Metagenome part 2

## Antibiotic resistance gene annotation - reads
Next step is the antibiotic resistance gene annotation.  

First convert all R1 files from fastq to fasta. You can use `fastq_to_fasta_fast` program that is included in the biokit.  
`fastq_to_fasta_fast TRIMMED.fastq > TRIMMED.fasta`  

Then add the sample name to each sequence header after the `>` sign. Keep the original fasta header and separate it with `-`. You can do this either separately for each file or write a batch script to go through all files .You can use the `sample_names.txt` as a input to your bash script.  
```
# Example header
>L3-M01457:76:000000000-BDYH7:1:1101:17200:1346
```

When all R1 files have been converted to fasta and renamed, combine them to one file. It will be the input file for antibiotic resistance gene annotation. The sample name in each fasta header will make it possible to count the gene abundances for each sample afterwards using the script you downloaded earlier.  

We will annotate the resistance genes using The Comprehensive Antibiotic Resistance Database, [CARD.](https://card.mcmaster.ca)  
Go to the CARD website and download the latest CARD release to folder called CARD under your user applications (`$USERAPPL`) folder and unpack the file.  

`mkdir CARD`

`cd CARD`

`bunzip2 broadstreet-v1.2.1.tar.bz2 && tar -xvf broadstreet-v1.2.1.tar `  

Then make a DIAMOND database file form the protein homolog model. It is a fasta file with amino acid sequences.  

`diamond makedb --in protein_fasta_protein_homolog_model.fasta -d CARD`

The annotation of resistance genes will be done as a batch job using DIAMOND against CARD. Make the batch script and submit it as previously. (This takes less than an hour + possible time in the queue).  
```
#!/bin/bash -l
#SBATCH -J diamond
#SBATCH -o diamond_out_%j.txt
#SBATCH -e diamond_err_%j.txt
#SBATCH -t 05:00:00
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH -p serial
#SBATCH --mem-per-cpu=2000
#

module load biokit
cd $WRKDIR/BioInfo_course
diamond blastx -d ~/appl_taito/CARD/CARD -q trimmed_data/birds_R1.fasta \
            --max-target-seqs 1 -o birds_R1_CARD.txt -f 6 --id 90 --min-orf 20 \
            -p $SLURM_CPUS_PER_TASK --masking 0
```

When the job is finished, inspect the results.  
They are in one file with each hit on one line and we would like to have a count table where we have different ARGs as rows and our samples as columns.  
This is a stupid script, but does the job.  

`python scripts/parse_diamond/parse_diamond.py -i birds_R1_CARD.txt -o birds_CARD.csv`

## Assembly quality assesment
Let's take a look at the assembly file from yesterday. From the log file at `$WRKDIR/BioInfo_course/trimmed_data/co-assembly` you can check how the assembly run and at the last rows how is the output. However, for more detailed analysis we run [MetaQUAST](http://bioinf.spbau.ru/metaquast) together with the assembly. Copy folder called "assembly_QC" to your computer. We will view the results in your favorite browser. 

## Taxonomic profiling with Metaxa2 continued...
When all Metaxa2 array jobs are done, we can combine the results to an OTU table. Different levels correspond to different taxonomic levels.  
When using any 16S rRNA based software, be cautious with species (and beyond) level classifications. Especially when using short reads.  
We will look at genus level classification.
```
# Genus level taxonomy
metaxa2_dc -o birds_metaxa6.txt *level_6.txt
```
