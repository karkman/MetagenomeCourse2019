# Metagenome analysis of infant gut metagenomes - part 3


Anvio is an analysis and visualization platform for omics data. You can read more from Anvio's [webpage](http://merenlab.org/software/anvio/).

![alt text](https://github.com/INNUENDOCON/MicrobialGenomeMetagenomeCourse/raw/master/Screen%20Shot%202017-12-07%20at%2013.50.20.png "Tom's fault")

Go to your course folder and make a new folder called ANVIO. All task on this section are to be done in this folder.

```
mkdir ANVIO
cd ANVIO
```
We need to do some tricks for the contigs from assembly before we can use them on Friday. You'll need to load bioconda and activate Anvio.  
Open a screen for Anvi'o workflow. `screen -S anvio`

```
# allocate resources and log in to a computing node
salloc -n 1 --mem=10000 -t 06:00:00 -p serial
srun --pty $SHELL
# load bioconda and activate Anvi'o environment
module load bioconda/3
source activate anvio5
```
## Rename the scaffolds and select those >2,500nt.
Anvio wants sequence IDs in your FASTA file as simple as possible. Therefore we need to reformat the headerlines to remove spaces and non-numeric characters. Also contigs shorter than 2500 bp will be removed.


```
anvi-script-reformat-fasta ../co-assembly/final.contigs.fa -l 2500 --simplify-names --prefix MEGAHIT_co_assembly -r REPORT -o MEGAHIT_co-assembly_2500nt.fa
````
Deattach from the screen with `Ctrl a+d`  

## Mapping the reads back to the assembly
Next thing to do is mapping all the reads back to the assembly. We use the renamed >2,500 nt contigs and do it sample-wise, so each sample is mapped separately using the trimmed R1 & R2 reads.  
We will need to two scripts for that, one for the actual mapping and another to run it as an array job. Save both scripts to your `scripts` folder.

But before doing that we have make a bowtie2 index from the contig file. Run the following command:  

```
bowtie2-build MEGAHIT_co-assembly_2500nt.fa co-assembly  
```

`co-assembly` is the base name for the resulting index files.  


The mapping script from Tom, modified to be used in Taito:
```
#!/bin/bash
# example run: ./bowtie2-map-batch.sh SAMPLE1_R1.fastq SAMPLE1_R2.fastq SAMPLE1 bt2_index

# $1: Forward reads
# $2: Reverse reads
# $3: Sample name
# $4: Bowtie2 index

set -e

bowtie2 --threads $SLURM_CPUS_PER_TASK -x $4 -1 $1 -2 $2 -S $3.sam --no-unal
samtools view -F 4 -bS $3.sam > $3-RAW.bam
samtools sort $3-RAW.bam -o $3.bam
samtools index $3.bam
rm $3.sam $3-RAW.bam
```
The array job script:
```
#!/bin/bash -l
#SBATCH -J array_map
#SBATCH -o array_map_out_%A_%a.txt
#SBATCH -e array_map_err_%A_%a.txt
#SBATCH -t 00:30:00
#SBATCH --mem-per-cpu=1000
#SBATCH --array=1-6
#SBATCH -n 1
#SBATCH --cpus-per-task=6
#SBATCH -p serial

cd $WRKDIR/Metagenomics2019/co-assembly
# we need Bowtie2 from the biokit
module load biokit
# each job will get one sample from the sample names file
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)
# run mapping script for each sample
bash ../scripts/bowtie2-map-batch.sh ../trimmed_data/$name"_R1_trimmed.fastq" \
     ../trimmed_data/$name"_R2_trimmed.fastq" \
     $name ../ANVIO/co-assembly
```
Then again submit the array job with `sbatch`.  

During luch break check what happens in the different steps in the mapping script.
```
bowtie2
samtools view
samtools sort
samtools index
```

## Back to Anvi'o

Reattach to your Anvi'o screen

```
screen -r anvio
```


## Generate CONTIGS.db

Contigs database (contigs.db) contains information on contig length, open reading frames (searched with Prodigal) and kmers. See [Anvio webpage](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#creating-an-anvio-contigs-database) for more information.  

```
anvi-gen-contigs-database -f MEGAHIT_co-assembly_2500nt.fa -o MEGAHIT_co-assembly_2500nt_CONTIGS.db -n MEGAHIT_co-assembly
```
## Run HMMs to identify single copy core genes for Bacteria and Archaea, plus rRNAs
```
anvi-run-hmms -c MEGAHIT_co-assembly_2500nt_CONTIGS.db -T 6
```

While the HMM identification is running, deattach from the screen and check if the mapping has been done.  
```
squeue -l -u $USER
```

When the mapping is done for all samples and the contigs database is ready, we can profile the samples using the DB and the mapping output. Write an array script for the profiling and submit it to the queue.

```
#!/bin/bash -l
#SBATCH -J array_profiling
#SBATCH -o array_profiling_out_%A_%a.txt
#SBATCH -e array_profiling_err_%A_%a.txt
#SBATCH -t 00:30:00
#SBATCH --mem-per-cpu=1000
#SBATCH --array=1-6
#SBATCH -n 1
#SBATCH --cpus-per-task=6
#SBATCH -p serial

cd $WRKDIR/Metagenomics2019/co-assembly
# we need to load Bioconda
module load bioconda/3
# then activate Anvi'o
source activate anvio5
# and also load the biokit, since Anvi'o uses samtools in the profiling
module load biokit
# each job will get one sample from the sample names file
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)
# then the actual profiling
anvi-profile -c ../ANVIO/MEGAHIT_co-assembly_2500nt_CONTIGS.db  -M 2500 -T $SLURM_CPUS_PER_TASK -i $name.bam -o $name
```
Submit the job with `sbatch` as previously.  

## Run COGs DOES NOT WORK; DO NOT SPEND TIME ON THIS

Next we annotate genes in  contigs database with functions from the NCBI’s Clusters of Orthologus Groups (COGs).
Again first reattach to your Anvio'o screen.  

```
anvi-run-ncbi-cogs -c MEGAHIT_co-assembly_2500nt_CONTIGS.db -T 6
```

## Export GENES
With this command we export the genecalls from Prodigal to gene-calls.fa and do taxonomic annotation against centrifuge database you installed on Wednesday

```
anvi-get-dna-sequences-for-gene-calls -o gene-calls.fa -c MEGAHIT_co-assembly_2500nt_CONTIGS.db
```

## Run centrifuge
“classification engine that enables rapid, accurate and sensitive labeling of reads and quantification of species on desktop computers”. Read more from [here](http://biorxiv.org/content/early/2016/05/25/054965).

Remember to set the environmental variable pointing to the centrifuge folder as shown in [MetagenomeInstallations](https://github.com/INNUENDOCON/MicrobialGenomeMetagenomeCourse/blob/master/MetagenomeInstallations.md).
# LINKKI MENEE MIRKON GITHUBBIIN

```
centrifuge -f -p 6 -x $CENTRIFUGE_BASE/p+h+v gene-calls.fa -S centrifuge_hits.tsv
```
## Import centrifuge results
```
anvi-import-taxonomy -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge -c MEGAHIT_co-assembly_2500nt_CONTIGS.db
```

## Antibiotic resistance gene annotation
```
# run BLASTN against ResFinder
blastn -query gene-calls.fa -subject ~/appl_taito/database/resfinder_db/resfinder.fasta -out RF.out -outfmt 6 -max_target_seqs 1 -perc_identity 90 -qcov_hsp_perc 80
# parse the output for Anvi'o  
printf "gene_callers_id\tsource\taccession\tfunction\te_value\n" > RF_functions.txt
awk '{print $1"\tCARD\t\t"$2"\t"$11}' RF.out >> RF_functions.txt
# and finally import the functions to Anvi'o contigs DB  
anvi-import-functions -c MEGAHIT_co-assembly_2500nt_CONTIGS.db -i RF_functions.txt
```

## Merging the profiles
When the profiling is done, you can merge them with one command.

```
anvi-merge ../co-assembly/*/PROFILE.db -o SAMPLES-MERGED -c MEGAHIT_co-assembly_2500nt_CONTIGS.db
```

## Visualization in the interface 
You don't need to specify any port when running Anvi'o on your own laptop.  
But when running the interactive interface from Taito, you will need your own port, because it is not possible to run two interfaces thru the same port.  
The available ports will assigned to each student on the course.

Remember to change the `XXXX` to the port you were given.  

Open a new ssh window. In mac:
```
ssh -L XXXX:localhost:XXXX YOUR_USERNAME@taito.csc.fi
```

in Windows with Putty:
In SSH tab select "tunnels". Add  

Source port: XXXX  
Destination: localhost:XXXX  

Click add and log in to Taito.

Activate anvio

```
anvi-interactive -c MEGAHIT_co-assembly_2500nt_CONTIGS.db -p SAMPLES-MERGED/PROFILE.db --server-only -P XXXX
```

Then open google chrome and go to address 

http://localhost:8080

**but change 8080 to your port number
