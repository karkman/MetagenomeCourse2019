#!/bin/bash
module load biokit
while read i
do
       name=($i)
       #Map fastq read pairs to Resfinder database and create bam files
        bowtie2 -x Resfinder -1 $WRKDIR/DONOTREMOVE/Metagenomics2019/fastqfiles/$name/$name"_1.fq" \
 -2  $WRKDIR/DONOTREMOVE/Metagenomics2019/fastqfiles/$name/$name"_2.fq" -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
       --threads 8 | samtools view -Sb - > $name.bam
       #Sort bam file and count reads which map in pairs and read pairs for which only one maps as one
       samtools view -h $name".bam" | awk '$7!="=" || ($7=="=" && and($2,0x40)) {print $0}' | samtools view -Su  - \
        | samtools sort -o $name"_sort.bam"
        #Index the bam file
        samtools index $name"_sort.bam"
        #Take the counts from column 3
        samtools idxstats $name"_sort.bam" | grep -v "*" | cut -f3 > $name"_counts"
        echo -e "GENE" > gene_names
        #Take the gene name from column 1
        samtools idxstats $name"_sort.bam" | grep -v "*" | cut -f1 >> gene_names
        #Add sample name to the counts
        echo -e "$name" > temp
        cat $name"_counts" >> temp
        mv -f temp $name"_counts"
        #Add sample name to progress file
        echo -e "$name" >> bowtie2_progress
done < names
       #Create the final ARG genematrix
        paste  gene_names *_counts > ARG_genemat.txt
