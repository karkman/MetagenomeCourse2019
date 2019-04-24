#!/bin/bash
module load biokit
while read i
do
       name=($i)
        bowtie2 -x Resfinder -1  <(zcat $WRKDIR/DONOTREMOVE/Metagenomics2019/fastqfiles/$name/$name"_1.fq.gz") \
 -2 <(zcat $WRKDIR/DONOTREMOVE/Metagenomics2019/fastqfiles/$name/$name"_2.fq.gz") -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
       --threads 8 | samtools view -Sb - > $name.bam
       samtools view -h $name".bam" | awk '$7!="=" || ($7=="=" && and($2,0x40)) {print $0}' | samtools view -Su  - \
        | samtools sort -o $name"_sort.bam"
        samtools index $name"_sort.bam"
        samtools idxstats $name"_sort.bam" | grep -v "*" | cut -f3 > $name"_counts"
        echo -e "GENE" > gene_names
        samtools idxstats $name"_sort.bam" | grep -v "*" | cut -f1 >> gene_names
        echo -e "$name" > temp
        cat $name"_counts" >> temp
        mv -f temp $name"_counts"
        echo -e "$name" >> bowtie2_progress
done < names
        paste  gene_names *_counts > ARG_genemat.txt
