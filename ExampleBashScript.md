# Example to convert from fastq to fasta

Save the bash script in your `scripts` folder.  

```
#!/bin/bash

while read i
do 
    fastq_to_fasta_fast trimmed_data/$i"_R1_trimmed.fastq" > trimmed_data/$i"_R1_trimmed.fasta"
    awk '/^>/{print ">'$i'-" ++i; next}{print}' < trimmed_data/$i"_R1_trimmed.fasta" > trimmed_data/$i"_R1_renamed.fasta"
done < $1
```

Then run it from the `BioInfo_course`
```
bash scripts/NAME_OF_THE_SCRIPT.sh sample_names.txt
```

Now you should have trimmed and renamed fasta files in your `trimmed_data` folder.  

