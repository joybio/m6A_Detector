#!/bin/python
#hisat2: 3'ligation: FR; dUTP: RF
#fr-secondstrand。意思就是read1在+链，相对的gene也同样在+链上，而read2在+链，相对的gene在-链上。
import os

os.system("ls *_1.fq.gz | while read id;do(cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 10 -e 0.1 -O 5 -m 20 -q 30,30 -o $(basename $id '_1.fq.gz').trimmed.R1.fq.gz -p $(basename $id '_1.fq.gz').trimmed.R2.fq.gz $id $(basename $id '_1.fq.gz')_2.fq.gz);done")

#remove PCR duplication
os.system("ls *.trimmed.R1.fq.gz | while read id;do(seqkit rmdup -s -i $id -o $(basename $id 'trimmed.R1.fq.gz')uniq.R1.fq.gz);done")
os.system("ls *uniq.R1.fq.gz | while read id;do(seqkit seq --name --only-id $id > $(basename $id 'uniq.R1.fq.gz')id.txt);done")

os.system("ls *.trimmed.R2.fq.gz | while read id;do(seqkit grep --pattern-file $(basename $id 'trimmed.R2.fq.gz')id.txt $id -o $(basename $id 'trimmed.R2.fq.gz')uniq.R2.fq.gz);done")
os.system("ls *uniq.R1.fq.gz | while read id;do(cutadapt -u 3 $id -o $(basename $id '.R1.fq.gz').cut3.R1.fq.gz);done")
os.system("ls *uniq.R2.fq.gz | while read id;do(cutadapt -u -3 $id -o $(basename $id '.R2.fq.gz').cut3.R2.fq.gz);done")

os.system("mkdir sfastqc_results")
#QC
os.system("ls *.cut3.R1.fq.gz | while read id; do(mkdir -p sfastqc_results/$(basename $id '.fq.gz');fastqc -o sfastqc_results/$(basename $id '.fq.gz') $id &);done")
os.system("ls *.cut3.R2.fq.gz | while read id; do(mkdir -p sfastqc_results/$(basename $id '.fq.gz');fastqc -o sfastqc_results/$(basename $id '.fq.gz') $id);done")
#multiQC
os.system("multiqc sfastqc_results/ -o multiqc_results")
#map 
os.system("ls *.uniq.cut3.R1.fq.gz | while read id;do(STAR --runThreadN 20 --quantMode TranscriptomeSAM GeneCounts --limitGenomeGenerateRAM 90488309989 --genomeDir /home/j/backup1/refgenome/mm/GRCm39/QY_spikein/STAR_149 --readFilesIn $id $(basename $id 'R1.fq.gz')R2.fq.gz  --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outFilterMismatchNoverReadLmax 0.06  --outFileNamePrefix $(basename $id '.uniq.cut3.R1.fq.gz').STAR --outSAMmultNmax 1 --outSAMmapqUnique 60);done")







