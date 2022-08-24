#!/bin/bash
#ls *.uniq.cut10.R1.fq.gz | while read id;do(STAR --runThreadN 20 --quantMode TranscriptomeSAM GeneCounts --limitGenomeGenerateRAM 90488309989 --genomeDir /home/j/backup1/refgenome/homo_sapiens/YSP_spikein/STAR_149/ --readFilesIn $id $(basename $id 'R1.fq.gz')R2.fq.gz  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outFilterMismatchNoverReadLmax 0.06  --outFileNamePrefix $(basename $id '.uniq.cut10.R1.fq.gz').STAR --outSAMmultNmax 1 --outSAMmapqUnique 60);done


#ls *.uniq.cut10.R1.fq.gz | while read id; do(hisat2 -p 20 --pen-noncansplice 12 --un-conc $(basename $id 'uniq.cut10.R1.fq.gz')unalign --al-conc $(basename $id 'uniq.cut10.R1.fq.gz')ambious -x /home/j/backup4/YSP-RT/ -1 $id -2 $(basename $id 'R1.fq.gz')R2.fq.gz -S $(basename $id '.uniq.cut10.R1.fq.gz').hisat2.sam); done

ls *.uniq.cut10.R1.fq.gz | while read id; do(hisat2 -p 20 --rna-strandness FR --pen-noncansplice 12 -x /home/j/backup4/YSP-RT/NN_m6A -1 $id -2 $(basename $id 'R1.fq.gz')R2.fq.gz -S $(basename $id '.uniq.cut10.R1.fq.gz').spikein_NN.sam); done

ls *.uniq.cut10.R1.fq.gz | while read id; do(hisat2 -p 20 --rna-strandness FR --pen-noncansplice 12 -x /home/j/backup4/YSP-RT/ref_format -1 $id -2 $(basename $id 'R1.fq.gz')R2.fq.gz -S $(basename $id '.uniq.cut10.R1.fq.gz').spikein.sam); done

