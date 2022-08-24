#!/bin/bash
ls *.bedgraph | while read id;do(python expand.py -i $id -o $(basename $id '.bedgraph')expand.bedgraph &);done
python compute_stop_coverage_by_R2_reads_start.py -f YSP_IP_1_9.R2.for.5expand.bedgraph -r YSP_IP_1_9.R2.rev.5expand.bedgraph -o  YSP_IP_1_9.Readsend.filter.xls -n 3 -s 2 
python compute_stop_coverage_by_R2_reads_start.py -f YSP_IP_2_10.R2.for.5expand.bedgraph -r YSP_IP_2_10.R2.rev.5expand.bedgraph -o YSP_IP_2_10.Readsend.filter.xls -n 3 -s 2
awk 'BEGIN{OFS="\t"}$2=$2-3' YSP_IP_1_9.R2.for.5expand.bedgraph > YSP_IP_1_9.R2.for.5expand.bedgraph.extend3.step1.bed
awk 'BEGIN{OFS="\t"}$3=$3+3' YSP_IP_1_9.R2.for.5expand.bedgraph.extend3.step1.bed >  YSP_IP_1_9.R2.for.5expand.bedgraph.extend3.bed

awk 'BEGIN{OFS="\t"}$2=$2-3' YSP_IP_1_9.R2.rev.5expand.bedgraph > YSP_IP_1_9.R2.rev.5expand.bedgraph.extend3.step1.bed
awk 'BEGIN{OFS="\t"}$3=$3+3' YSP_IP_1_9.R2.rev.5expand.bedgraph.extend3.step1.bed >  YSP_IP_1_9.R2.rev.5expand.bedgraph.extend3.bed

awk 'BEGIN{OFS="\t"}$2=$2-3' YSP_IP_2_10.R2.for.5expand.bedgraph > YSP_IP_2_10.R2.for.5expand.bedgraph.extend3.step1.bed
awk 'BEGIN{OFS="\t"}$3=$3+3' YSP_IP_2_10.R2.for.5expand.bedgraph.extend3.step1.bed >  YSP_IP_2_10.R2.for.5expand.bedgraph.extend3.bed

awk 'BEGIN{OFS="\t"}$2=$2-3' YSP_IP_2_10.R2.rev.5expand.bedgraph > YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.step1.bed
awk 'BEGIN{OFS="\t"}$3=$3+3' YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.step1.bed >  YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.bed

rm *.step1.bed
samtools depth -aa -b YSP_IP_1_9.R2.for.5expand.bedgraph.extend3.bed -f depth.for.txt > YSP_IP_1_9.R2.for.5expand.bedgraph.extend3.depth &
samtools depth -aa -b YSP_IP_2_10.R2.for.5expand.bedgraph.extend3.bed -f depth.for.txt > YSP_IP_2_10.R2.for.5expand.bedgraph.extend3.depth &
samtools depth -aa -b YSP_IP_1_9.R2.rev.5expand.bedgraph.extend3.bed -f depth.rev.txt > YSP_IP_1_9.R2.rev.5expand.bedgraph.extend3.depth &
awk '{if($2<0){print $1"\t"0"\t"$3"\t"$4"\t"$5}}{if($2>0){print}}' YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.bed > YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.bed1
mv YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.bed1 YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.bed
samtools depth -aa -b YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.bed -f depth.rev.txt > YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.depth &

python compute_arrest_rate_of_filtered_stop.py -i YSP_IP_1_9.Readsend.filter.xls -f YSP_IP_1_9.R2.for.5expand.bedgraph.extend3.depth -r YSP_IP_1_9.R2.rev.5expand.bedgraph.extend3.depth -o YSP_IP_1_9.Readsend.filter.Arrest.xls -l YSP_IP_1_9.Readsend.filter.Arrest.filter.xls -n 4 -a 0.2 &
python compute_arrest_rate_of_filtered_stop.py -i YSP_IP_2_10.Readsend.filter.xls -f YSP_IP_2_10.R2.for.5expand.bedgraph.extend3.depth -r YSP_IP_2_10.R2.rev.5expand.bedgraph.extend3.depth -o YSP_IP_2_10.Readsend.filter.Arrest.xls -l YSP_IP_2_10.Readsend.filter.Arrest.filter.xls -n 5 -a 0.2 &

