#!/bin/bash
ls ../../../../*.R2.coord.bam | while read id;do(bedtools genomecov -ibam $id -bg -scale 1 -strand + -5 > $(basename $id '.R2.coord.bam').R2.for.5.bedgraph &);done
#ls ../../../../*.R2.coord.bam | while read id;do(bedtools genomecov -ibam $id -bg -scale 1 -strand - -5 > $(basename $id '.R2.coord.bam').R2.rev.5.bedgraph &);done
#samtools在计算depth的过程中对reads进行了过滤，那么它过滤的原则是什么呢？
#可以看到，其中有一个SECONDARY比对，对应数字355，这样的reads有4条，去除这4条之后，就是samtools计算的最终结果了。
#作为SAM文件的官方操作工具，samtools内置了对flag的过滤， 而bedtools默认没有进行这样的过滤，只是简单统计了该区域比对上的reads数目。

samtools depth -aa -b IP.R2.for.5.expand.bedgraph.extend50.bed -f depth.for.txt > IP.R2.for.5.expand.bedgraph.extend50.depth.bed &
samtools depth -aa -b IP.R2.rev.5.expand.bedgraph.extend50.bed -f depth.rev.txt > IP.R2.rev.5.expand.bedgraph.extend50.depth.bed &




