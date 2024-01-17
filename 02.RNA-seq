#mapping
  STAR --runMode alignReads \
	--runThreadN 15 \
	--genomeDir susScr11 \
	--readFilesIn 1.clean.fq.gz 2.clean.fq.gz \
	--readFilesCommand zcat \
	--outSJfilterOverhangMin 30 16 16 16 \
	--outSJfilterCountUniqueMin 4 2 2 2 \
	--alignSJoverhangMin 6 \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical \
	--alignIntronMax 500000 \
	--alignMatesGapMax 500000 \
	--outFilterMismatchNoverReadLmax 0.07 \
	--outSAMtype BAM SortedByCoordinate \
	--waspOutputMode SAMtag /
	--outSAMattributes vA vG
	--outBAMsortingThreadN 15

#unique mapped reads
samtools view -q 255 -@ 10 -bS -O bam -o ${name}_sort.bam ${bamfile}
samtools index -@ 10 ${bamfile}

#RseQC
geneBody_coverage.py -i ${bamfile} -r 2020-3-30-UCSC-pig_annot.bed -f pdf -o ${name}


#gene expression
stringtie ${bamfile} -b ${OUTPUT} -e -G Sus_scrofa.gtf -C ${individual}_cov_ref.gtf -p 5 -o ${individual}
stringtie --merge -p 8 -G Sus_scrofa.gtf -o stringtie_merged.gtf merge.txt 
featureCounts -T 15 -F GTF -t exon -g gene_id -s 2 -Q 20 -C -p -D 1000 -a stringtie_merged.gtf -o ${individual}.fcounts ${bamfile}











