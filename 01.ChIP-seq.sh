#Mapping single-end reads with aln algorithm
REF=Sscrofa11.1.fa
LIBRARY="ChIP" 
LANE="L1"
RG="@RG\\tID:${NAME}_${LIBRARY}_${LANE}\\tLB:${NAME}_${LIBRARY}\\tPL:ILLUMINA\\tSM:${NAME}"

BWA aln -t 10 -q 15 ${REF} ${NAME}_ac_R1.fq.gz > ${NAME}.sai
BWA samse -r $RG ${REF} ${NAME}.sai ${NAME}_ac_R1.fq.gz > ${NAME}.sam 
Samtools view -q 1 -bSu ${NAME}.sam > ${NAME}_tmp.bam
Samtools sort ${NAME}_tmp.bam -O bam -o ${NAME}_sorted.bam  
Samtools index ${NAME}_sorted.bam  

#unique mapped reads
sambamba view -h -t 2 -f bam \
-F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${Bam_file} > ${name}_input_R1_sorted_unique.bam

#remove duplicate reads
${PICARD}/MarkDuplicates.jar INPUT=${Bam_file} OUTPUT=${name}_dupmark.bam METRICS_FILE=${name}_dup.qc \
VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
samtools view -F 1804 -b ${name}_dupmark.bam -o ${name}_final.bam
samtools sort ${name}_final.bam -O bam -o ${name}_final_sort.bam
samtools index ${name}_final_sort.bam

#call peak
MACS filterdup -i ${Sample}_ac_final_sort.bam -f BAM -g hs --keep-dup=1 --verbose=3 -o ${name}_filterdup.bed
MACS predictd -i ${name}_filterdup.bed -g 2.50e9 > ${name}_temp.txt 2>&1
MACS callpeak -t ${Sample}_ac_final_sort.bam -c ${Sample}_input_final_sort.bam \
-n ${Sample}_ac -g 2.50e9  -p 1e-2 -f BAM --nomodel --extsize ${name}_temp --keep-dup all -B --SPMR

#peak activity
bedtools multicov -bams ${Sample}_ac_final_sort.bam -bed consensus_peak.txt > name_depth.bed

