#reference
snp2h5 --chrom chromInfo \
         --format vcf \
         --haplotype haplotypes.h5 \
         --snp_index snp_index.h5 \
         --snp_tab   snp_tab.h5 \
        all_snp1.vcf.gz


#Flagging mapped reads with WASP
python3 find_intersecting_snps.py \
--output_dir out_path --snp_index snp_index.h5 \
--snp_tab snp_tab.h5 --haplotype haplotypes.h5 \
--is_sorted ac_final_sort.bam

#remap
NAME=Fastq
LIBRARY="ChIP" 
LANE="L1"
RG="@RG\\tID:${NAME}_${LIBRARY}_${LANE}\\tLB:${NAME}_${LIBRARY}\\tPL:ILLUMINA\\tSM:${NAME}"
BWA aln -t 10 -q 15 Sscrofa11.1.fa ac_final_sort.remap.fq.gz > name.remap.sai
BWA samse -r $RG Sscrofa11.1.fa name.remap.sai name_ac_final_sort.remap.fq.gz > name.remap.sam 
Samtools view -q 1 -bSu name.remap.sam > name_tmp.remap.bam
Samtools sort name_tmp.remap.bam -O bam -o name_sorted.remap.bam  
Samtools index name_sorted.remap.bam  

#Filtering reads with WASP
python3 filter_remapped_reads.py ac_final_sort.to.remap.bam \
name_sorted.remap.bam \
name_sorted.remap.filtered.bam

#merge and sort
samtools merge name.merged.bam ac_final_sort.keep.bam name_sorted.remap.filtered.bam
samtools sort name.merged.bam -O bam -o name_sorted.merged.bam
samtools index name_sorted.merged.bam











