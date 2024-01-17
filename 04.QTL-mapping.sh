#cis-eQTL/cis-acQTL mapping
python3.6 acqtl_prepare_expression.py FPM.gct counts.gct Sus.sorted.gtf \
    sample.txt all_chr_list.txt ac-TMM \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm

Rscript run_PEER.R -m 100 -o ac-TMM.bed.gz ac 45

python3.6 run_FastQTL_threaded.py all_vcf.gz \
ac_adjust_covariants.TMM.bed.gz fastq-permu-ac \
--window 1e6 --permute 1000 10000 --maf_threshold 0.01 \
--ma_sample_threshold 10 --chunks 100 --threads 16

python3.6 run_FastQTL_threaded.py all_vcf.gz \
ac_adjust_covariants.TMM.bed.gz fastq-nominal-ac \
--window 1e6 --ma_sample_threshold 10 \
--maf_threshold 0.01 --chunks 100 --threads 16


#trans-eQTL/trans-acQTL mapping
QTLtools trans --vcf all_vcf.gz --bed ac_adjust_covariants.TMM.bed.gz --window 1e6 \
--sample 1000  --out trans-peak-permutations.txt

QTLtools trans --vcf all_vcf.gz --bed ac_adjust_covariants.TMM.bed.gz \
--adjust trans-peak-permutations.txt.best.txt.gz --threshold 0.1 --out trans-peak.adjust 

Rscript runFDR_atrans.R trans-peak.adjust.best.txt.gz \
trans-peak.adjust.hits.txt.gz 0.05 peak-trans-final.txt



