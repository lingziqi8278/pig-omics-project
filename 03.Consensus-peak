#generating consensus peak

library(DiffBind)
samples <- data.frame(SampleID=name,
Tissue="liver",Factor=covariate[name,"sex"],Condition=covariate[name,"sex"],
Treatment=covariate[name,"sex"],Replicate="1",ControlID="input",
bamReads=name_ac_final_sort.bam,
bamControl=name_input_final_sort.bam,
Peaks=name_ac_peaks.narrowPeak,
PeakCaller="macs")
dbObj <- dba(sampleSheet=samples,minOverlap=3)
consensus_peaks <- dba.peakset(dbObj, bRetrieve=TRUE,)
consensus_peak <- as.data.frame(consensus_peaks)
