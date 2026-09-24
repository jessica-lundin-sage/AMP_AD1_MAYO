## Preprocessing RNASeq data for AMP-AD1.0
# J Lundin
# July 14 2026

## WORKFLOW FOR MAYO


DE_MAYO01_synapse_metadata.R
   Loads and cleans metadata
   output: final md_all file (synxx) (n=xx)

DE_MAYO02_synapse_counts_tpms_from_synapse.R
   Loads and filters counts and tpms 
   counts and tpms filtered for both counts and tpms thresholds and were output as separate files 
   output: MAYO_all_counts_filtered.txt (synxx)
   output: MAYO_all_tpm_filtered.txt (synxx)

   Check for additional count and tpm filter stratified by tissue and sex
   output: list of gene_ids to remove based on tissue and sex specific filtering (synxx)

Technical variables
   technical_stats_multiqc_star.R Output: MAYO_multiqc_star_technical_stats.csv
   technical_stats_fastqc.R [run on AWS] (MAYO_fq_stats.rds contains: basic_stats.txt, phred_per_base.txt, base_content.txt) 

MAYO_QC.Rmd
    Run checks on RIN thresholds, inferred sex, fast-qc, multiqc, calcuate technical covariate for RNA metrics (plus technical outliers check) -- filter final md and count files
    output: MAYO_md_counts.rds ("metadata", "counts") (synxxx)

    Ran CQN normalization on QC file
    output: MAYO_md_counts_cqn_final.rds ("metadata", "counts", "dge_cqn", "dge_cqn_df") (synxxx)  

#here for reference but not longer used
DE_QC_MAYO_xx.Rmd  
     Outlier detection using PCA (crude)
     Calculate SVs - no longer using bc of PCA of RNA metrics
     Variance partitioning visualization - figure saved
     Correlation of model PCs with covariates

