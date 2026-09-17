## Reprossing RNASeq data for AMP-AD1.0
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

MAYO_QC.Rmd
    Run checks on RIN thresholds, inferred sex, fast-qc, multiqc, calcuate technical covariate for RNA metrics (plus technical outliers check) -- filter final md and count files
    output: MAYO_md_counts.rds ("metadata", "counts") (synxxx)

    Ran CQN normalization on QC file
    output: MAYO_md_counts_cqn_final.rds ("metadata", "counts", "dge_cqn", "dge_cqn_df") (synxxx)  

On AWS:
  DE_QC_MAYO_xx.Rmd  #here for reference but not longer used
     Outlier detection using PCA (crude)
     Calculate SVs - no longer using bc of PCA of RNA metrics
     Variance partitioning visualization - figure saved
     Correlation of model PCs with covariates

On AWS:
  MAYO_DE_FINAL_models.R  Models cqn normalized counts using dream weights.
    input: MAYO_md_counts_cqn_FINAL.rds ("metadata", "counts", "dge_cqn" and "dge_cqn_df")  
    output MAYO_DE_final.rds (synid below) contains:
       "metadata" = md_sv, 
       "vobj_expr" = voom_gene_expression, 
       "dge_cqn" = dge_cqn, 
       "ebayes" = dream.cont.ebayes6, 
       "fit_contrasts" = fit_contrasts6,
       "males_CER6" = males_CER6, 
       "females_CER6" = females_CER6,  
       "males_TCX6" = males_TCX6, 
       "females_TCX6" = females_TCX6



On AWS:
  DE_Residuals_for_sharing.R   Models technical variables only on cqn normalized counts using dream weights. formula: ~ RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|flowcell)
    input: MAYO_md_counts_cqn_FINAL.rds ("metadata", "counts", "dge_cqn" and "dge_cqn_df")  
    output: MAYO_DE_res.rds (synid below) contains:   
           "metadata" = md_sv
           "vobj_res" = voom_res (voomwithDreamWeights output), 
           "dge_cqn" = dge_cqn (counts with dge_can$E from CQN normalization), 
           "fit_res" = fit_res (dream() with counts, CQN offset (E), and dream weights), 
           "residual_gene_expression" = residual_gene_expression (residuals from dream())
      
    output: MAYO_DE_res2.rds (synid below) contains: 
           "metadata_res3" = md_sv, 
           "fit_res2.dream" = fit_res2 (dream() with residualized counts, formula with additional vars (form_ck), and tissue X sex X diagnosis contrasts), 
           "fit_res3.ebayes" = fit_res3 (eBayes of fit_res2), 
           "males_CER6_res3" = males_CER6, (topTable of fit_res3 for each contrast), 
           "females_CER6_res3" = females_CER6,  
           "males_TCX6_res3" = males_TCX6, 
           "females_TCX6_res3" = females_TCX6



Technical variables
technical_stats_multiqc_star.R Output: MAYO_multiqc_star_technical_stats.csv
technical_stats_fastqc.R [run on AWS] (MAYO_fq_stats.rds contains: basic_stats.txt, phred_per_base.txt, base_content.txt) 
