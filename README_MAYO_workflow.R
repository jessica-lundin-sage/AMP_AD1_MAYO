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
output: MAYO_md_counts_cqn.rds ("metadata", "counts", "dge_cqn", "dge_cqn_df") (synxxx)  

On AWS:
  DE_QC_MAYO_xx.Rmd  #here for reference but not longer used
     Outlier detection using PCA (crude)
     Calculate SVs - no longer using bc of PCA of RNA metrics
     Variance partitioning visualization - figure saved
     Correlation of model PCs with covariates

On AWS:
  DE_with_built_interaction_no_svs.R 
     output "MAYO_DE_final.rds" (synxx) contains:
       "metadata" = md_sv, 
       "vobj_expr" = voom_gene_expression, 
       "dge_cqn" = dge_cqn, 
       "ebayes" = dream.cont.ebayes6, 
       "fit_contrasts" = fit_contrasts6,
       "males_CER6" = males_CER6, 
       "females_CER6" = females_CER6,  
       "males_TCX6" = males_TCX6, 
       "females_TCX6" = females_TCX6
Brings in ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds ("metadata", "counts", "dge_cqn", "dge_cqn_df") (syn76553509)  


On AWS:
  DE_Residuals_for_sharing.R
     Models technical variables only on cqn normalized counts using dream weights. formula: ~ RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|flowcell)
     "MAYO_DE_res2.rds" (synxx) contains: 
     "metadata_res3" = md_sv, 
     "fit_res2.dream" = fit_res2
     "fit_res3.ebayes" = fit_res3
     "males_CER6_res3" = males_CER6, 
     "females_CER6_res3" = females_CER6,  
     "males_TCX6_res3" = males_TCX6, 
     "females_TCX6_res3" = females_TCX6
 
     ## may also save to MAYO_DE_res.rds    
     "vobj_res" = voom_res (voomwithDreamWeights output), 
     "dge_cqn" = dge_cqn (counts with dge_can$E from CQN normalization), 
     "fit_res" = fit_res (dream() with counts, CQN offset (E), and dream weights), 
     "residual_gene_expression" = residual_gene_expression (residuals from dream())


form_ck <- ~ 0 + group  + apoe4Status + age_cat + (1|individualID)
How group was defined: md_sv$group <- interaction(md_sv$diagnosis, md_sv$sex, md_sv$tissue2, drop = TRUE)
group <- interaction(md_sv$diagnosis, md_sv$sex, md_sv$tissue2, drop = TRUE)
"ROSMAP_DE_res2.rds" (syn76565771) contains: 
  "metadata_res3" = md_sv, 
"fit_res2.dream" = fit_res2 (dream() with residualized counts, formula with additional vars (form_ck), and tissue X sex X diagnosis contrasts), 
"fit_res3.ebayes" = fit_res3 (eBayes of fit_res2), 
"males_DLPFC6_res3" = males_DLPFC6 (topTable of fit_res3 for each contrast), 
"females_DLPFC6_res3" = females_DLPFC6,  
"males_PCC6_res3" = males_PCC6, 
"females_PCC6_res3" = females_PCC6, 
"males_CN6_res3" = males_CN6, 
"females_CN6_res3" = females_CN6

On local (moved to be part of QC Markdown):
  DE_residualized_plots_by_final_batch.R
Run PCA on residualized data with overlay markers for sex and diagnosis - plots output in Rmd



Technical variables
technical_stats_multiqc_star.R (syn76227881) Output: MSBB_multiqc_star_technical_stats.csv
technical_stats_fastqc.R (syn76283403) [run on AWS] (ROSMAP_fq_stats.rds contains: basic_stats.txt, phred_per_base.txt, base_content.txt) 




Helpful links:
  https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
https://diseaseneurogenomics.github.io/variancePartition/articles/dream.html
https://github.com/DiseaseNeuroGenomics/variancePartition/blob/HEAD/R/dream.R
https://github.com/Sage-Bionetworks/amp-rnaseq
https://github.com/Sage-Bionetworks/ampad-rnaseq-reprocessing/tree/main/code/metadata_preprocessing
https://github.com/Sage-Bionetworks/sageseqr/blob/master/R/functions.R
https://github.com/Sage-Bionetworks/ampad-DiffExp/blob/master/gene_level_analysis/MAYO_geneLevel_TMM.Rmd
https://github.com/Sage-Bionetworks/sageseqr/blob/master/R/functions.R




