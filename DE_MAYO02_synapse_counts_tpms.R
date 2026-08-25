## Differential expression analyses on RNAseq data from MAYO
# STEP 02: count and tpm data from synapse
# J Lundin
# June 16 2026


pacman::p_load(tidyverse, limma, edgeR, biomaRt, DESeq2, vsn, sva, pamr)
pacman::p_load(synapser,dplyr,purrr,readr,lubridate,stringr,tibble,ggplot2)
synLogin()

work_dir <- ("C:/Users/jlundin/OneDrive - Sage Bionetworks/RNASeq_Harm/MAYO")
setwd(work_dir)

source("C:/Users/jlundin/OneDrive/git_code/RNASeq_DE/AMP_AD1/RNASeq_DE/AMP_AD1/functions/functions_filter_low_count_genes.R")

### pulling count and tpm data from synapse ----
mayo_counts_tc_temp <- synapser::synGet("syn69369067") 
  mayo_counts_tc <- read.csv(mayo_counts_tc_temp$path, sep="\t", header=T, check.names = FALSE)
mayo_tpm_tc_temp <- synapser::synGet("syn69369068") 
  mayo_tpm_tc <- read.csv(mayo_tpm_tc_temp$path, sep="\t", header=T, check.names = FALSE)

mayo_counts_cb_temp <- synapser::synGet("syn69076198") 
  mayo_counts_cb <- read.csv(mayo_counts_cb_temp$path, sep="\t", header=T, check.names = FALSE)
mayo_tpm_cb_temp <- synapser::synGet("syn69076199") 
  mayo_tpm_cb <- read.csv(mayo_tpm_cb_temp$path, sep="\t", header=T, check.names = FALSE)

mayo_md_temp <- synapser::synGet("syn76141636") #mayo_md_all
  mayo_md <- read.csv(mayo_md_temp$path, header=T)

## filtering on low counts ----
filtered_genes_mayo_tc <- filter_gene_expression(
  tpm_file   = mayo_tpm_tc,
  reads_file  = mayo_counts_tc,
  metadata  = mayo_md,
  synid_outfile = c("syn75817798"),
  study_var = c("mayo_tcx_all")
  )

filtered_genes_mayo_cb <- filter_gene_expression(
  tpm_file   = mayo_tpm_cb,
  reads_file  = mayo_counts_cb,
  metadata  = mayo_md,
  synid_outfile = c("syn75817798"),
  study_var = c("mayo_cer_all")
)
