
## Differential expression analyses on RNAseq data from MAYO
# STEP 01: metadata from synapse
# J Lundin
# June 16 2026


pacman::p_load(tidyverse, limma, edgeR, biomaRt, DESeq2, vsn, sva, pamr)
pacman::p_load(synapser,dplyr,purrr,readr,lubridate,stringr,tibble,ggplot2)
synLogin()

work_dir <- ("C:/Users/jlundin/OneDrive - Sage Bionetworks/RNASeq_Harm/MAYO")
setwd(work_dir)

### pulling metadata from synapse ----
MAYO_meta_ind <- read.csv(synapser::synGet('syn73713766')$path, stringsAsFactors = F, check.names = FALSE) #from metadata harmonization study
MAYO_meta_biosp <- read.csv(synapser::synGet('syn20827192')$path, stringsAsFactors = F, check.names = FALSE)
MAYO_meta_assay <- read.csv(synapser::synGet('syn20827193')$path, stringsAsFactors = F, check.names = FALSE)

MAYO_meta_biosp <- MAYO_meta_biosp %>% filter(assay == "rnaSeq"& tissue != 'blood')
#length(unique(MAYO_meta_biosp$individualID)) #620

table(MAYO_meta_biosp$tissue)
table(MAYO_meta_assay$specimenID %in% MAYO_meta_biosp$specimenID)

metadata_temp <- merge(MAYO_meta_assay, MAYO_meta_biosp, by="specimenID")
metadata_temp2 <- merge(metadata_temp, MAYO_meta_ind, by = "individualID")
comb <- metadata_temp2

table(comb$libraryPrep)
table(comb$libraryPreparationMethod)
table(comb$tissue, comb$sex, comb$exclude)
table(comb$exclude)
table(comb$exclude,comb$diagnosis)
table(comb$exclude,comb$tissue)
table(comb$excludeReason)
table(comb$diagnosis)
summary(comb$RIN)
summary(comb$PMI)

#comb2 <- comb %>% filter(exclude == FALSE & (diagnosis == "Alzheimer Disease" | diagnosis == "control"))
comb2 <- comb %>% filter(exclude == FALSE)

table(comb2$diagnosis, comb2$sex, comb2$tissue)
table(comb2$diagnosis, comb2$tissue)
table(comb2$tissue)
summary(comb2$RIN)
summary(comb2$PMI)
table(comb2$apoe4Status)

#apoe4, is now #apoeGenotype (22, 23, 24, 33, 34, 44) and #apoe4Status (yes (if any 4) or no) with clinical_harmonized dataset
table(comb$apoeGenotype, comb$apoe4Status)

comb2 <- comb2 %>% mutate(tissue2 = recode(tissue, 
                                             "temporal cortex" = "TCX",
                                             "cerebellum"  = "CER"))  

comb2$age_cat <- NA
comb2$age_cat[comb2$ageDeath == "90+"] <- "90+"
comb2$age_cat[as.numeric(comb2$ageDeath[comb2$ageDeath != "90+"])<85] <- "<85"
comb2$age_cat[as.numeric(comb2$ageDeath[comb2$ageDeath != "90+"])>=85] <- "ge85lt90"
comb2$age_cat[comb2$ageDeath == "90+"] <- "90+"


# Sequencing Statistics - UPDATED 7-10-2026
metrics <- read.csv(synapser::synGet('syn76816214')$path, stringsAsFactors = F)

metrics2 <- metrics %>% dplyr::select(specimenID, "picard_UNPAIRED_READS_EXAMINED"  ,          
                                      "picard_READ_PAIRS_EXAMINED" ,                  
                                      "picard_UNMAPPED_READS"  ,                         
                                      "picard_PERCENT_DUPLICATION"   ,                   
                                      "picard_READS_UNMAPPED"     ,                      
                                      "rsem_alignable_percent" ,                         
                                      "rsem_uniquely_aligned_percent"  ,   
                                      "samtools_reads_mapped", 
                                      "samtools_reads_duplicated",
                                      "samtools_error_rate",                             
                                      "samtools_average_length"  ,                
                                      "samtools_average_quality"  ,                       
                                      "samtools_insert_size_average" ,       
                                      "samtools_percentage_of_properly_paired_reads_...",
                                      "samtools_reads_mapped_percent"  ,                  
                                      "samtools_reads_mapped_and_paired_percent"  ,      
                                      "samtools_reads_properly_paired_percent"   ,       
                                      "samtools_reads_duplicated_percent"     ,                 
                                      "samtools_percent_mapped_X",                        
                                      "samtools_percent_mapped_Y"      ,                 
                                      "cutadapt_percent_trimmed_R2" ,                     
                                      "cutadapt_percent_trimmed_R1",                     
                                      "cutadapt_mean_percent_trimmed" ,          
                                      "Average.input.read.length"  ,                     
                                      "Uniquely.mapped.reads.."   ,                      
                                      "Number.of.splices..GT.AG" ,                       
                                      "Number.of.splices..GC.AG" ,                        
                                      "Number.of.splices..AT.AC" ,                       
                                      "Mismatch.rate.per.base..." ,                   
                                      "Deletion.rate.per.base"   ,                                      
                                      "Insertion.rate.per.base"  ,                        
                                      "X..of.reads.mapped.to.multiple.loci"  ,           
                                      "X..of.reads.mapped.to.too.many.loci"    ,         
                                      "X..of.reads.unmapped..too.many.mismatches"   ,    
                                      "X..of.reads.unmapped..too.short" ,                
                                      "X..of.reads.unmapped..other" ,                    
                                      "X..of.chimeric.reads" ) 

md <- merge(comb2, metrics2, by.x="specimenID", all.x=T)


# add sequencing statistics below ASAP - until then use this
file_path <- "mayo_md_all.csv"
write.csv(md, file = file_path, row.names = FALSE)
file <- synapser::synStore(synapser::File(path = file_path, name="mayo_md_all.csv", parent = "syn75817798"))
