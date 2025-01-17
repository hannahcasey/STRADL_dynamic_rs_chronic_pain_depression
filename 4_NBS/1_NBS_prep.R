library(R.matlab)
library(dplyr)
library(brainGraph)

## Phenotypes ----
STRADL_phenotype_QC_covariates <-  read.csv("/Volumes/GenScotDepression/users/hcasey/rs_fMRI_STRADL_CP_MDD/resources/STRADL_phenotype_QC_covariates.csv", header = T)

## Get list of participants to analyse
ID <- list.files("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/BOLD/LEiDA/power/ICA/", pattern = ".mat")
ID <- gsub("_.*.mat", "", ID)

phenotypes <- left_join(as.data.frame(ID), STRADL_phenotype_QC_covariates, by = "ID")

depression <- data.frame(case = ifelse(phenotypes$current_depression == 1, 1 , 0),
                         control = ifelse(phenotypes$current_depression == 0, 1 , 0),
                         age = phenotypes$Age,
                         sex = ifelse(phenotypes$Sex == "M", 1 , 0),
                         site = ifelse(phenotypes$Site == "DU", 1 , 0)
                         )

pain <- data.frame(case = ifelse(phenotypes$current_pain == 1, 1 , 0),
                   control = ifelse(phenotypes$current_pain == 0, 1 , 0),
                   age = phenotypes$Age,
                   sex = ifelse(phenotypes$Sex == "M", 1 , 0),
                   site = ifelse(phenotypes$Site == "DU", 1 , 0))

comorbidity_Pain_MDD <- data.frame(Pain_Dep = ifelse(phenotypes$comorbidity == "Pain+MDD+", 1 , 0),
                                   noPain_noDep = ifelse(phenotypes$comorbidity == "Pain-MDD-", 1 , 0),
                                   Pain_noDep = ifelse(phenotypes$comorbidity == "Pain+MDD-", 1 , 0),
                                   noPain_Dep = ifelse(phenotypes$comorbidity == "Pain-MDD+", 1 , 0),
                          age = phenotypes$Age,
                          sex = ifelse(phenotypes$Sex == "M", 1 , 0),
                          site = ifelse(phenotypes$Site == "DU", 1 , 0))


write.table(depression, "~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/NBS1.2/STRADL_depression/test_depression_designMatrix.txt", quote = F, row.names = F, col.names = F)
write.table(pain, "~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/NBS1.2/STRADL_depression/pain_designMatrix.txt", quote = F, row.names = F, col.names = F)
write.table(comorbidity_Pain_MDD, "~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/NBS1.2/STRADL_depression/comorbidity_designMatrix.txt", quote = F, row.names = F, col.names = F)


## Power Atlas ----
atlas <- brainGraph::power264 %>%
  select(name)

coords <- brainGraph::power264 %>%
  select(x.mni, y.mni, z.mni)

write.table(atlas, "~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/NBS1.2/power/nodeLabels.txt", quote = F, row.names = F, col.names = F)
write.table(coords, "~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/NBS1.2/power/COG.txt", quote = F, row.names = F, col.names = F)


## Timeseries ----
BOLD_corr_list <- list.files(path="/Volumes/GenScotDepression/users/hcasey/rs_fMRI_STRADL_CP_MDD/correlation_matrix_brainnetome/",
                                   pattern="corrMatrix_atlas-power2011Dseg_desc-correlation_matrix.tsv", full.names = T)

BOLD_corr_df_list = lapply(BOLD_corr_list, function(x)read.table(x, header=F)) 

## Use STRADL IDs as list names
names(BOLD_corr_df_list) <-  gsub("/Volumes/GenScotDepression/users/hcasey/rs_fMRI_STRADL_CP_MDD/correlation_matrix_brainnetome//sub-|_task-resting_feature-corrMatrix_atlas-power2011Dseg_desc-correlation_matrix.tsv|_task-Resting_feature-corrMatrix_atlas-power2011Dseg_desc-correlation_matrix.tsv","", BOLD_corr_list)

## Save timeseries
ts_dir <- "~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/BOLD/NBS/"

for (participant in STRADL_phenotype_QC_covariates$ID){
  
  ## Remove NaN values from matrix
  #BOLD_corr_df_list[[participant]] <- BOLD_corr_df_list[[participant]][,!colnames(BOLD_corr_df_list[[participant]]) %in% NaN_regions]
  
  ## Print "writeMat" commands to be copy and pasted - cannot be run in loop
  cat(paste0("write.table(BOLD_corr_df_list[['", participant, "']]", ",'", ts_dir, participant,  ".txt', quote = F, row.names = F, col.names = F, na = 'NaN')"), "\n")
}

