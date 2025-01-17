library(dplyr)


## Get STRADL phenotypes
data <-  read.csv("/Volumes/GenScotDepression/users/hcasey/rs_fMRI_STRADL_CP_Dep/resources/STRADL_phenotype_QC_covariates.csv", header = T)

## Descriptive statistics of age and sex in outcomes ----
## Entire sample
n_full <- nrow(data)
age_full <- round(mean(data$Age, na.rm = T), 2)
n_females_full <- sum(data$Sex == "F", na.rm = T)
prop_females_full <- round(n_females_full/nrow(data), 4) *100

## No Pain
n_no_pain <- sum(data$current_pain == 0, na.rm = T)
prop_no_pain <- round(n_no_pain/sum(!is.na(data$current_pain), na.rm = T), 4) *100
age_no_pain <- round(mean(data$Age[data$current_pain == 0], na.rm = T), 2)
n_females_no_pain <- sum(data$current_pain == 0 & data$Sex == "F", na.rm = T)
prop_females_no_pain <- round(n_females_no_pain/sum(data$current_pain == 0, na.rm = T), 4) *100

## Pain
n_pain <- sum(data$current_pain == 1, na.rm = T)
prop_pain <- round(n_pain/sum(!is.na(data$current_pain), na.rm = T), 4) *100
age_pain <- round(mean(data$Age[data$current_pain == 1], na.rm = T), 2)
n_females_pain <- sum(data$current_pain == 1 & data$Sex == "F", na.rm = T)
prop_females_pain <- round(n_females_pain/sum(data$current_pain == 1, na.rm = T), 4) *100

## No Depression
n_no_depression <- sum(data$QIDS < 6 , na.rm = T)
prop_no_depression <- round(n_no_depression/sum(!is.na(data$QIDS), na.rm = T), 4) *100
age_no_depression <- round(mean(data$Age[data$QIDS < 6], na.rm = T), 2)
n_females_no_depression <- sum(data$QIDS < 6 & data$Sex == "F", na.rm = T)
prop_females_no_depression <- round(n_females_no_depression/sum(data$QIDS < 6, na.rm = T), 4) *100

## Depression
n_depression <- sum(data$QIDS > 5 , na.rm = T)
prop_depression <- round(n_depression/sum(!is.na(data$QIDS), na.rm = T), 4) *100
age_depression <- round(mean(data$Age[data$QIDS > 5], na.rm = T), 2)
n_females_depression <- sum(data$QIDS > 5 & data$Sex == "F", na.rm = T)
prop_females_depression <- round(n_females_depression/sum(data$QIDS> 5, na.rm = T), 4) *100

## Pain-Dep-
n_noCPnoDep <- sum(data$comorbidity == "Pain-Dep-", na.rm = T)
prop_noCPnoDep <- round(n_noCPnoDep/length(na.omit(data$comorbidity)), 4) *100
age_noCPnoDep <- round(mean(data$Age[data$comorbidity == "Pain-Dep-"], na.rm = T), 2)
n_females_noCPnoDep <- sum(data$comorbidity == "Pain-Dep-" & data$Sex == "F", na.rm = T)
prop_females_noCPnoDep <- round(n_females_noCPnoDep/sum(data$comorbidity == "Pain-Dep-", na.rm = T), 4) *100
## Pain+Dep-
n_CPnoDep <- sum(data$comorbidity == "Pain+Dep-", na.rm = T)
prop_CPnoDep <- round(n_CPnoDep/length(na.omit(data$comorbidity)), 4) *100
age_CPnoDep <- round(mean(data$Age[data$comorbidity == "Pain+Dep-"], na.rm = T), 2)
n_females_CPnoDep <- sum(data$comorbidity == "Pain+Dep-" & data$Sex == "F", na.rm = T)
prop_females_CPnoDep <- round(n_females_CPnoDep/sum(data$comorbidity == "Pain+Dep-", na.rm = T), 4) *100
## Pain-Dep+
n_noCPDep <- sum(data$comorbidity == "Pain-Dep+", na.rm = T)
prop_noCPDep <- round(n_noCPDep/length(na.omit(data$comorbidity)), 4) *100
age_noCPDep <- round(mean(data$Age[data$comorbidity == "Pain-Dep+"], na.rm = T), 2)
n_females_noCPDep <- sum(data$comorbidity == "Pain-Dep+" & data$Sex == "F", na.rm = T)
prop_females_noCPDep <- round(n_females_noCPDep/sum(data$comorbidity == "Pain-Dep+", na.rm = T), 4) *100
## Pain+Dep+
n_CPDep <- sum(data$comorbidity == "Pain+Dep+", na.rm = T)
prop_CPDep <- round(n_CPDep/length(na.omit(data$comorbidity)), 4) *100
age_CPDep <- round(mean(data$Age[data$comorbidity == "Pain+Dep+"], na.rm = T), 2)
n_females_CPDep <- sum(data$comorbidity == "Pain+Dep+" & data$Sex == "F", na.rm = T)
prop_females_CPDep <- round(n_females_CPDep/sum(data$comorbidity == "Pain+Dep+", na.rm = T), 4) *100

## Combine descriptive statistics
descriptive_statistics <- data.frame(Group = c("Full Sample", 
                                                 "No Pain", "Pain", 
                                                 "No Depression", "Depression",
                                                 "No Pain + No Depression","No Pain + Depression",  
                                                 "Pain + No Depression", "Pain + Depression"),
                                     "N (% of sample with data)" = c(
                                       paste0(prettyNum(n_full, big.mark = ","), " (100%)"),
                                       paste0(prettyNum(n_no_pain, big.mark = ","), " (", prop_no_pain, "%)"),
                                       paste0(prettyNum(n_pain, big.mark = ","), " (", prop_pain, "%)"),
                                       paste0(prettyNum(n_no_depression, big.mark = ","), " (", prop_no_depression, "%)"),
                                       paste0(prettyNum(n_depression, big.mark = ","), " (", prop_depression, "%)"),
                                       paste0(prettyNum(n_noCPnoDep, big.mark = ","), " (", prop_noCPnoDep, "%)"),
                                       paste0(prettyNum(n_noCPDep, big.mark = ","), " (", prop_noCPDep, "%)"),
                                       paste0(prettyNum(n_CPnoDep, big.mark = ","), " (", prop_CPnoDep, "%)"),
                                       paste0(prettyNum(n_CPDep, big.mark = ","), " (", prop_CPDep, "%)")),
                                     "Mean Age" = c(age_full, 
                                                    age_no_pain, age_pain, 
                                                    age_no_depression, age_depression,
                                                    age_noCPnoDep, age_noCPDep, age_CPnoDep, age_CPDep),
                                     "N females (% of cases)" = c(
                                       paste0(prettyNum(n_females_full, big.mark = ","), " (", prop_females_full, "%)"),
                                       paste0(prettyNum(n_females_no_pain, big.mark = ","), " (", prop_females_no_pain, "%)"),
                                       paste0(prettyNum(n_females_pain, big.mark = ","), " (", prop_females_pain, "%)"),
                                       paste0(prettyNum(n_females_no_depression, big.mark = ","), " (", prop_females_no_depression, "%)"),
                                       paste0(prettyNum(n_females_depression, big.mark = ","), " (", prop_females_depression, "%)"),
                                       paste0(prettyNum(n_females_noCPnoDep, big.mark = ","), " (", prop_females_noCPnoDep, "%)"),
                                       paste0(prettyNum(n_females_noCPDep, big.mark = ","), " (", prop_females_noCPDep, "%)"),
                                       paste0(prettyNum(n_females_CPnoDep, big.mark = ","), " (", prop_females_CPnoDep, "%)"),
                                       paste0(prettyNum(n_females_CPDep, big.mark = ","), " (", prop_females_CPDep, "%)"))
)


# ## Add column indicating statistical significance ----
# CP_Dep_results <- CP_Dep_results %>%
#   mutate(significant = ifelse(CP_Dep_results$p_adjust < 0.05,
#                               "Yes", "No"))
# 
# CPDep_results <- CPDep_results %>%
#   mutate(significant = ifelse(CPDep_results$p_adjust < 0.05,
#                               "Yes", "No"))

## Save results
write.csv(descriptive_statistics,"/Users/hannahcasey/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/phenotype/descriptive_statistics.csv")

