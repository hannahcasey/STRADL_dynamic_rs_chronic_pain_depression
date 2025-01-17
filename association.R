library(R.matlab)


LEiDA_results_dir <- "/Users/hannahcasey/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/extracted_results/"

LEiDA_data_info <-  read.csv(paste0(LEiDA_results_dir, "data_info.csv"), header = T)
STRADL_LEiDA_IDs <- gsub("_control.mat|_case.mat", "", LEiDA_data_info$name)

dwell_time <- readMat("/Users/hannahcasey/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/Power/ICA/run_1/LEiDA_Stats_DwellTime.mat")
fractional_occupancy <- readMat("/Users/hannahcasey/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/Power/ICA/run_1/LEiDA_Stats_FracOccup.mat")

STRADL_phenotype_QC_covariates <-  read.csv("/Volumes/GenScotDepression/users/hcasey/rs_fMRI_STRADL_CP_MDD/resources/STRADL_phenotype_QC_covariates.csv", header = T)


## Function to analyze LEiDA outputs between binary groups
binary_assoc <- function(dependent, independent, covariates, dataset, type){
  
  ## Input:
  ## dependent: names of dependent variables
  ## independent: name of independent variable
  ## covariates: variables to adjust for
  ## dataset: df contating data
  ## category: type of outcome
  
  ## Create dataframe to store association analysis output
  output <-  data.frame(dependent_var=character(), dependent_cat=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                        p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric(), n=numeric())
  
  for (dependent_index in 1:length(dependent)){ ## Iterate through dependent variables
    
    ## Define model
    f <- paste0(dependent[dependent_index], " ~ ", independent,  " + " ,  paste(covariates, collapse = " + "))
    
    res <- glm(as.formula(f),
               data = dataset,
               na.action = na.exclude)
    
    output[dependent_index,"dependent_var"] <- dependent[dependent_index]
    output[dependent_index, "beta"] <- summary(res)[["coefficients"]][independent, "Estimate"]
    output[dependent_index, "std"] <- summary(res)[["coefficients"]][independent, "Std. Error"]
    output[dependent_index, "p.value"] <- summary(res)[["coefficients"]][independent, "Pr(>|t|)"]
    output[dependent_index, c("Lower_95CI", "Upper_95CI")] <- confint(res)[independent,]
    output[dependent_index, c("n")] <- nobs(res)
  }
  
  ## Indicate category of LEiDA output
  output$dependent_cat <- type
  
  ## Adjust p-value
  output$p.adjust <- p.adjust(output$p.value, method = "fdr")
  
  return(output)
}


depression_DT_results <-  data.frame(dependent_var=character(), dependent_cat=character(), groups=character(), beta=numeric(), std=numeric(), p.value=numeric(),p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric(), n=numeric(), K=numeric())
depression_FO_results <-  data.frame(dependent_var=character(), dependent_cat=character(), groups=character(), beta=numeric(), std=numeric(), p.value=numeric(),p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric(), n=numeric(), K=numeric())
pain_DT_results <-  data.frame(dependent_var=character(), dependent_cat=character(), groups=character(), beta=numeric(), std=numeric(), p.value=numeric(),p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric(), n=numeric(), K=numeric())
pain_FO_results <-  data.frame(dependent_var=character(), dependent_cat=character(), groups=character(), beta=numeric(), std=numeric(), p.value=numeric(),p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric(), n=numeric(), K=numeric())

for (K in 2:20){
  
  DT <- dwell_time$LT[,K - 1,1:K]; ## Get dwell times for selected K
  FO <- fractional_occupancy$P[,K - 1,1:K]; ## Get fractional occupancy for selected K
  
  DT_df <- cbind(as.data.frame(STRADL_LEiDA_IDs), DT)
  FO_df <- cbind(as.data.frame(STRADL_LEiDA_IDs), FO)
  
  names(DT_df) <- c("ID", paste0("dwell_time_", 1:K))
  names(FO_df) <- c("ID", paste0("fractional_occupancy_", 1:K))
  
  joined_DT_df <- left_join(DT_df, STRADL_phenotype_QC_covariates, by = "ID")
  joined_FO_df <- left_join(FO_df, STRADL_phenotype_QC_covariates, by = "ID")
  
  
  depression_DT <- binary_assoc(
    dependent = paste0("dwell_time_", 2:K),
    independent = "current_depression",
    covariates = c("Sex", "Site", "Age", "Age_squared"),
    dataset = joined_DT_df,
    type = "Fractional Occupancy"
  )
  
  pain_DT <- binary_assoc(
    dependent = paste0("dwell_time_", 2:K),
    independent = "current_pain",
    covariates = c("Sex", "Site", "Age", "Age_squared"),
    dataset = joined_DT_df,
    type = "Fractional Occupancy"
  )
  
  depression_FO <- binary_assoc(
    dependent = paste0("fractional_occupancy_", 2:K),
    independent = "current_depression",
    covariates = c("Sex", "Site", "Age", "Age_squared"),
    dataset = joined_FO_df,
    type = "Fractional Occupancy"
  )
  
  pain_FO <- binary_assoc(
    dependent = paste0("fractional_occupancy_", 2:K),
    independent = "current_pain",
    covariates = c("Sex", "Site", "Age", "Age_squared"),
    dataset = joined_FO_df,
    type = "Fractional Occupancy"
  )
  
  depression_DT$K <- K
  depression_FO$K <- K
  
  pain_DT$K <- K
  pain_FO$K <- K
  
  depression_DT_results <- rbind(depression_DT_results, depression_DT)
  depression_FO_results <- rbind(depression_FO_results, depression_FO)
  
  pain_DT_results <- rbind(pain_DT_results, pain_DT)
  pain_FO_results <- rbind(pain_FO_results, pain_FO)
}

depression_DT_results$p.adjust_all <- p.adjust(depression_DT_results$p.value, method = "fdr")
depression_FO_results$p.adjust_all <- p.adjust(depression_FO_results$p.value, method = "fdr")
pain_DT_results$p.adjust_all <- p.adjust(pain_DT_results$p.value)
pain_FO_results$p.adjust_all <- p.adjust(pain_FO_results$p.value)





## New analysis plan: look at all associations and focus of ones significant across all Ks