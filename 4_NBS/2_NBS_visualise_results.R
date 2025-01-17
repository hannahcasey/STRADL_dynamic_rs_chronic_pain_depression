#remotes::install_github("sidchop/brainconn")
library(brainconn)
library(readxl)
library(dplyr)
library(R.matlab)


## For each region in power atlas, map to Yeo brain network
## http://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/tutorial_packages/OSX/freesurfer/average/Yeo_JNeurophysiol11_MNI152/
network_atlas <- readNifti("~/Desktop/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii")
mni_coords <- as.matrix(brainGraph::power264[,c("x.mni", "y.mni", "z.mni")])

get_network <- function(x){
  x_mni <- x[1]
  y_mni <- x[2]
  z_mni <- x[3]
  # Convert MNI coordinates to voxel indices using the sform/qform matrix
  voxel_coords <- worldToVoxel(c(x_mni, y_mni, z_mni), network_atlas)
  # Round the voxel coordinates to nearest integers since voxel indices must be integers
  voxel_coords <- round(voxel_coords)
  # Lookup the voxel value at the computed indices
  voxel_value <- network_atlas[voxel_coords[1], voxel_coords[2], voxel_coords[3],1]
}

power_networks <- apply(mni_coords, 1, get_network)

## Load in LEiDA clusters
cluster_regions <- read_excel("/Users/hannahcasey/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/writeup/power_atlas_info.xlsx")

## Plot each cluster by network
power_atlas <- cluster_regions %>%
  rename(ROI.Name = brainGraph,
         x.mni = X,
         y.mni = Y,
         z.mni = Z) %>% 
  mutate(Network = recode(power_networks,
                          `0` = "Subcortical/Unknown",
                          `1` = "Visual",
                          `2` = "Somatomotor",
                          `3` = "Dorsal Attention",
                          `4` = "Ventral Attention",
                          `5` = "Limbic",
                          `6` = "Frontoparietal",
                          `7` = "Default Mode")) %>%
  mutate_at(vars(x.mni, y.mni, z.mni), funs(as.integer(.)))



## Load in NBS results
NBS_depression_decrease_con_mat_3_5 <- readMat("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/NBS/depression_decrease_3_5.mat", )
NBS_depression_decrease_con_mat_3_5 <- summary(NBS_depression_decrease_con_mat_3_5[["nbs"]][[2]][[4]][[1]][[1]])

binary_matrix_3_5 <- as.data.frame(matrix(0, nrow = 264, ncol = 264))

for (i in 1:length(NBS_depression_decrease_con_mat_3_5$i)) {
  
    binary_matrix_3_5[NBS_depression_decrease_con_mat_3_5$i[i],NBS_depression_decrease_con_mat_3_5$j[i]] <- 1
    binary_matrix_3_5[NBS_depression_decrease_con_mat_3_5$j[i],NBS_depression_decrease_con_mat_3_5$i[i]] <- 1
}

NBS_depression_decrease_con_mat_3_7 <- readMat("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/NBS/depression_decrease_3_7.mat", )
NBS_depression_decrease_con_mat_3_7_test_stats <- NBS_depression_decrease_con_mat_3_7[["nbs"]][[2]][[6]]
NBS_depression_decrease_con_mat_3_7 <- summary(NBS_depression_decrease_con_mat_3_7[["nbs"]][[2]][[4]][[1]][[1]])

binary_matrix_3_7 <- as.data.frame(matrix(0, nrow = 264, ncol = 264))

for (i in 1:length(NBS_depression_decrease_con_mat_3_7$i)) {
  
  binary_matrix_3_7[NBS_depression_decrease_con_mat_3_7$i[i],NBS_depression_decrease_con_mat_3_7$j[i]] <- 1
  binary_matrix_3_7[NBS_depression_decrease_con_mat_3_7$j[i],NBS_depression_decrease_con_mat_3_7$i[i]] <- 1
}

NBS_depression_decrease_con_mat_3_8 <- readMat("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/NBS/depression_decrease_3_8.mat", )
NBS_depression_decrease_con_mat_3_8 <- summary(NBS_depression_decrease_con_mat_3_8[["nbs"]][[2]][[4]][[1]][[1]])

binary_matrix_3_8 <- as.data.frame(matrix(0, nrow = 264, ncol = 264))

for (i in 1:length(NBS_depression_decrease_con_mat_3_8$i)) {
  
  binary_matrix_3_8[NBS_depression_decrease_con_mat_3_8$i[i],NBS_depression_decrease_con_mat_3_8$j[i]] <- 1
  binary_matrix_3_8[NBS_depression_decrease_con_mat_3_8$j[i],NBS_depression_decrease_con_mat_3_8$i[i]] <- 1
}


## Plot findings
node_colour <- recode(power_atlas$Network[sort(unique(c(NBS_depression_decrease_con_mat_3_5$i, NBS_depression_decrease_con_mat_3_5$j)))],
                      "Subcortical/Unknown" = "grey",
                      "Visual" = "#a251ae",
                      "Somatomotor" = "#789ac0",
                      "Dorsal Attention" = "#409832",
                      "Ventral Attention" = "#e166ff",
                      "Limbic" = "#f7feca",
                      "Frontoparietal" = "#f1b845",
                      "Default Mode" = "#da727d")

NBS_depression_3_5_plot <- brainconn(atlas = power_atlas, 
          conmat=binary_matrix_3_5, 
          view="ortho",
          node.size = 5,
          edge.alpha = 0.8,
          edge.width = 0.2,
          show.legend = T,
          node.color = c(node_colour),
          background.alpha = 0.5)

node_colour <- recode(power_atlas$Network[sort(unique(c(NBS_depression_decrease_con_mat_3_7$i, NBS_depression_decrease_con_mat_3_7$j)))],
                      "Subcortical/Unknown" = "grey",
                      "Visual" = "#a251ae",
                      "Somatomotor" = "#789ac0",
                      "Dorsal Attention" = "#409832",
                      "Ventral Attention" = "#e166ff",
                      "Limbic" = "#f7feca",
                      "Frontoparietal" = "#f1b845",
                      "Default Mode" = "#da727d")

NBS_depression_3_7_plot <- brainconn(atlas =power_atlas, 
          conmat=binary_matrix_3_7, 
          view="ortho",
          node.size = 5,
          edge.alpha = 0.8,
          edge.width = 0.2,
          node.color = node_colour,
          label.size = 2,
          background.alpha = 0.5)

NBS_depression_3_7_labels_plot <- brainconn(atlas =power_atlas, 
                                     conmat=binary_matrix_3_7, 
                                     view="ortho",
                                     node.size = 5,
                                     edge.alpha = 0.8,
                                     edge.width = 0.2,
                                     show.legend = T,
                                     labels = T, 
                                     label.size = 4,
                                     node.color = node_colour,
                                     background.alpha = 0.5)


node_colour <- recode(power_atlas$Network[sort(unique(c(NBS_depression_decrease_con_mat_3_8$i, NBS_depression_decrease_con_mat_3_8$j)))],
                      "Subcortical/Unknown" = "grey",
                      "Visual" = "#a251ae",
                      "Somatomotor" = "#789ac0",
                      "Dorsal Attention" = "#409832",
                      "Ventral Attention" = "#e166ff",
                      "Limbic" = "#f7feca",
                      "Frontoparietal" = "#f1b845",
                      "Default Mode" = "#da727d")

NBS_depression_3_8_plot <- brainconn(atlas =power_atlas, 
                                     conmat=binary_matrix_3_8, 
                                     view="ortho",
                                     node.size = 5,
                                     edge.alpha = 0.8,
                                     edge.width = 0.2,
                                     show.legend = T,
                                     node.color = node_colour,
                                     background.alpha = 0.5)

## Save plots

# jpeg(file="~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/NBS/depression_decrease_3_5.jpeg", width = 2000, height = 2000, res = 200)
# NBS_depression_3_5_plot
# dev.off()
# 
# jpeg(file="~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/NBS/depression_decrease_3_7.jpeg",  width = 2000, height = 2000, res = 200)
# NBS_depression_3_7_plot
# dev.off()
# 
# jpeg(file="~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/NBS/depression_decrease_3_7_labels.jpeg",  width = 2500, height = 2500, res = 200)
# NBS_depression_3_7_labels_plot
# dev.off()
# 
# jpeg(file="~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/NBS/depression_decrease_3_8.jpeg",  width = 2000, height = 2000, res = 200)
# NBS_depression_3_8_plot
# dev.off()

## Save test statistics
test_stats <- NBS_depression_decrease_con_mat_3_7
test_stats$test_stat <- NA

for (row in 1:nrow(test_stats)){
  
  ## Get test statistic
  test_stats$test_stat[row] <- NBS_depression_decrease_con_mat_3_7_test_stats[as.numeric(test_stats[row,]["j"]), as.numeric(NBS_depression_decrease_con_mat_3_7[row,]["i"])]
  
  ## Replace atlas index with region names
  test_stats$i[row] <- cluster_regions$`Node Names`[as.numeric(test_stats$i[row])] 
  test_stats$j[row] <- cluster_regions$`Node Names`[as.numeric(test_stats$j[row])] 
  
}

names(test_stats) <- c("Node 1", "Node 2","x", "t")

write_csv(test_stats, "~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/writeup/NBS_stats.csv")







