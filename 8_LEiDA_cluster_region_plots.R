library(ggplot2)
library(readxl)
library(brainconn)
library(patchwork)
library(pals)
library(dplyr)
library(RNifti)

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


PL_names <- grep("PL State", names(cluster_regions), value = T)

for (i in 2:length(PL_names)){
  
  ## Make binary matrixes of regions in clusters
  PL_index <- cluster_regions[PL_names[i]] > 0
  PL_index[is.na(PL_index)] <- F
  
  binary_matrix <-  as.data.frame(matrix(0, nrow = 264, ncol = 264))
  first_index <- which(PL_index)[1]
  binary_matrix[PL_index,first_index] <- 1
  binary_matrix[first_index,PL_index] <- 1
  binary_matrix[first_index,first_index] <- 0
  node_size <- as.vector(scale(cluster_regions[PL_index, PL_names[i]], center = F) * 7)
  node_colour <- recode(power_atlas$Network[PL_index],
                        "Subcortical/Unknown" = "grey",
                        "Visual" = "#a251ae",
                        "Somatomotor" = "#789ac0",
                        "Dorsal Attention" = "#409832",
                        "Ventral Attention" = "#e166ff",
                        "Limbic" = "#f7feca",
                        "Frontoparietal" = "#f1b845",
                        "Default Mode" = "#da727d")

  ## Plot clusters
  p1 <- brainconn(atlas = power_atlas, 
            conmat=binary_matrix, 
            view="top",
            edge.alpha = 0,
            show.legend = T, 
            node.color = node_colour,
            labels = F,
            label.size = 5,
            background.alpha = 0.5,
            node.size = node_size)
  
  
  power_atlas_filtered <- power_atlas %>%
    filter(PL_index)
  
  new_pal <- recode(sort(unique(power_atlas_filtered$Network)),
                    "Subcortical/Unknown" = "grey",
                    "Visual" = "#a251ae",
                    "Somatomotor" = "#789ac0",
                    "Dorsal Attention" = "#409832",
                    "Ventral Attention" = "#e166ff",
                    "Limbic" = "#f7feca",
                    "Frontoparietal" = "#f1b845",
                    "Default Mode" = "#da727d")
  
  p2 <- ggplot(power_atlas_filtered, aes(x = reorder(ROI.Name, get(PL_names[i])), y = get(PL_names[i]), fill = Network)) +
    geom_bar(stat = "identity") +
    coord_flip() + 
    theme_classic()+
    labs(x = "", 
         y = "",
         title = "") +
    scale_fill_manual(values = new_pal)
  
  assign(paste0("brain_plot_PL_", i),p1)
  assign(paste0("effect_plot_PL_", i),p2)
}



jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/brain_plot_PL_2.jpeg"), width = 1500, height = 1000, res = 200)
brain_plot_PL_2
dev.off()
jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/effect_plot_PL_2.jpeg"), width = 1500, height = 1000, res = 200)
effect_plot_PL_2
dev.off()

jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/brain_plot_PL_3.jpeg"), width = 1500, height = 1000, res = 200)
brain_plot_PL_3
dev.off()
jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/effect_plot_PL_3.jpeg"), width = 1500, height = 2000, res = 200)
effect_plot_PL_3
dev.off()

jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/brain_plot_PL_4.jpeg"), width = 1500, height = 1000, res = 200)
brain_plot_PL_4
dev.off()
jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/effect_plot_PL_4.jpeg"), width = 1500, height = 1000, res = 200)
effect_plot_PL_4
dev.off()

jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/brain_plot_PL_5.jpeg"), width = 1500, height = 1000, res = 200)
brain_plot_PL_5
dev.off()
jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/effect_plot_PL_5.jpeg"), width = 1500, height = 1500, res = 200)
effect_plot_PL_5
dev.off()

jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/brain_plot_PL_6.jpeg"), width = 1500, height = 1000, res = 200)
brain_plot_PL_6
dev.off()
jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/effect_plot_PL_6.jpeg"), width = 1500, height = 1500, res = 200)
effect_plot_PL_6
dev.off()

jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/brain_plot_PL_7.jpeg"), width = 1500, height = 1000, res = 200)
brain_plot_PL_7
dev.off()
jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/effect_plot_PL_7.jpeg"), width = 1500, height = 1000, res = 200)
effect_plot_PL_7
dev.off()

jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/brain_plot_PL_8.jpeg"), width = 1500, height = 1000, res = 200)
brain_plot_PL_8
dev.off()
jpeg(file=paste0("~/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/results/LEiDA/effect_plot_PL_8.jpeg"), width = 1500, height = 2000, res = 200)
effect_plot_PL_8
dev.off()


df <- data_frame(x = 1:8,
                 y = c("Subcortical/Unknown", "Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention", "Limbic", 
                       "Frontoparietal", "Default Mode"))
ggplot(df, aes(x = x, y = y, fill = y)) +
  geom_bar(stat="identity", position = position_stack(reverse = F)) +
  scale_fill_manual(values = c("#da727d", "#409832", "#f1b845", "#f7feca", "#789ac0", "grey", "#e166ff", "#a251ae"))

