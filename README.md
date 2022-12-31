# brainstem_celltype
Understanding the molecular similarities between brainstem and spinal motor circuits will be key to understanding the cellular features that define motor circuits in the brainstem. By integrating scRNA-seq datasets of the mouse brainstem and spinal cord, we try to determine molecular similarities at both the population level and at the level of individual cell types. This repository contains codes for processing scRNA-seq data of brainstem obtained from mice (in-house). We make use of the Russ et al (2021) scRNA-seq spinal cord dataset to learn about the different cell types.

The different scripts used are:

1. levine_explore.R - This script helps to reduce the entire dataset (~53,000 cells) to create a dataset which only contains neurons (~20,000 cells). These neuronal cells are clustered and visualized in a UMAP space. The final dataset is saved as 'integrate2.rds'.

2. integrate_levine.R - This script integrated the ~40,000 brainstem cells dataset with the spinal cord neuronal dataset using the Seurat integration functions. After integration, the combined dataset is treated as a single unit upon which one can perform clustering and saved as 'levine_integrate1_clustered.rds'. A threshold of 5% is set for which clusters are to be used for further analysis (clusters should have atleast 5% of brainstem or spinal cord cells). These clusters are then saved in another file 'levine_integrate1_clustered_thresh5.rds'.

3. mapping_bs.R - This script maps the entire brainstem dataset of ~40,000 cells to the entire spinal cord dataset (~53,000 cells). Essentilaly, one is able to project the cell labels assigned to the spinal cord dataset to the brainstem cells. The brainstem dataset with teh additinal labels of cells in the metadata field is stored as 'bs_levine_label_raw.rds'. These labels are also stored in all the integrated whole and threshold dataset. Now that every cell has a label in the threshold dataset, we can split the dataset by the broad types of Excitatory and Inhibitory types, as well whether they are from brainstem or the spinal cord. Four files are essentially saved, which will be used for the second phase of integration.

4. excit_levine.R - This script integrates the 'excit_bs.rds' and 'excit_sc.rds' into a single dataset which can be clustered and visualized on a UMAP space. The file is labeled as 'excit_integrate1_clustered_raw.rds'.

5. inhib_levine.R - This script integrates the 'inhib_bs.rds' and 'inhib_sc.rds' into a single dataset which can be clustered and visualized on a UMAP space. The file is labeled as 'inhib_integrate1_clustered_raw.rds'.

6. excit_marker_analysis.R - This script renames the clusters as a majority of the cells found and performs hierarchal clustering. A particular node along with the leaves below it can be submitted as a collective for differential expression testing (marker gene analysis). This script is used only for the excitatory dataset. 

7. inhib_marker_analysis.R - This script performs the same functions as the previous script, but for the inhibitory dataset. 

Scripts for figures: 

8. bs_fig.R: This script creates the plots for all brainstem figures: 
  a. histograms of UMIs per cell and genes per cell
  b. clustering the brainstem data and grouping by neurotransmitter status
  c. applying the labels to the brainstem dataset after mapping - contains lamina information as well as value of predicted scores

9. sc_fig.R: This script creates the plots for all the spinal cord figures: 
  a. histograms of UMIs per cell and genes per cell
  b. clustering the spinal cord data and grouping by neurotransmitter status, lamina and lineage information.

10. integrated_fig.R: This script creates all the figures for integrated datasets: 
  a. barplot representing the proportion of the brainstem and spinal cord cells in each cluster
  b. excitatory dataset clustered by cell types and lamina information 
  c. inhibitory dataset clustered by cell types and lamina information 

11. heatmap.R: This script plots the gene expression of the 16 marker genes chosen for the clusters in the ‘integrated_levine_label_thresh5.rds’. 


Files for additional functions: 

12. functions.R: Contains functions for different normalisation methods for gene expression, heat map construction function and a function to provide discrete colors for clusters in UMAP/tSNE plots. 

13. tree_functions_seurat.R: Contains functions taken from Seurat and a function written to get all terminal leaf nodes which is useful in marker analysis scripts when performing hierarchical clustering.  


