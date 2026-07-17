## Script overview (execution order)

All analysis scripts are located in the `src/` directory.  
Each script corresponds to a specific step in the human–macaque–marmoset comparative analysis workflow.

### 01_Preprocess_Human_Macaque.Rmd
Human and macaque preprocessing pipeline (QC, normalization, initial filtering).

### 02_Preprocess_Marmoset.Rmd
Marmoset preprocessing pipeline (QC, normalization, initial filtering).

### 03_Marge_Human_Marmoset.Rmd
Integration of human and marmoset datasets.

### 04_Make_Figures.Rmd
Generation of figures used in the manuscript (Fig.1, Fig.5).

### 05_Integrate_ATACdata.Rmd
Integration of ATAC-seq data (human/macaque).

### 06_Integrate_ATAC_Marmoset.Rmd
Integration of marmoset ATAC-seq data.

### 07_Pseudotime_PredRealtime_cellnum.Rmd
Pseudotime analysis and BrdU comparison.

### 08_Trajectory_Combined_DEGdetection.Rmd
Differential expression analysis along the combined trajectory.

### 09_Trajectory_CombinedData_Rezvani_2022_PANDO.Rmd
Integrated pseudotime analysis including Rezvani 2022 PANDO dataset.

### Function_list_combined_pando.r
Utility functions used for the integrated pseudotime analysis.

### conf/
Configuration files used during script execution (paths, parameters, settings).
