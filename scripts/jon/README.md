# README.md

Overall Goals of Analysis:

  1 - Pseudobulk comparisons and ranking
  
  2 - Identify Trajectories that are comparable to ground truth 
  
  3 - Identify gene changes associated with each trajectory and different between Healthy and T21
  
  4 - Overlay ChromVAR transcription factor binding accessibility onto the trajectories 
  
  5 - Use MIRA to identify how accessibility and gene expression vary between Healthy and T21


Primary Tools Used:
  
  - scanpy
  - pegasus
  - CellRank
  - MIRA
  
  
Description of scripts:

  "GIT_Initial_Data_Analysis_Trajectory_Calculation.ipynb"
  
    - Data processing and generation of PAGA, UMAP, and FLE layouts 
    
   "GIT_Myeloid_Cell_Rank_Trajectories.ipynb"
   
    - Using Cell Rank to identify trajectories 
    - Overlay ChromVAR data onto the trajectories
