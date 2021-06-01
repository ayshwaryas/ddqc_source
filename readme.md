# Primes ddqc analysis code
### Code files description:
- Config/
  - read_info/ - csv files describing the datasets and their locations, which are used by pegasus aggregate_matrices function. 
      Each folder inside corresponds to one project, in it each csv corresponds to one tissue describing sample name and relative location of each dataset.  
      All tissues are also included into projects.csv, which contains project and tissue name, human or mouse indication, and path to annotations.
  - config - constants and paths. More on that in running instructions
- figures_code/ - code used to create figures for the paper
  - annotation_comparison.R - finds difference between PanglaoDB annotations and author provided annotations and calculates accuracy.
  - config.R - imports and paths. More on that in running instructions
  - count_cells.R - counts total human and mouse cells
  - ddqc_vs_paper.R - compares ddqc results vs. paper results in TM
  - examples_search_new_joint_pg.R - searches for unique clusters (compared to standard cutoff)
  - figure1_R_code.R, figure2_code.R, figure3_remake.R, figure4.R - code for corresponding figure
  - figure3_table.R - table with detailed info and percentages for each tissue (Table S4)
  - paper_vs_ddqc_UMAP.R - umap for comparing ddqc and paper
  - percent_comparison_table.R - table with percentage difference in number of cells between cutoff and ddqc
  - scrublet_scater_plots.R, signatures_scatter_code.R, signatures_scatter_high_mito_code.R - scrublet and signature scores based plots
  - tm_lung_heatmaps.ipynb - score violins and heatmaps for TM lung
  - trends_percentile_search.R, trends_search.R - trends tables
- Main pipeline files:
  - filtering.py - pipeline implementation of ddqc and standard cutoff, and filtering statistics recording
  - joint_clustering_old.py - joint clustering of cutoff and ddqc results
  - method_comparison.py - main pipeline that runs for each tissue using basic filtering only, cutoff, and ddqc filtering methods
  - plotting.R - figure3 plots that are generated right after pipeline runs
  - run_all.py - runs all tissues in all projects
  - run_tissue.py - runs specific tissue
  - submit_tissue.sh - schedules run of a specific tissue as an UGER task
  - utils.py - clustering and annotation functions.

    
### Output files description

### Running instructions
