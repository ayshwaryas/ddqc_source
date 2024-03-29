# Primes ddqc analysis code
### Code files description:
- Config/
  - read_info/ - csv files describing the datasets and their locations, which are used by pegasus aggregate_matrices function. 
      Each folder inside corresponds to one project, in it each csv corresponds to one tissue describing sample name and relative location of each dataset.  
      All tissues are also included into projects.csv, which contains project and tissue name, human or mouse indication, and path to annotations.
  - config - constants and paths. More on that in running instructions
- figures_code/ - code used to create figures for the paper
  - annotation_comparison.R - finds difference between PanglaoDB annotations and author provided annotations and calculates accuracy.
  - config.R - imports and paths
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
- misc/
  - signatures/ signature genes CSVs
  - cluster_plots.ipynb, cluster_plots-joint.ipynb, dotplots_generator.ipynb, dotplots_generator-joint.ipynb - annotation plots (dotplots and umaps on common genes and markers)
  - convert_gene_names.py - conversion of genes.tsv from ENSG gene codes to gene names (was used for EBI)
  - genes_dec.tsv - matches gene names with ENSG codes
  - ct_markers.csv - common celltype markers
  - markers.tsv - PanglaoDB markers database
  - gene_search.py - tool that takes marker genes from !clusters.csv and finds their mentions in ct_markers.csv and prints possible cell types
  - rds_to_mtx.R - converts RDS object to mtx matrix format
  - top_de_check.py - checks top DE genes for specific tissue and determines whether they are low percent or high percent (was used for annotation)
- R_code/ - old implementation of pipeline in R and Seurat. Probably does not work properly.
- Main pipeline files:
  - filtering.py - pipeline implementation of ddqc and standard cutoff, and filtering statistics recording
  - joint_clustering_old.py - joint clustering of cutoff and ddqc results
  - method_comparison.py - main pipeline that runs for each tissue using basic filtering only, cutoff, and ddqc filtering methods
  - plotting.R - figure3 plots that are generated right after pipeline runs
  - run_all.py - runs all tissues in all projects
  - run_tissue.py - runs specific tissue
  - submit_tissue.sh - schedules run of a specific tissue as an UGER task
  - utils.py - clustering and annotation functions.

    
### Output and Data files description
Original data files are located in the `/ahg/regevdata/projects/scqc/data/`. They are organized by organism, then by project, then by tissue.

All the output generated by the pipeline is located in `/ahg/regevdata/projects/scqc/output_pg`. They are organized by project, then by tissue. Each folder than will have the following directories in it:
- 1.4-cutoff-10 - cutoff approach results
- 1.4-joint_clustering_old - joint clustering results
- 1.4-mad-2 - ddqc results
- 1.4-none-0 - clustering after basic filtering

Each method was clustered with pegasus using resolution 1.4
Each of the method folders will have the following files:
- !cells.csv - metadata for all the cells retained by the approach:
  - barcodekey - cell barcodekey
  - Channel - sample name
  - annotations - cell annotation from publisher (if present)
  - n_genes - number of genes
  - n_counts - number of counts
  - percent_mito - percent of mitochondrial genes
  - percent_ribo - percent of ribosomal genes
  - louvain_labels - cluster labels
  - apoptosis, G1/S, G2/M, cycle_diff, predicted_phase, mito_genes, mito_ribo, ribo_genes, apoptosis1, apoptosis2, apoptosis3, cycling_g2m, cycling_s, dissociation, er_stress, go_mito_resp_chain_human - signature scores
  - pca1, pca2, umap1, umap2 - pca and umap coordinates
- !clusters.csv - information about the clusters
  - cluster - cluster number
  - annotation - most common annotation from publisher in the cluster
  - annotation2 - 2nd most common annotation from publisher in the cluster
  - %annotation, %annotation2 - percent of cells with annotation and annotation2 in the cluster
  - cell_type - cell type determined by PanglaoDB
  - cell_type2 - 2nd most probable cell type as determined by PanglaoDB
  - n_cells - number of cells in the cluster
  - genes_mean, genes_median, mito_mean, mito_median, ribo_mean, ribo_median - cluster QC metrics statistics
  - markers - cluster marker genes separated by ;
  - score_genes - marker genes used by PanglaoDB to call the cell type
- !markers.csv - information about the markers. Result of pegasus de_analysis function.
- box_%metric_name%.pdf - boxplot of the qc metric
- density_%metric_name%.pdf - density plot of the qc metric
- umap_%metric_name%.pdf - umap colored with qc metric
- umap_clusters.pdf - umap colored by cluster number
- violin_%metric_name%.pdf - violin plot of the qc metric

In addition to all of those files 1.4-joint_clustering_old will have the following 2 plots:
- !p_barplot.pdf - barplot of cluster composition based on which cells passed which qc method
- !p_filterplot.pdf - UMAP colored based on which cells passed which qc method

### Running instructions
#### Requirements
Python requirements:
- numpy>1.20
- matplotlib>=3.4.0
- pandas>=1.2.0
- pegasusio
- pegasuspy>=1.3
- seaborn>=0.11

R requirements:
- ggplot
- Seurat
- dplyr

#### Configuration
First, in config/ create a file named local_config.py, with the following code:
`local = False`
Then adjust config files as needed:

##### config/config.py:
- Paths under `if not local:`:
  - DATA_DIR - path to the directory with raw files
  - OUTPUT_DIR - path to the output directory
- resolution - clustering resolution
- do_batch_correction - whether to do batch correction
- basic_genes_filter - basic nGenes filter (performed for all methods)
- basic_mito_filter - same, but for mito
- MITO_PREFIXES, RIBO_PREFIXES - mito and ribo prefixes for human and mouse
- do_%metric_name% - whether to do ddqc on a given metrics
- do not modify parameters under METHOD COMPARISON SCRIPTS

##### Reading the tissue files and adding metadata from the output_pg folder 
Refer to `reading_data_tutorial.ipynb`.

#### Running main pipeline:
##### method_comparison.py:
Launch the script. It will prompt you to enter project. Then enter -1 to enter tissue name afterwards. After that enter method ID (you can find all methods in order in the config file)
You can also run it by calling mc_main function and providing the arguments to the function

##### joint_clustering_old.py:
Same as above, except cutoff and mad methods for the given tissue need to be run prior to running this script.

