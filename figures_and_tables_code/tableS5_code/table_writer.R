path <- "Z:\\revisions\\revisions2\\cell_classification\\krasnow_lung\\"

summary <- read.csv(paste0(path, "summary.csv"))
write.csv(summary, paste0(path, "summary.csv"))

tbl <- table(summary$ddqc_cluster, summary$cell_typist)
write.csv(tbl, paste0(path, "table_ct_vs_ddqc.csv"))

cl_mad_results <- read.csv(paste0(path, "classification_mad_results.csv"))
write.csv(cl_mad_results, paste0(path, "classification_mad_results.csv"))

summary.filtered <- cl_mad_results[cl_mad_results$ddqc_cluster_passed_qc == "True", ]
tbl.filtered <- table(summary.filtered$ddqc_cluster, summary.filtered$cell_typist)
write.csv(tbl.filtered, paste0(path, "table_ct_vs_ddqc_filtered.csv"))

ct.unique <- cl_mad_results[cl_mad_results$cell_typist_passed_qc == "True" & cl_mad_results$ddqc_cluster_passed_qc == "False",]
ddqc.unique <- cl_mad_results[cl_mad_results$cell_typist_passed_qc == "False" & cl_mad_results$ddqc_cluster_passed_qc == "True",]
write.csv(ct.unique, paste0(path, "cell_typist_unique.csv"))
write.csv(ddqc.unique, paste0(path, "ddqc_unique.csv"))

tmp <- cl_mad_results[cl_mad_results$cell_typist_passed_qc == "False" & cl_mad_results$ddqc_cluster_passed_qc == "True" & cl_mad_results$n_genes < 200,]