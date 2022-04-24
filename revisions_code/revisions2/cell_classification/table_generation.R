df <- read.csv("Z:\\cell_classification\\pbmc_summary.csv")
df$ddqc_cluster[is.na(df$ddqc_cluster)] <- -1
write.csv(table(df$ddqc_cluster, df$single_r), "Z:\\tmp_single_r.csv")
write.csv(table(df$ddqc_cluster, df$azimuth), "Z:\\tmp_azimuth.csv")
write.csv(table(df$ddqc_cluster, df$cell_typist), "Z:\\tmp_cell_typist.csv")
write.csv(df, "Z:\\tmp_summary.csv")