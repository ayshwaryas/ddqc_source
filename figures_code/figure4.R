source("figures_code/config.R")


results.dir <- paste0(PATH, "figure4_plots/")
dir.create(results.dir)
results.dir1 <- paste0(results.dir, "panel1/")
results.dir2 <- paste0(results.dir, "panel2/")
results.dir3 <- paste0(results.dir, "panel3/")
dir.create(results.dir1)
dir.create(results.dir2)
dir.create(results.dir3)

plot.cols <- c("Cutoff" = "#8AC4D0", "ddqc" = "#2B2E4A")

data <- read.csv(paste0(PATH, "figure3_table.csv"))
plot.data <- data %>% select(Directory, n_cells.after.ddqc, n_cells.after.cutoff) %>% gather("method", "n.cells", -Directory)
plot.data$Directory <- gsub("/ahg/regevdata/projects/scqc/output_pg/", "", plot.data$Directory)
plot.data$project <- sapply(strsplit(plot.data$Directory, "/"), "[[", 1)
plot.data$tissue <- sapply(strsplit(plot.data$Directory, "/"), "[[", 2)
plot.data[plot.data=="n_cells.after.ddqc"] <- "ddqc"
plot.data[plot.data=="n_cells.after.cutoff"] <- "Cutoff"


for (prj in unique(plot.data$project)) {
  project.data = subset(plot.data, project == prj) %>% arrange(n.cells)
  project.data$tissue = with(project.data, reorder(tissue, -n.cells, sum))
  pl <- ggplot(project.data, aes(x=tissue, y=n.cells, fill=method)) + 
    geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)
  pl_log <- ggplot(project.data, aes(x=tissue, y=log2(n.cells), fill=method)) + 
    geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)
  
  n.tissuess <- length(project.data$tissue)
  ggsave1(paste0(results.dir1, prj, ".pdf"), pl, n.tissuess)
  ggsave1(paste0(results.dir1, prj, "_log.pdf"), pl_log, n.tissuess)
  
  if (prj == "human_other") {
    project.data = subset(project.data, tissue != "Heart_Circulation") %>% arrange(n.cells)
    pl <- ggplot(project.data, aes(x=tissue, y=n.cells, fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)
    pl_log <- ggplot(project.data, aes(x=tissue, y=log2(n.cells), fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)
    
    n.tissuess <- length(project.data$tissue)
    ggsave1(paste0(results.dir1, prj, "_no_hc.pdf"), pl, n.tissuess)
    ggsave1(paste0(results.dir1, prj, "_log_no_hc.pdf"), pl_log, n.tissuess)
  }
  
  if (prj == "ebi_sc_experiment") {
    project.data = subset(project.data, tissue != "E-MTAB-6701-fetal-maternal_interface") %>% arrange(n.cells)
    pl <- ggplot(project.data, aes(x=tissue, y=n.cells, fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)
    pl_log <- ggplot(project.data, aes(x=tissue, y=log2(n.cells), fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)
    
    n.tissuess <- length(project.data$tissue)
    ggsave1(paste0(results.dir1, prj, "_no_fmi.pdf"), pl, n.tissuess)
    ggsave1(paste0(results.dir1, prj, "_log_no_fmi.pdf"), pl_log, n.tissuess)
  }
  
  if (prj == "htapp") {
    project.data[project.data=="lung_cancer"] <- "lung-cancer"
    project.data[["tumor"]] <- sapply(strsplit(as.character(project.data$tissue), "_"), "[[", 1)
    project.data[["type"]] <- sapply(strsplit(as.character(project.data$tissue), "_"), "[[", 2)
    tumor.table <- table(project.data$tumor)
    tumor.table["ovarian"] <- 2
    tumor.table["sarcoma"] <- 2
    project.data.both <- subset(project.data, tumor.table[project.data$tumor] / 2 > 1)
    
    pl <- ggplot(project.data, aes(x=type, y=n.cells, fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols) + 
      facet_wrap(tumor~.,scales="free_y")
    pl_log <- ggplot(project.data, aes(x=type, y=log2(n.cells), fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols) +
      facet_wrap(tumor~.,scales="free_y")
    
    pl.both <- ggplot(project.data.both, aes(x=type, y=n.cells, fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols) + 
      facet_wrap(tumor~.,scales="free_y")
    pl_log.both <- ggplot(project.data.both, aes(x=type, y=log2(n.cells), fill=method)) + 
      geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols) +
      facet_wrap(tumor~.,scales="free_y")
    
    n.tissuess <- length(project.data$tissue)
    ggsave1(paste0(results.dir2, prj, "_facet.pdf"), pl, n.tissuess)
    ggsave1(paste0(results.dir2, prj, "_log_facet.pdf"), pl_log, n.tissuess)
    ggsave1(paste0(results.dir2, prj, "_both_facet.pdf"), pl.both, n.tissuess)
    ggsave1(paste0(results.dir2, prj, "_both_log_facet.pdf"), pl_log.both, n.tissuess)
  }
}

project.data = subset(plot.data, project %in% c("tabula_muris", "tabula_muris_smartseq2", "mouse_cell_atlas")) %>% arrange(n.cells)
project.data[project.data=="tabula_muris"] <- "10x"
project.data[project.data=="tabula_muris_smartseq2"] <- "smartseq2"
project.data[project.data=="mouse_cell_atlas"] <- "microwellseq"

pl <- ggplot(project.data, aes(x=project, y=n.cells, fill=method)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols) + 
  facet_wrap(tissue~.,scales="free_y")
pl_log <- ggplot(project.data, aes(x=project, y=log2(n.cells), fill=method)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols) + 
  facet_wrap(tissue~.,scales="free_y")

project.data <- project.data %>% group_by(project, method) %>% summarize_at("n.cells", sum)

pl_agg <- ggplot(project.data, aes(x=project, y=n.cells, fill=method)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)
pl_agg_log <- ggplot(project.data, aes(x=project, y=log2(n.cells), fill=method)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme_horizontal_with_legend + scale_fill_manual(values = plot.cols)

ggsave1(paste0(results.dir3, "mouse_by_technology.pdf"), pl, 30)
ggsave1(paste0(results.dir3, "mouse_by_technology_log.pdf"), pl_log, 30)
ggsave1(paste0(results.dir3, "mouse_by_technology_aggregated.pdf"), pl_agg, 30)
ggsave1(paste0(results.dir3, "mouse_by_technology_aggregated_log.pdf"), pl_agg_log, 30)
