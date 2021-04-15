source("figures_code/config.R")


human.projects <- c("ebi_sc_experiment", "htapp", "human_other", "human_tissue_atlas",
                    "PanglaoDB")
mouse.projects <- c("mouse_cell_atlas", "tabula_muris", "tabula_muris_smartseq2",
                    "tabula_senis_24m", "tabula_senis_30m", "mouse_other")

hc <- 0 
mc <- 0
for (prj in list.dirs(OUTPUT.DIR, full.names = FALSE, recursive = FALSE)) {
  dataset <- read.csv(paste0(PATH, "figure1_plots/", prj, "/", "!cells.csv"))
  if (prj %in% human.projects) {
    hc <- hc + length(dataset$barcodekey)
  } else  if (prj %in% mouse.projects) {
    mc <- mc + length(dataset$barcodekey)
  } else {
    warning(paste0("Unknown project ", prj))
  }
}

print(paste0("human cells: ", hc, "; mouse cells: ", mc))


