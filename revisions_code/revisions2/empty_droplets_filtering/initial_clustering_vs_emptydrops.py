import pandas as pd


FDR_CUTOFF = 0.01

species = "mouse"
project = "tabula_muris"
tissue = "Heart_and_Aorta"
data_path = f"Z:\\data\\{species}\\{project}\\{tissue}"
output_path = f"Z:\\output_pg\\{project}\\{tissue}\\1.4-mad-2"

sample = "Heart_and_Aorta-10X_P7_4"
emptydrops_results = pd.read_csv(f"{data_path}\\emptydrops\\{sample}.tsv", sep="\t", index_col=0)
emptydrops_results.index = [sample + "-" + t[:-2] for t in emptydrops_results.index]

initial_clustering_results = pd.read_csv(f"{output_path}\\!cells_initial.csv", index_col=0)
initial_clustering_results["FDR"] = emptydrops_results.FDR
initial_clustering_results["is_retained_by_emptydrops"] = [t < FDR_CUTOFF for t in initial_clustering_results.FDR]

emptydrops_results = emptydrops_results[emptydrops_results.FDR < FDR_CUTOFF]
emptydrops_results["is_retained_by_iqc"] = [t in initial_clustering_results.index for t in emptydrops_results.index]

print("Cells retained by iqc:", len(initial_clustering_results.index))
print("Cells retained by emptydrops:", len(emptydrops_results.index))
print("iqc cells breakdown:\n", initial_clustering_results.is_retained_by_emptydrops.value_counts(), "\n")
print("emptydrops cells breakdown:\n", emptydrops_results.is_retained_by_iqc.value_counts())
