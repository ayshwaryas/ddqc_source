import numpy as np
import pandas as pd

INF = 10 ** 10  # infinity for cases no filtering is required


def mad(data, axis=0, constant=1.4826):
    return constant * np.median(np.absolute(data - np.median(data, axis)), axis)


def metric_filter(data, classification, method, param, metric_name, do_lower_co=False, do_upper_co=False,
                  lower_bound=INF, upper_bound=-INF):
    data[metric_name + "_qc_pass"] = False  # T/F array to tell whether the cell is filtered
    data[metric_name + "_lower_co"] = None  # array recording lower cutoff for cell (if exists)
    data[metric_name + "_upper_co"] = None  # array recording upper cutoff for cell (if exists)
    if method == "mad":
        data[metric_name + "_median"] = None
        data[metric_name + "_mad"] = None
    for ct in data[classification].unique():  # iterate though all cell_types
        lower_co = -INF
        upper_co = INF
        cell_type = data[data[classification] == ct]  # subset adata based on cell_type number

        if method == "mad":  # calculate MAD cutoffs, which are median ± param * MAD
            data.loc[data[classification] == ct, metric_name + "_median"] = np.median(cell_type[metric_name])
            data.loc[data[classification] == ct, metric_name + "_mad"] = mad(cell_type[metric_name])
            if do_lower_co:
                lower_co = min(np.median(cell_type[metric_name]) - param * mad(cell_type[metric_name]), lower_bound)
            if do_upper_co:
                upper_co = max(np.median(cell_type[metric_name]) + param * mad(cell_type[metric_name]), upper_bound)

        if method == "outlier":  # calculate Outlier cutoffs, which are Q1 - 1.5 * IQR or Q3 + 1.5 * IQR
            q75, q25 = np.percentile(cell_type[metric_name], [75, 25])
            if do_lower_co:
                lower_co = min(q25 - 1.5 * (q75 - q25), lower_bound)
            if do_upper_co:
                upper_co = max(q75 + 1.5 * (q75 - q25), upper_bound)

        filters = [
            data[classification] == ct,
            data[metric_name] >= lower_co,
            data[metric_name] <= upper_co
        ]  # filtering condition
        if do_upper_co:
            data.loc[data[classification] == ct, metric_name + "_upper_co"] = upper_co
        if do_lower_co:
            data.loc[data[classification] == ct, metric_name + "_lower_co"] = lower_co
        # for cells that satisfy the condition set the value to true
        data.loc[np.logical_and.reduce(filters), metric_name + "_qc_pass"] = True
    return data


def filter_cells(data, classification, method, threshold, n_genes_lower_bound=200, percent_mito_upper_bound=10,
                 do_counts=True, do_genes=True, do_mito=True, do_ribo=False, return_full=False):
    data_copy = data.copy()

    if do_counts:
        data_copy = metric_filter(data_copy, classification, method, threshold, "n_counts", do_lower_co=True)
    else:
        data_copy["n_counts_qc_pass"] = True
    if do_genes:
        data_copy = metric_filter(data_copy, classification, method, threshold, "n_genes", do_lower_co=True,
                                  lower_bound=n_genes_lower_bound)
    else:
        data_copy["n_genes_qc_pass"] = True
    if do_mito:
        data_copy = metric_filter(data_copy, classification, method, threshold, "percent_mito", do_upper_co=True,
                                  upper_bound=percent_mito_upper_bound)
    else:
        data_copy["percent_mito_qc_pass"] = True
    if do_ribo:
        data_copy = metric_filter(data_copy, classification, method, threshold, "percent_ribo", do_upper_co=True)
    else:
        data_copy["percent_ribo_qc_pass"] = True

    filters = [
        data_copy["n_counts_qc_pass"],
        data_copy["n_genes_qc_pass"],
        data_copy["percent_mito_qc_pass"],
        data_copy["percent_ribo_qc_pass"],
    ]  # filtering condition
    data_copy["passed_qc"] = False  # cumulative T/F array for to determine cells that passed the filtering
    # for cells that satisfy the condition set the value to true
    data_copy.loc[np.logical_and.reduce(filters), "passed_qc"] = True

    if not return_full:
        data[classification + "_passed_qc"] = data_copy.passed_qc  # transfer array from the copy to actual object
        return data
    else:
        data_copy[classification + "_passed_qc"] = data_copy.passed_qc  # transfer array from the copy to actual object
        return data_copy


if __name__ == '__main__':
    path = "Z:\\revisions\\revisions2\\cell_classification\\"
    tissue = "krasnow_lung"
    classification_data = pd.read_csv(f"{path}{tissue}\\summary.csv", index_col=0)
    for c in ["cell_typist", "ddqc_cluster"]:  # "single_r", "azimuth",
        classification_data = filter_cells(classification_data, classification=c, method="mad", threshold=2,
                                           do_counts=False)
    with open(f"{path}{tissue}\\classification_mad_results.csv", "w") as file:
        file.write(classification_data.to_csv())
