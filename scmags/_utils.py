import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad


def convert_form_anndata(adata: ad.AnnData, cell_annotation_col: str):
    """Function to retrieve the expression data, cell labels and gene names from an AnnData object.

    Args:
            adata (ad.AnnData): AnnData that contains .obs with the cell labels and .var with the gene names.
            cell_annotation_col (str): Column form adata.obs that contains the cell type annotation

    Returns:
            exp_data: expression data in the form of a numpy array (n_cells x n_genes)
            labels: 1D numpy array with the cell type labels
            gene_names: 1D numpy array with the gene names
    """

    adata.var_names_make_unique()
    exp_data = pd.DataFrame(
        data=adata.X.todense(), columns=adata.var_names, index=adata.obs.index
    ).to_numpy()
    labels = adata.obs[cell_annotation_col].to_numpy()
    gene_names = adata.var_names.to_numpy()

    return exp_data, labels, gene_names
