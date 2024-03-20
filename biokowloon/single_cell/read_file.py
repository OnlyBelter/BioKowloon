import os
import numpy as np
import pandas as pd
import anndata as an
from typing import Union
from scipy.sparse import csr_matrix
from sklearn import preprocessing as pp
from .pub_func import (log_exp2cpm, read_df, non_log2log_cpm, non_log2cpm)


class ReadH5AD(object):
    """
    Read .h5ad file, usually the values are log2 transformed

    :param file_path: the file path of .h5ad file, samples by genes, log2cpm1p format
    :param show_info: whether to show the information of the dataset after reading
    """
    def __init__(self, file_path: str, show_info: bool = False):
        """
        """
        self.dataset = an.read_h5ad(file_path)
        if show_info:
            print(self.dataset)

    def get_df(self, result_file_path: str = None, convert_to_tpm: bool = False,
               scaling_by_sample: bool = False) -> pd.DataFrame:
        """
        Convert to DataFrame, samples by genes, log space (log2cpm1p)

        :param result_file_path:
        :param convert_to_tpm: whether to convert log2cpm1p to TPM
        :param scaling_by_sample: whether to scale the expression values of each sample to [0, 1] by 'min_max'
        """
        if type(self.dataset.X) == csr_matrix:
            x_data = self.dataset.X.A.astype(np.float32)  # convert sparse matrix to dense matrix
        else:
            x_data = self.dataset.X.astype(np.float32)

        if convert_to_tpm:
            x_data = log_exp2cpm(x_data)
        if scaling_by_sample:
            scaler = pp.MinMaxScaler(feature_range=(0, 1), copy=True)
            x_data = scaler.fit_transform(x_data.T).T

        df = pd.DataFrame(data=x_data, index=self.dataset.obs.index,
                          columns=self.dataset.var.index).round(3)
        if result_file_path is not None:
            df.to_csv(result_file_path, float_format='%.3f')
        return df

    def get_cell_fraction(self) -> Union[None, pd.DataFrame]:
        """
        Get cell fraction, cells by cell types
        """
        if self.dataset.obs.shape[1] > 0:
            return self.dataset.obs.round(3)
        else:
            print('   There is no cell fraction in this .h5ad file')
            return None

    def get_h5ad(self):
        """
        Get the .h5ad file
        """
        return self.dataset


class ReadExp(object):
    """
    Read gene expression file, and convert to specific format (TPM / CPM, log2cpm1p)

    - TPM: transcript per million

    - CPM: UMI reads per million (3' end sc-RNA seq), same as TPM in the full-length RNA-seq of bulk cells

    - log_space: log2(CPM + 1), or log2(TPM + 1)

    - non_log: non log space, could be normalized to TPM

    - Data from full-length protocols may benefit from normalization methods that take into account gene length
      (e.g. Patelet al, 2014; Kowalczyket al,2015; Soneson & Robinson, 2018), while 3' enrichment data do not.

    - A commonly used normalization method for full-length scRNA-seq data is TPM normalization (Liet al, 2009),
      which comes from bulk RNA-seq analysis. (Luecken, M. D. & Theis, F. J., Mol. Syst. Biol. 15, e8746 (2019))

    :param exp_file: file path or DataFrame, samples by genes
    :param exp_type: TPM / CPM, log_space, non_log
    :param transpose: transpose if exp_file formed as genes (index) by samples (columns)
    """
    def __init__(self, exp_file, exp_type='TPM', transpose: bool = False):
        """
        """
        assert exp_type in ['TPM', 'CPM', 'log_space', 'non_log']
        self.file_type = exp_type
        self.exp = read_df(exp_file)
        if transpose:
            self.exp = self.exp.T
        self.scaled_by_sample = False

    def to_tpm(self):
        """
        Convert to TPM
        """
        if self.file_type == 'non_log':
            self.exp = non_log2cpm(self.exp)
        elif self.file_type == 'TPM' or self.file_type == 'CPM':
            pass
        elif self.file_type == 'log_space':
            self.exp = log_exp2cpm(self.exp)
        self.file_type = 'TPM'

    def to_log2cpm1p(self):
        """
        Convert to log2(TPM + 1)
        """
        if self.file_type != 'log_space':
            self.exp = non_log2log_cpm(self.exp, transpose=False)
        else:
            print('   This file has already log2 transformed.')
        self.file_type = 'log_space'

    def get_file_type(self) -> str:
        """
        Get the file type
        """
        return self.file_type

    def get_exp(self) -> pd.DataFrame:
        """
        Get the expression matrix
        """
        return self.exp.round(3)

    def save(self, file_path, sep=',', transpose: bool = False):
        """
        Save the expression matrix to file

        :param file_path: file path
        :param sep: separator, default is ','
        :param transpose: transpose index and columns
        """
        if transpose:
            self.exp = self.exp.T.copy()
        self.exp.to_csv(file_path, sep=sep, float_format='%.3f')

    def do_scaling(self):
        """
        Scaling GEPs by sample to [0, 1], same as Scaden
        """
        if not self.scaled_by_sample:
            scaler = pp.MinMaxScaler(feature_range=(0, 1), copy=True)
            x_scaled = scaler.fit_transform(self.exp.T).T  # scaling by column (sample), so T is needed here
            self.scaled_by_sample = True
            self.exp = pd.DataFrame(data=x_scaled, index=self.exp.index,
                                    columns=self.exp.columns).round(3)
        # return self.exp

    def do_scaling_by_constant(self, divide_by=20):
        """
        Scaling GEPs by dividing a constant in log space, so all expression values are in [0, 1)
        """
        if self.file_type != 'log_space':
            raise ValueError('   This file is not in log space')
        if np.any(self.exp.values > 1.0):
            self.exp = self.exp / divide_by

    def align_with_gene_list(self, gene_list: list = None, fill_not_exist=False, pathway_list: bool = False):
        """
        Align the expression matrix with a gene list and rescale to TPM or log2(TPM + 1)

        :param gene_list: gene list
        :param fill_not_exist: fill 0 if gene not exist in the provided gene_list when True
        :param pathway_list: gene list contains pathway names, so TPM normalization is not suitable
        """
        common_genes = [i for i in gene_list if i in self.exp.columns]
        not_exist_in_gene_list = [i for i in gene_list if i not in common_genes]
        removed_genes = [i for i in self.exp.columns if i not in common_genes]
        print(f'   {len(common_genes)} common genes will be used, {len(removed_genes)} genes will be removed.')
        self.exp = self.exp.loc[:, common_genes].copy()
        if fill_not_exist and (len(not_exist_in_gene_list) != 0):
            print(f'   {len(not_exist_in_gene_list)} genes are not in current dataset, 0 will be filled')
            _not_exist_exp = pd.DataFrame(np.zeros((self.exp.shape[0], len(not_exist_in_gene_list))), index=self.exp.index,
                                          columns=not_exist_in_gene_list)
            self.exp = pd.concat([self.exp, _not_exist_exp], axis=1)
            self.exp = self.exp.loc[:, gene_list].copy()
        if not pathway_list:
            if self.file_type == 'log_space':  # scaling to TPM after alignment
                self.to_tpm()
                self.to_log2cpm1p()
            else:
                self.file_type = 'non_log'
                self.to_tpm()


def read_gene_set(gene_set_file_path: list, max_n_genes: int = 300) -> pd.DataFrame:
    """
    read gene set from .gmt files and convert to DataFrame with genes as index and gene sets as columns,
     1 for a gene in a gene set and 0 for not
    :param gene_set_file_path: the file path of gene set
    :param max_n_genes: the maximum number of genes in a gene set
    :return: DataFrame of gene set with genes as index and gene sets as columns
    """
    gs2genes = {}
    all_genes = set()
    for gs_file in gene_set_file_path:
        if not os.path.exists(gs_file):
            raise FileNotFoundError(f'gene set file {gs_file} not found')
        with open(gs_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    gs, _, *genes = line.split('\t')
                    if len(genes) > max_n_genes:
                        genes = genes[:max_n_genes]
                    gs2genes[gs] = genes
                    all_genes.update(genes)
    gene_set_df = pd.DataFrame(index=list(all_genes), columns=list(gs2genes.keys()))
    for gs, genes in gs2genes.items():
        gene_set_df.loc[genes, gs] = 1
    gene_set_df.fillna(0, inplace=True)
    return gene_set_df


def read_cancer_purity(cancer_purity_file_path, sample_names: list):
    """
    Tumor purity estimates for TCGA samples
    Aran, D., Sirota, M. & Butte, A. Systematic pan-cancer analysis of tumour purity. Nat Commun 6, 8971 (2015).
    https://doi.org/10.1038/ncomms9971
    :param cancer_purity_file_path:
    :param sample_names: all sample names need to compare
    :return:
    """
    cancer_purity = pd.read_csv(cancer_purity_file_path, index_col=0)
    cancer_purity = cancer_purity.loc[~cancer_purity['CPE'].isnull()].copy()
    # print_df(cancer_purity)
    sample_name_mapping = {i[0:16]: i for i in sample_names}
    common_sample = list(set(list(sample_name_mapping.keys())) & set(cancer_purity.index))
    cancer_purity = cancer_purity.loc[common_sample, :].copy()
    cancer_purity.index = cancer_purity.index.map(sample_name_mapping)
    return cancer_purity


def read_marker_gene(marker_gene_file_path: str, include_t_cell: bool = False,
                     include_cd8_nk_marker: bool = False, use_cancer_cell: bool = False,
                     add_top_corr_gene: bool = False, corr_mean: bool = False) -> dict:
    """
    read marker genes for each cell type
    :param marker_gene_file_path: file path of selected marker genes for each cell type
    :param include_t_cell: if include the marker genes of T Cells for both CD4 and CD8 T
    :param include_cd8_nk_marker: if remain "CD8 T / NK" marker genes
    :param use_cancer_cell: if use "Cancer Cells" to replace "Epithelial Cells"
    :param add_top_corr_gene: add top correlated genes from the corr between cell fraction and gene expression value
    :param corr_mean: only two genes for CD4 T Cells and two genes for B Cells, others are same as core_marker
    :return: a dict of cell type to marker genes, {'': []}
    """
    cell_type2marker = {}
    marker_gene = pd.read_csv(marker_gene_file_path)
    if 'corr_mean' in marker_gene.columns and corr_mean:
        marker_gene = marker_gene.loc[marker_gene['corr_mean'] == 1, :].copy()
    if 'core_marker' in marker_gene.columns:
        if add_top_corr_gene:
            marker_gene = marker_gene.loc[marker_gene['core_marker'].isin([1, 2]), :].copy()
        else:
            marker_gene = marker_gene.loc[marker_gene['core_marker'] == 1, :].copy()
        if not include_cd8_nk_marker:
            marker_gene = marker_gene.loc[marker_gene['cell_type'] != 'CD8 T / NK', :].copy()
    cell_types = sorted(marker_gene['cell_type'].unique())
    for ct in cell_types:
        if include_t_cell:
            cell_type2marker[ct] = list(marker_gene.loc[marker_gene['cell_type'] == ct, 'marker_gene'].unique())
        else:
            if ct != 'T Cells':
                cell_type2marker[ct] = list(marker_gene.loc[marker_gene['cell_type'] == ct, 'marker_gene'].unique())
    if include_cd8_nk_marker and cell_type2marker.get('CD8 T / NK', ''):
        cell_type2marker['NK'] += cell_type2marker['CD8 T / NK']
        cell_type2marker['CD8 T'] += cell_type2marker['CD8 T / NK']
    if use_cancer_cell and 'Epithelial Cells' in cell_type2marker:
        cell_type2marker['Cancer Cells'] = cell_type2marker['Epithelial Cells']
        del cell_type2marker['Epithelial Cells']
    return cell_type2marker
