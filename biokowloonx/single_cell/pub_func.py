import os
import time
import umap
from typing import Union
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
from pathlib import Path
import scipy.stats as stats
from joblib import dump, load
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from anndata import AnnData, read_h5ad
from sklearn.metrics import mean_squared_error, r2_score
import gzip
import shutil


def print_df(df):
    assert type(df) == pd.DataFrame
    print('  >>  <<  ')
    print(df.shape)
    print(df.head(2))


def set_fig_style(font_family=None, font_size=None):
    fig, ax = plt.subplots()
    sns.set_style("white")
    try:
        # need to install the package of "SciencePlots" first, see https://github.com/garrettj403/SciencePlots
        plt.style.use(['science', 'no-latex'])
    except:
        print('No science style, please install the package of "SciencePlots" first, '
              'see https://github.com/garrettj403/SciencePlots')
        sns.set(palette='muted', font_scale=1.5)

    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'
    if font_family:
        mpl.rcParams['font.family'] = font_family
    if font_size:
        mpl.rcParams['font.size'] = font_size
    # print('figure.dpi will be set to', mpl.rcParams['figure.dpi'])
    plt.close('all')


def log2_transform(df) -> pd.DataFrame:
    """
    log2 transform expression values (plus 1)

    :param df: a dataframe with expression values, cells by genes

    :return: a dataframe with log2 transformed values
    """
    df = df.astype(np.float64)
    df = np.log2(df + 1)
    df = df.astype(np.float32)
    return df


def center_value(df, return_mean=False):
    """
    exp - mean(exp) for each genes
    :param df: expression dataframe, gene x sample
    :param return_mean: if return df_mean
    :return:
    """
    df_mean = df.mean(axis=1)
    if return_mean:
        return df - np.vstack(df_mean), df_mean.to_frame('gene_mean')
    return df - np.vstack(df_mean)


def cal_relative_error(y_true, y_pred, max_error=None, min_error=None):
    """
    calculate relative error between two dataFrame
    relative error = (y_pred - y_true) / y_true
    :param y_true: dataframe
    :param y_pred: dataframe, two dataframe have same index and columns (also same order)
    :param max_error: float
        all errors should <= this value
    :param min_error: float
        all errors should >= this value
    :return:
    """
    relative_error = (y_pred - y_true) / y_true
    if max_error:
        relative_error[relative_error > max_error] = max_error
    if min_error:
        relative_error[relative_error < min_error] = min_error
    return relative_error


def calculate_rmse(y_true: pd.DataFrame, y_pred: pd.DataFrame):
    """
    https://scikit-learn.org/stable/modules/generated/sklearn.metrics.mean_squared_error.html
    calculate the RMSE of each cell types by columns
    :param y_true: a dataFrame
        shape: number of samples x number of cell types
    :param y_pred: a dataFrame
    :return:
    """
    if y_true.shape[1] == 1:  # only one feature
        multioutput = 'uniform_average'
    else:  # multiple cell types
        multioutput = 'raw_values'
    return mean_squared_error(y_true=y_true, y_pred=y_pred, multioutput=multioutput, squared=False)


def calculate_r2(y_true, y_pred):
    """
    https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html
    calculate the R^2 (coefficient of determination) of each cell types by columns
    :param y_true:
    :param y_pred:
    :return:
    """
    return r2_score(y_true=y_true, y_pred=y_pred, multioutput='raw_values')


def check_dir(path):
    """
    check if a path exist, create if not exist
    :param path:
    :return:
    """
    if not os.path.exists(path):
        os.makedirs(path)


def parse_log_file(log_file_path, search_type='sample'):
    """
    parse log file, @sample xxx, ...
    :param log_file_path:
    :param search_type: search gene or sample in log file
    :return: list
        a list of sample names
    """
    sample_names = []
    if os.path.exists(log_file_path):
        with open(log_file_path, 'r') as f:
            for i in f:
                i = i.strip()
                _sample_name = i.split(',')[0].replace(f'@{search_type} ', '')
                sample_names.append(_sample_name)
    return sample_names


def write_to_log(log_file_path, template, sample_names, n_round):
    """

    :param log_file_path:
    :param template: '@sample {}, removed, ..., round {}', only two {} for sample name and n_round
    :param sample_names:
    :param n_round:
    :return:
    """
    sample_names_in_log_file = parse_log_file(log_file_path)
    with open(log_file_path, 'a') as f:
        for s in sample_names:
            _info = template.format(s, n_round)
            if s not in sample_names_in_log_file:
                f.write(_info + '\n')
            print(_info)


def extract_gz_file(file_dir: str, file_name: str):
    file_path = os.path.join(file_dir, file_name)
    file_path_out = os.path.join(file_dir, file_name.replace('.gz', ''))
    if os.path.exists(file_path) and file_name.endswith('.gz'):
        if not os.path.exists(file_path_out):
            with gzip.open(file_path, 'rb') as f_in:
                with open(file_path_out, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)


def create_h5ad_dataset(simulated_bulk_exp_file_path, cell_fraction_file_path,
                        dataset_info, result_file_path: str,
                        filtering: bool = False, merge_t_cell: bool = False, gep_type='bulk'):
    """
    create .h5ad file according to cell fraction and simulated bulk expression profiles
    https://anndata.readthedocs.io/en/latest/index.html
    :param simulated_bulk_exp_file_path: simulated bulk expression profile, samples by genes
        .csv file or .h5ad file
    :param cell_fraction_file_path: .csv file, samples by cell types
    :param dataset_info: str
    :param result_file_path:
    :param filtering: if filtered by marker ratio
    :param merge_t_cell: whether merge the cell fraction of CD8 T and CD4 T cells
    :param gep_type: bulk (mixture contains >= 2 cell types) / sct (single cell type)
    :return:
    """
    try:
        if filtering:
            simulated_bulk_exp_file_path = simulated_bulk_exp_file_path.replace('_log2cpm1p.csv',
                                                                                '_log2cpm1p_filtered.csv')
            cell_fraction_file_path = cell_fraction_file_path.replace('.csv', '_filtered.csv')
        simu_bulk_exp_raw = read_df(simulated_bulk_exp_file_path)
    except UnicodeDecodeError:
        simu_bulk_exp_raw = read_h5ad(simulated_bulk_exp_file_path)
    # simu_bulk_exp.index = simu_bulk_exp.index.astype(str)
    cell_frac = read_df(cell_fraction_file_path)
    if merge_t_cell:
        if 'T Cells' not in cell_frac.columns:
            cell_frac['T Cells'] = cell_frac.loc[:, ['CD4 T', 'CD8 T']].sum(axis=1)
        cell_frac.drop(columns=['CD4 T', 'CD8 T'], inplace=True)
    sample_list = cell_frac.index.to_list()
    # cell_frac.index = cell_frac.index.astype(str)
    uns = {'cell_types': cell_frac.columns.to_list(),
           'dataset_info': dataset_info}
    simu_bulk_exp = pd.DataFrame()
    if type(simu_bulk_exp_raw) == pd.DataFrame:
        simu_bulk_exp = simu_bulk_exp_raw.loc[sample_list, :]
    elif type(simu_bulk_exp_raw) == AnnData:
        simu_bulk_exp_df = pd.DataFrame(simu_bulk_exp_raw.X, index=simu_bulk_exp_raw.obs.index,
                                        columns=simu_bulk_exp_raw.var.index)
        simu_bulk_exp_df = simu_bulk_exp_df.loc[sample_list, :]
        simu_bulk_exp = simu_bulk_exp_df.copy()
    if (cell_frac.shape[0] == simu_bulk_exp.shape[0]) and not np.all(simu_bulk_exp.index == cell_frac.index):
        cell_frac = cell_frac.loc[simu_bulk_exp.index, :].copy()
    elif cell_frac.shape[0] != simu_bulk_exp.shape[0]:
        print(f'   The shape of simu_bulk_exp: {simu_bulk_exp.shape}, the shape of cell_frac: {cell_frac.shape}')
        cell_frac = cell_frac[~cell_frac.index.duplicated(keep='first')].copy()
        simu_bulk_exp = simu_bulk_exp[~simu_bulk_exp.index.duplicated(keep='first')].copy()
        intersection_inx = [i for i in simu_bulk_exp.index if i in cell_frac.index]
        print(f'   The length of intersections in both cell_frac and simu_bulk_exp: {len(intersection_inx)}')
        cell_frac = cell_frac.loc[intersection_inx, :].copy()
        simu_bulk_exp = simu_bulk_exp.loc[intersection_inx, :].copy()
    # print(cell_frac.index)
    # print(simu_bulk_exp.index)
    if np.all(simu_bulk_exp.index == cell_frac.index):
        var = pd.DataFrame(index=simu_bulk_exp.columns, columns=[f'in_{gep_type}'])
        var[f'in_{gep_type}'] = 1
        adata = AnnData(X=simu_bulk_exp.values.astype('float32'), obs=cell_frac, uns=uns,
                        var=var, dtype=np.dtype('float32'))
        adata.write_h5ad(filename=Path(result_file_path), compression='gzip')
    else:
        raise KeyError('simu_bulk_exp and cell_frac file should have same sample order')


def log_exp2cpm(exp_df: Union[pd.DataFrame, np.array], log_base=2, correct=1) -> Union[pd.DataFrame, np.array]:
    """
    Convert log2(CPM + 1) to non-log space values (CPM / TPM)

    :param exp_df: samples by genes

    :param log_base: the base of log transform

    :param correct: plus 1 for avoiding log transform 0

    :return: counts per million (CPM) or transcript per million (TPM)
    """
    exp = np.power(log_base, exp_df) - correct
    # exp = exp.astype(np.float64)
    cpm = exp / np.vstack(exp.sum(axis=1)) * 1e6
    # cpm = cpm.astype(np.float32)
    return cpm


def non_log2log_cpm(input_file_path: Union[str, pd.DataFrame], result_file_path: str = None,
                    transpose: bool = True, correct: int = 1):
    """
    Convert non-log expression data to log2(CPM + 1) or log2(TPM + 1)

    :param input_file_path: non-log space expression file, genes by samples

    :param result_file_path: file path, samples by genes

    :param transpose: if input file is samples by genes, set to False, otherwise set to True

    :param correct: plus 1 for avoiding log transform 0

    :return: log2(CPM + 1) or save result to file, samples by genes if transpose is True, otherwise genes by samples
    """

    bulk_exp = pd.DataFrame()
    if type(input_file_path) is str:
        sep = get_sep(input_file_path)
        bulk_exp = pd.read_csv(input_file_path, index_col=0, sep=sep)
    elif type(input_file_path) is pd.DataFrame:
        bulk_exp = input_file_path
    if transpose:
        bulk_exp = bulk_exp.T  # transpose to samples by genes
    bulk_exp = non_log2cpm(bulk_exp)  # CPM/TPM
    bulk_exp = np.log2(bulk_exp + correct)
    if result_file_path is not None:
        bulk_exp.round(3).to_csv(result_file_path)
    else:
        return bulk_exp.round(3)


def non_log2cpm(exp_df, sum_exp=1e6) -> pd.DataFrame:
    """
    Normalize gene expression to CPM / TPM for non-log space

    :param exp_df: gene expression profile in non-log space, sample by gene

    :param sum_exp: sum of gene expression for each sample, default is 1e6

    :return: counts per million (CPM) or transcript per million (TPM)
    """
    return exp_df / np.vstack(exp_df.sum(axis=1)) * sum_exp


def get_corr(df_col1, df_col2, return_p_value=False) -> Union[float, tuple]:
    """
    calculate the correlation between two columns of dataframe
    :param df_col1: series, column1
    :param df_col2: series, column2
    :param return_p_value: if return p-value
    :return:
    """
    # correlation = np.corrcoef(df_col1, df_col2)
    corr, p_value = stats.pearsonr(df_col1, df_col2)
    if return_p_value:
        return corr, p_value
    else:
        return corr
    # return correlation[0, 1]


def get_sep(file_path, comment: str = None):
    """
    check the separater (`\t` or `,`) in this file
    :param file_path:
    :param comment: if remove comment lines start with ''
    :return:
    """
    sep = '\t'
    with open(file_path, 'r') as f_handle:
        # first_line = ''
        while True:
            first_line = f_handle.readline()
            if comment is not None:
                if first_line.startswith(comment):
                    continue
                else:
                    break
            break
        if ',' in first_line:
            sep = ','
    return sep


def cal_exp_by_gene_list(exp_df, gene_list, min_exp_value=0, method='mean'):
    """
    calculate mean expression (or max) of a gene list (for single cell type)
    :param exp_df: sample by gene, usually in CPM / TPM
    :param gene_list: a list of genes which all are included in exp_df
    :param min_exp_value: min mean expression of the gene list
    :param method: mean or max
    :return: mean or max expression value of marker genes for each sample
    """
    for gene in gene_list:
        if gene not in exp_df.columns:
            raise KeyError(f'Gene {gene} not include in exp_df')
    current_exp = exp_df.loc[:, gene_list].copy()
    if method == 'mean':
        current_value = current_exp.mean(axis=1)
    elif method == 'max':
        current_value = current_exp.max(axis=1)
    else:
        raise KeyError('Only "mean" or "max" allowed')
    current_value[current_value < min_exp_value] = min_exp_value
    if exp_df.shape[0] == 1:
        return round(float(current_value), 3)  # float
    return current_value.round(3)  # pd.Series


def print_msg(p_str, log_file_path=None):
    print()
    current_info = f'---->>> {p_str} <<<----'
    current_time = time.ctime()
    print(current_info)
    print(current_time)
    if log_file_path is not None:
        with open(log_file_path, 'a') as f_handle:
            f_handle.write(current_info + '\n')
            f_handle.write(current_time + '\n')
            f_handle.write('\n')


def read_df(df_file: Union[str, pd.DataFrame, np.ndarray]) -> Union[pd.DataFrame, np.ndarray]:
    """
    check the type of df_file
    - if df_file is a file path, read this file
    - if df_file is a DataFrame, return directly
    """
    if type(df_file) == str:  # the file path of current file
        sep = get_sep(df_file)  # separated by '\t' or ','
        df = pd.read_csv(df_file, index_col=0, sep=sep)
    elif type(df_file) == pd.DataFrame:
        return df_file  # return DataFrame directly
    elif type(df_file) == np.ndarray:
        return df_file  # return np.ndarray directly
    else:
        raise TypeError(
            f'Only file path or pd.DataFrame was supported by df_file, {type(df_file)} is not supported.')
    return df


def do_pca_analysis(exp_df, n_components=5, pca_result_fp=None, save_model: bool = False):
    """
    PCA analysis
    :param exp_df:
    :param n_components:
    :param pca_result_fp:
    :param save_model:
    :return: fitted PCA model
    """
    if os.path.exists(pca_result_fp):
        print(f'Loading PCA result from file: {pca_result_fp}')
        pca = load(pca_result_fp)
    else:
        pca = PCA(n_components=n_components)
        pca.fit(exp_df)
        if save_model:
            dump(pca, pca_result_fp)
    return pca


def do_umap_analysis(exp_df, n_components=5, n_neighbors=15, min_dist=0.1,
                     umap_model_result_fp=None, save_model: bool = False):
    """
    t-SNE analysis
    :param exp_df:
    :param n_components:
    :param n_neighbors:
    :param min_dist:
    :param umap_model_result_fp:
    :param save_model:
    :return:
    """
    if os.path.exists(umap_model_result_fp):
        print(f'Loading UMAP result from file: {umap_model_result_fp}')
        umap_model = load(umap_model_result_fp)
    else:
        umap_model = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components)
        umap_model.fit(exp_df)
        if save_model:
            dump(umap_model, umap_model_result_fp)
    return umap_model


def get_ccc(x, y):
    # Concordance Correlation Coefficient(CCC), https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    vx, cov_xy, cov_xy, vy = np.cov(x, y, bias=True).flatten()
    mx, my = x.mean(), y.mean()
    return 2*cov_xy / (vx + vy + (mx-my)**2)


if __name__ == '__main__':
    pass
