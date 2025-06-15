
#!/usr/bin/env python
# coding: utf-8

# %matplotlib inline
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib
import statsmodels.stats.multitest as smt
from scipy.stats import mannwhitneyu, ttest_ind
from os.path import join
from glob import glob
import warnings

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['figure.dpi'] = 100
warnings.filterwarnings("ignore")

# ==== Utility Functions ====

def load_sample_info(file_path):
    df = pd.read_csv(file_path, sep='\t')[[
        'sample_ID', 'rookery-region', 'Year', 'group-Lkadj',
        'Management units - based on DNA', 'Post_dupe_depth_2_all_LKGATK_USE',
        'Exclude-all_GL_ROH_analyses', 'Sequencing_source-type', 'Stage_Class'
    ]]
    return df.rename(columns={
        'sample_ID': 'sample',
        'rookery-region': 'rookery',
        'Year': 'year',
        'group-Lkadj': 'group',
        'Management units - based on DNA': 'MU',
        'Post_dupe_depth_2_all_LKGATK_USE': 'depth',
        'Exclude-all_GL_ROH_analyses': 'exclude',
        'Sequencing_source-type': 'batch',
        'Stage_Class': 'stage'
    })

def apply_color_palette(group_col):
    palettes = {
        'rookery': rookery_palette_dict,
        'MU': mu_palette_dict,
        'batch': batch_palette_dict,
        'group': population_palette_dict
    }
    return palettes.get(group_col, None)

def plot_box_swarm(
    df, x, y, hue=None, order=None, palette=None,
    figsize=(6, 4), title="", ylabel="", xlabel="", stat_test=False
):
    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(data=df, x=x, y=y, hue=hue, order=order, palette=palette, ax=ax, showfliers=False)
    sns.stripplot(data=df, x=x, y=y, hue=hue, order=order, palette=palette, ax=ax,
                  dodge=True if hue else False, linewidth=0.5, edgecolor='gray', alpha=0.7)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    # Statistical annotation
    if stat_test and order and len(order) == 2:
        group1 = df[df[x] == order[0]][y]
        group2 = df[df[x] == order[1]][y]
        stat, pval = mannwhitneyu(group1, group2)
        ax.text(0.5, 1.05, f"Mann–Whitney p = {pval:.2e}", transform=ax.transAxes,
                ha='center', va='bottom', fontsize=10)

    if hue:
        ax.legend_.remove()

    plt.tight_layout()
    return ax

# ==== Color Palettes ====

rookery_palette_dict = {
    'Mexico': '#1f77b4',
    'Atl Costa Rica': '#005a32',
    'Indonesia': '#d62728',
    'Pacific Costa Rica': '#7570b3',
    'Papua New Guinea': '#dd3497',
    'Solomon Islands': '#fde0dd',
    'Malaysia': '#2ca02c',
    'South Africa': '#17becf',
    'Ghana': '#ff7f0e',
    'Gabon': '#fdd0a2',
    'French Guiana-Suriname': '#8c564b',
    'Virgin Islands': '#d6eaf8',
    'Florida': '#d9f0a3',
    'na': '#000000',
    'unknown': '#7f7f7f'
}

mu_palette_dict = {
    'Eastern Pacific': '#1f77b4',
    'Western Pacific': '#d62728',
    'Northwest Caribbean': '#005a32',
    'Northeast Caribbean': '#6a51a3',
    'Northern Caribbean': '#d9f0a3',
    'South Africa': '#17becf',
    'SE Atlantic/West Africa': '#ff7f0e',
    'South-east Caribbean': '#8c564b',
    'Indo-Western Pacific': '#2ca02c',
    'na': '#000000',
    'unknown': '#7f7f7f'
}

batch_palette_dict = {
    'novo_UMass_Illumina_SR': '#1f77b4',
    'novo42_Illumina_SR': '#ff7f0e',
    'LR-getinfofromJH': '#2ca02c',
    'Duffy_Illumina_SR': '#d62728'
}

population_palette_dict = {
    'larger_declining': '#d62728',
    'small_stable':'#1f77b4',
    'small_declining': '#005a32',
    'larger_stable': '#ff7f0e',
    'small_recovering': '#6a51a3'
}

# ==== Data Loading ====

dir_path = '/Users/osipova/Documents/LabDocs/Leatherback_popgen/'
file_name = 'dc_222_samples_info.tsv'
file_data = glob(join(dir_path, file_name))[0]
INFO = load_sample_info(file_data)

# ==== You can now use the functions like this in your workflow ====
# plot_box_swarm(df=INFO, x='group', y='depth', order=['small_stable', 'larger_declining'], palette=apply_color_palette('group'))

