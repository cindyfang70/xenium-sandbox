import scanpy as sc
import os,csv,re
import pandas as pd
import numpy as np
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
from io import StringIO
from scipy.io import mmread
import sys

adata = sc.read_10x_h5(
    filename="/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-09-06_Xen/output-XETG00089__0005434__Region_1__20230831__172339/cell_feature_matrix.h5"
)
df = pd.read_csv("/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2023-09-06_Xen/output-XETG00089__0005434__Region_1__20230831__172339/cells.csv.gz")


df.set_index(adata.obs_names, inplace=True)
adata.obs = df.copy()
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
adata.obs

br8667m_ids = pd.read_csv("./Br8667_Mid_5434_cell_ids.csv")
br8667m_ids["x"]
br8667m = adata[adata.obs["cell_id"].isin(br8667m_ids["x"])]

sys.setrecursionlimit(9999)

br8667m.var_names_make_unique()
spg.prefilter_specialgenes(br8667m) # TODO: need to remove negative control genes
sc.pp.normalize_total(br8667m)
sc.pp.log1p(br8667m)

x_array=br8667m.obs["x_centroid"].tolist()
y_array=br8667m.obs["y_centroid"].tolist()
spg_adj = spg.calculate_adj_matrix(x=x_array, y=y_array,histology=False)

print("adata shape")
print(br8667m.shape)
print("adjacency mat shape")
print(spg_adj.shape)

# set the hyperparameters
p=0.8 # they recommend setting a higher p for MERFISH-like data
l=spg.search_l(p, spg_adj, start=0.01, end=1000, tol=0.01, max_run=100)

n_clusters=7
#Set seed
r_seed=t_seed=n_seed=100
#Seaech for suitable resolution
res=spg.search_res(br8667m, spg_adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

