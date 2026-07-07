#01
import scanpy as sc
import anndata as ad
import scipy.sparse as sp

files = [
    "GSM5333151_SCG_rep1.dge.txt",
    "GSM5333152_SCG_rep2.dge.txt",
    "GSM5333153_SCG_rep3.dge.txt",
    "GSM5333154_SCG_rep4.dge.txt",
    "GSM5333155_SCG_rep5.dge.txt",
]

adatas = []

for f in files:
    print("Loading", f)
    a = sc.read_csv(f, delimiter="\t").transpose()
    a.obs["sample"] = f.replace(".dge.txt", "")

    sc.pp.filter_cells(a, min_counts=1500)
    sc.pp.filter_cells(a, max_counts=45000)

    a.X = sp.csr_matrix(a.X)
    print("After cell filter:", a)
    adatas.append(a)

adata = ad.concat(adatas, join="inner")
print("Combined:", adata)

adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)
adata = adata[adata.obs["pct_counts_mt"] < 10].copy()

# keep genes expressed in at least 1 cell
sc.pp.filter_genes(adata, min_cells=1)

print("After QC:", adata)

# normalize full gene matrix
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

# save full expression for plotting
adata.raw = adata.copy()

# compute UMAP using HVGs only, but keep all genes in adata.raw
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50)

adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

adata.write("SCG_Mapps_fullgenes_processed.h5ad")
print("Saved SCG_Mapps_fullgenes_processed.h5ad")
