import scanpy as sc

adata = sc.read_h5ad("SCG_Mapps_fullgenes_processed.h5ad")

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    frameon=False,
    save="_clusters.png"
)
