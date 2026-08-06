#03: Loads the neuron-only dataset and produces Leiden cluster figures. 
import scanpy as sc

adata = sc.read_h5ad("SCG_Mapps_neurons_only.h5ad")

print("Loaded neuron-only object:")
print(adata)
print("Number of cells:", adata.n_obs)

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    frameon=False,
    save="_neuron_clusters.png"
)
