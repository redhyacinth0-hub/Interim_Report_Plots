#02: Loads the data from step 1 and selects only sympathetic neurons. Recomputes neuron-only UMAP. 
import scanpy as sc
import matplotlib.pyplot as plt
import os

adata = sc.read_h5ad("SCG_Mapps_fullgenes_processed.h5ad")

# TH/DBH-positive sympathetic neuron clusters from your Leiden plot
neuron_clusters = ["1", "2", "3", "4"]

neurons = adata[adata.obs["leiden"].isin(neuron_clusters)].copy()
print("Neuron-only object:", neurons)

# Recompute neuron-only UMAP
sc.pp.highly_variable_genes(neurons, min_mean=0.0125, max_mean=3, min_disp=0.5)
neurons_hvg = neurons[:, neurons.var["highly_variable"]].copy()

sc.pp.scale(neurons_hvg, max_value=10)
sc.tl.pca(neurons_hvg, n_comps=30)

neurons.obsm["X_pca"] = neurons_hvg.obsm["X_pca"]

sc.pp.neighbors(neurons, n_neighbors=20, n_pcs=20)
sc.tl.umap(neurons)
sc.tl.leiden(neurons, resolution=0.4)

neurons.write("SCG_Mapps_neurons_only.h5ad")

genes = ["Th", "Dbh", "Esr1", "Pgr", "Gper1", "Oxtr", "Rxfp1"]

os.makedirs("neuron_only_umaps", exist_ok=True)

for gene in genes:
    sc.pl.umap(
        neurons,
        color=gene,
        use_raw=True,
        cmap="viridis",
        size=25,
        frameon=False,
        vmax="p99",
        show=False,
    )

    ax = plt.gca()
    ax.set_title(gene, fontstyle="italic", fontsize=24)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig(f"neuron_only_umaps/{gene}.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"neuron_only_umaps/{gene}.pdf", bbox_inches="tight")
    plt.close()

print("Done. Saved neuron-only figures in neuron_only_umaps/")

