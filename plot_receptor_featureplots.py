#05: Loads neuron-only dataset and produces cleaned-up figures. 
import scanpy as sc
import matplotlib.pyplot as plt
import os

adata = sc.read_h5ad("SCG_Mapps_fullgenes_processed.h5ad")

genes = [
    "Th",
    "Dbh",
    "Esr1",
    "Pgr",
    "Gper1",
    "Oxtr",
    "Rxfp1"
]

os.makedirs("publication_umaps", exist_ok=True)

for gene in genes:
    sc.pl.umap(
        adata,
        color=gene,
        use_raw=True,
        cmap="viridis",
        size=18,
        frameon=False,
        show=False
    )

    plt.savefig(
        f"publication_umaps/{gene}.pdf",
        dpi=600,
        bbox_inches="tight"
    )

    plt.savefig(
        f"publication_umaps/{gene}.png",
        dpi=600,
        bbox_inches="tight"
    )

    plt.close()

print("Done!")

