#04: Loads neuron-only dataset and produces receptor UMAPs using viridis (default Scanpy style). 
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

adata = sc.read_h5ad("SCG_Mapps_light_processed.h5ad")

genes = ["Esr1", "Esr2", "Pgr", "Pgrmc1", "Pgrmc2", "Gper1", "Rxfp1", "Oxtr", "Prlr", "Th"]
Path("receptor_umaps").mkdir(exist_ok=True)

for gene in genes:
    if gene not in adata.raw.var_names:
        print(f"{gene} not found")
        continue

    sc.pl.umap(
        adata,
        color=gene,
        use_raw=True,
        cmap="viridis",
        size=12,
        frameon=False,
        vmax="p99",
        show=False,
    )

    ax = plt.gca()
    ax.set_title(gene, fontstyle="italic", fontsize=22)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig(f"receptor_umaps/{gene}_UMAP.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"receptor_umaps/{gene}_UMAP.pdf", bbox_inches="tight")
    plt.close()

print("Done. Figures saved in receptor_umaps/")
