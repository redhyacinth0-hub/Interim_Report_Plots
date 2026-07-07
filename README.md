# Interim_Report_Plots
Data source:
GSE175421, SCG Mapps dataset from GEO.

Input files:
GSM5333151_SCG_rep1.dge.txt
GSM5333152_SCG_rep2.dge.txt
GSM5333153_SCG_rep3.dge.txt
GSM5333154_SCG_rep4.dge.txt
GSM5333155_SCG_rep5.dge.txt

Environment:
conda create -n sc python=3.10
conda activate sc
pip install scanpy anndata pandas numpy scipy matplotlib harmonypy igraph leidenalg

Run order:

build_scg_umap_fullgenes.py        

build_neuron_only_umaps.py        

plot_clusters.py                   

plot_receptor_umaps.py            

plot_receptor_featureplots.py      

