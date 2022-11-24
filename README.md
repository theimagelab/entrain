# Entrain

Entrain is an R package made for finding environmental signals that regulate a cell differentiation trajectory. 
If you have a dataset of cells that lie along a differentiation trajectory, and/or exhibit robust RNA velocity dynamics, then Entrain will help you generate hypotheses for identifying the environmental signals that are driving your observed dynamics. 

While existing ligand-receptor tools exist, these tools (e.g., TraSig, CellChat) are reliant solely on the expression level of ligand and receptor genes,  ignoring regulons and marker genes of a trajectory which constitute a large proportion of transcriptional information. Here, we score ligands based on their contribution towards regulating genes implicated in driving pseudotime-based trajectories or genes with significant RNA velocity vectors.

Entrain is based on, and can be considered a unification, of existing tools: Monocle3, NicheNet and RNA Velocity. Furthermore, Entrain works directly with the native single-cell data types, namely SeuratObject and anndata.

## Documentation
Full documentation, including example analysis, can be found at https://theimagelab.github.io/entrain/
