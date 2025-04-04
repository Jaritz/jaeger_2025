# Jaeger et al., submitted
Snakemake/R/Pyhton code used for total RNA-Seq analysis of Jaeger et al. "Chemically-induced protein degradation reveals 1 inflammation-dependent requirement for TREG lineage-defining transcription factor FOXP3". Submitted, 2025.
Christina Jäger, Polina Dimitrova1, Qiong Sun, Jesse Tennebroek, Elisa Marchiori, Markus Jaritz, Rene Rauschmeier, Guillem Estivill, Anna Obenauf, Meinrad Busslinger, Joris van der Veeken,  [van der Veeken lab](https://www.imp.ac.at//groups/joris-van-der-veeken)

* **Snakefile** contains the *Snakemake* file to run the individual analyses.
* **rTregaTreg** contains snakemake config file and sample tables for analysing differential gene expression in steady-state or activated GFP+ CD4 T cells isolated from pooled secondary lymphoid tissues of animals treated with 5-Ph-IAA or PBS, all normalized together.
* **rTregaTreg_byCellType** contains snakemake config file and sample tables for analysing differential gene expression in steady-state or activated GFP+ CD4 T cells isolated from pooled secondary lymphoid tissues of animals treated with 5-Ph-IAA or PBS, with celltypes normalized individually.
* **Treg_Il2aIl2** contains snakemake config file and sample tables for analysing differential gene expression for differentially expressed genes in GFP+ TREG cells isolated from mice injected with IL-2/anti-IL2 complexes with or without 5-Ph-IAA.
* **code** contains additinonal R and python code for the individual processing steps.

Note: these scripts are optimized for runnig on the Vienna Biocenter Computing cluster (CBE) using SLURM.

Original implementation: Kimon Froussios

Updates/changes (adding all-noexon and UMI handling): Markus Jaritz, [van der Veeken lab](https://www.imp.ac.at//groups/joris-van-der-veeken)
