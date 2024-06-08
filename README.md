# Generation of spermatogonia from pluripotent stem cells in humans and non-human primates

# Abstract

Failures in germline development drive male infertility, but the lack of model systems that recapitulate human spermatogenesis hamper therapeutic development. Here, we develop a system to differentiate human induced pluripotent stem cells (iPSCs) into primordial germ cell-like cells that self-organize with mouse fetal testicular cells into proper seminiferous tubule architecture within xenogeneic reconstituted testes (xrTestes). Subsequent transplant of xrTestes into immunodeficient mice results in efficient generation of spermatogonial stem cells and differentiating spermatogonia as well as preleptotene spermatocytes with striking transcriptomic and phenotypic similarities to their in vivo counterparts. As future clinical application will require testing in non-human primates, we utilize a similar strategy to differentiate rhesus iPSCs through all fetal germ cell stages into spermatogonia-like cells. Together these models will provide an unparalleled opportunity to manipulate fetal and adult primate germ cell lineages and represent a substantial step towards the goal of in vitro gametogenesis.

# bioRxiv link

https://www.biorxiv.org/content/10.1101/2024.05.03.592203v1

# Methods

Sequence data was demultiplexed and fastqs produced using Cell Ranger (v.7.1.0). For human samples, reads were mapped to the GRCh38 human reference and independently mapped to the GRCm38 mouse reference. Similarly, macaque samples were aligned to the Mmul_10 reference, created using the Mmul_10 genome assembly and gtf file using Cell Ranger’s mkref function, as well as to the mouse reference. This produced two separate alignments for each cell, one for the human/macaque reference (as appropriate for human or macaque samples) and another for the mouse reference. Both sets of UMI counts were read into R (v.4.2.2) and secondary analysis was performed using Seurat (v.4.4.0). UMI count tables were first loaded into R by using the Read10X function, and Seurat objects were built from each sample, one per alignment. The sum of UMI counts for each alignment was totaled for each cell for both human/macaque and mouse counts. The total UMIs for human/macaque alignment were then compared with the mouse alignment. Any cell with greater than twofold assignment of UMIs to human/macaque than mouse was designated as human or macaque as appropriate for the sample. Similarly, any cell with twofold higher mouse counts than primate would be designated mouse. Any cell that failed to meet either threshold was labeled “undetermined” and excluded from the analysis. Human (iPSC and in vivo 12-year-old testis) samples and mouse (two independent adult testes) samples were used as controls. No human or mouse control cells were assigned to be the wrong species (Figure XX). In this fashion, all cells with “mouse” or “undetermined” identities were excluded, and only human cells aligned to the human reference or macaque cells aligned to the macaque reference were kept for downstream analysis. For each sample a minimum number of genes = 100 and maximum genes = 3000 were set as well as a maximum mitochondrial content of 20% and cells outside those bounds were excluded. A minimum UMI was set on a per-sample basis to exclude low-quality cells by assigning a minimum cutoff based on an inspection of a ranked plot of UMI per sample. For the sample NCG8, a cutoff of 103.5 was set, for NCG5, NCG19 and NCG_c41 a cutoff of 103.7 was used, while for NCG4, NCG9 and NCG_182_2 a cutoff of 103.8 was applied, and for all other samples a minimum of 104 UMI/cell was set. Samples were integrated using Seurat’s IntegrateData function using 3000 highly-variable genes and the first 30 dimensions. Cell cycle state was assigned using Seurat’s CellCycleScoring function and samples were split by phase and re-integrated to regress out cell cycle effects. Doublets were identified using DoubletFinder with an expected doublet rate of 4%. This identified three clusters which showed a high rate of doublets and with an intermediate identity when looking at key marker genes, so these clusters were removed, and the dataset was re-clustered using 30 principle components. Following these filtering steps, 10,298 cells were brought forward for further analysis, grouped via Seurat into 18 clusters. RNA velocity was determined using velocyto with the parameters deltaT = 10, kCells = 25,  fit.quantile = 0.01. Pseudotime was generated using monocle3. As PGCLCs showed a clear disjunction from the unbroken progression of all other germ cells due to a lack of intermediate cell types captured, they were excluded before running pseudotime, which was generated using monocle’s learn_graph function with ncenter=50, minimal_branch_len=1. Using previously-identified marker genes for PGCs, prospermatogonia and postnatal germ cell types, clusters were assigned to one of ten cell types, ranging from early PGCs to preleptotene spermatocytes. A random selection of MEIOB transcripts were selected by searching fastq files for a 3’ MEIOB sequence and aligned to the human MEIOB and mouse Meiob genes, and the cells containing these transcripts were identified via cell barcode. These showed 100% match to human MEIOB and only 35% identity with mouse Meiob via BLAST for the 3’ region of the gene. When the cells that contain these human-aligned MEIOB molecules were mapped back onto the UMAP, they corresponded to the region of MEIOB expression.

DEGs between groups were determined using Wilcoxon Rank Sum test via Seurat’s FindMarkers function. For pseudobulk analysis, mean normalized counts were used by cell type. Gene ontology enrichment was analyzed via DAVID (v2023q4).  Ingenuity Pathway Analysis (QIAGEN) was used for all pathway analyses. 

To examine the heterogeneity of the spermatogonial compartment, T1.prosp., spermatogonia and Diff.spg.E were subset, reclustered, and pseudotime run as described above. Four clusters were identified that corresponded to T1 prospermatogonia, S0 spermatogonia, S1 spermatogonia and early differentiating spermatogonia. 

To compare xrTestis with normal in vivo spermatogenesis we used an atlas of human germ cells we have derived previously (how do we cite this?). Briefly, testis cells derived from between gestational week (GW) 7 and GW 17 fetuses, along with cells from patients from shortly after birth to 12 years old were analyzed via scRNAseq. In vivo samples were first integrated together in the same manner as described for xrTestes. Somatic cells, identified by key gene expression, were eliminated and the remainder reclustered and cell cycle phase state assigned via Seurat as described previously. In vivo germ cells were then split by cell cycle phase and integrated with xrTestis data, also split by phase. The integrated object was clustered by Monocle 3 and a UMAP generated by Seurat and cell types were assigned to these new clusters anew, replacing prior designations. Similarities and differences between in vivo and xrTestis clusters were therefore based on shared clusters assigned together. Differentially expressed markers were defined using FindMarkers and pseudobulk DEGs as described previously. 

Analysis of macaque xrTestis-derived germ cells was performed in the same manner as human  germ cells. For comparison of macaque and human datasets, data transfer was employed using macaque dataset as a query against the human datasets, using a subset of homologous genes with identical genes. Seurat’s MapQuery and TransferData functions were used to project macaque cells onto human dimensional reduction/UMAP.


# Software used

R version 4.2.2 (2022-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.2.1

dplyr version 1.1.4
Seurat version 5.1.0
SeuratWrappers version 0.3.1
patchwork version 1.2.0
gplots version 3.1.3.1
ggplot2 version 3.5.1
SeuratObject version 5.0.2
monocle3 version 1.3.1
cowplot version 1.1.3
viridis version 0.6.5
tidyverse version 2.0.0
forcats version 1.0.0
scCustomize version 2.1.2
stringr version 1.5.1
ComplexHeatmap version 2.14.0
circlize version 0.4.16
lattice version 0.22-6
scales version 1.3.0
ggridges version 0.5.6
reshape2 version 1.4.4
ggrastr version 1.0.2
tibble version 3.2.1
scico version 1.5.0
RColorBrewer version 1.1-3
DoubletFinder version 2.0.3
