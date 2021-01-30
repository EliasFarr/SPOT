### Workflow
![](Figures/Figure1.PNG)


### Files

File | Use
--- | ---
ui.R | user interface (UI)
server.R | server code
helper_module.R | shiny modules and functions for UI and server
sc_P_berghei_counts.csv | *P. berghei* scRNAseq data from Howick *et al.* link XXX
sc_P_berghei_averaged.csv | *P. berghei* scRNAseq data, averaged for every lifecycle stage (annotation shortenedLifestage4 from Howick *et al.* 2019 link XXX)
sc_P_berghei_dotplot.csv | *P. berghei* scRNAseq data in dot plot file format 
sc_P_berghei_UMAP.csv | UMAP coordinates *P. berghei* scRNAseq data
bulk_H_sapiens_averaged.csv | *H. sapiens* bulkRNAseq data from Cardoso *et al.*, averaged for every developmental stage 


### Uploadfile formats

Files for upload should have entities as first row and a gene identifier as first column (see below and Example file). Counts should be in TPM normalization.
![](Figures/Figure2.PNG)

### Requirements

All code was excecuted with R 4.0.3. Package requirements are listed below:

package | version
--- | ---
DESeq2	| 1.30.0
DT	| 0.16
EnvStats	| 2.4.0
Limma| 3.46.0
MAST	| 1.16.0
matrixStats	| 0.57.0
Plotly	| 4.9.2.1
Plyr	| 1.8.6
Reshape2	| 1.4.4
Rlist| 	0.4.6.1
Scales	| 1.1.1
Seurat	| 3.2.2
Shiny	| 1.5.0
Shinybusy	| 0.2.2
Shinycssloaders	| 1.0.0
shinyEffects	| 0.1.0
Shinywidgets	| 0.5.4
Sortable	| 0.4.4
Stringr	| 1.4.0
Tidyr	| 1.3.0



