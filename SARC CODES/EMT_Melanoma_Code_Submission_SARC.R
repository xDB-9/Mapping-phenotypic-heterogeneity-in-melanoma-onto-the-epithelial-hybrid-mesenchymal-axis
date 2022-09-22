library(matrixStats)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(ggprism)
library(patchwork)
library(magrittr)
library(matrixStats)
library(gplots)
library(RColorBrewer)
library(fgsea)
library(GSEABase)
library(enrichplot)
library(GSVA)
library(plotly)
library(scatterplot3d)
library(wesanderson)
library(paletteer)
library(readxl)
#Loading all the necessary Libraries

makeExpressionSet <- function(dat, state=colnames(dat)){
    dat <- data.matrix(dat)
    pdata <- as.data.frame(state)
    rownames(pdata) <- colnames(dat)
    
    metadata <- data.frame(labelDescription=c("state"), row.names=c("state"))
    phenoData <- new("AnnotatedDataFrame", data=pdata, varMetadata=metadata)
    dataExp <- ExpressionSet(assayData=dat,  phenoData=phenoData)
    dataExp
} #Function to convert a data frame to expression class

TCGA_SARC <- read_csv("~/Desktop/EMT-Melanoma/GDC TCGA Sarcoma (SARC)/TCGA-SARC.csv") #Reading the TCGA Dataset
EM_gene_signature_tumor_KS <- read_excel("~/Desktop/EMT-Melanoma/DataSets/EM_gene_signature_tumor_KS.xlsx",
    col_names = FALSE) #Reading the Epithelial and Mesenchymal Genelist
mapping_file<- read.delim(file.choose(), stringsAsFactor = FALSE) #Reading gencode.v22.annotation.gene.probeMap file

colnames(mapping_file)<-c("Ensembl_ID","gene","chrom","chromStart","chromEnd","strand") #Renaming the columns for mapping the ENSEMBL ID with Gene Name
TCGA_SARC_mapped<-merge(x = TCGA_SARC, y = mapping_file, by="Ensembl_ID") #Mapping the ENSEMBL ID with Gene Name

TCGA_SARC_mapped<-TCGA_SARC_mapped[,-c(1,268,269,270,271)] #Data Cleaning and processing

TCGA_SARC_discrete<-aggregate(TCGA_SARC_mapped,by=list(TCGA_SARC_mapped$gene), FUN=mean,na.action=na.omit) #Aggregating the similar gene name expression levels and storing them once

TCGA_SARC_discrete<-TCGA_SARC_discrete[,-267] #Data Cleaning and processing

row.names(TCGA_SARC_discrete)<-TCGA_SARC_discrete$Group.1 #Column to RowName
TCGA_SARC_discrete<-TCGA_SARC_discrete[,-1] #Column to RowName

EM_gene_signature_tumor_KS<-na.omit(EM_gene_signature_tumor_KS) #Data Cleaning and processing
EM_gene_signature_tumor_KS_Epi <- EM_gene_signature_tumor_KS[EM_gene_signature_tumor_KS$...2 == "Epi", ] #EM_gene_signature_tumor_KS_Mes_list<-list(c("GAS1", "CXCL12", "ZEB1", "GLYR1","FHL1", "FERMT2", "C1S", "FYN", "WIPF1", "CYP1B1", "SERPING1","SERPINF1", "VCAM1", "MAP1B", "TCF4", "SRPX", "EMP3", "DPT","CALD1", "PTGIS", "VIM", "CD163", "C1R", "FBN1", "FN1", "FXYD6","IGF1", "NAP1L3", "MRC1", "QKI", "MS4A4A", "DCN", "LOX", "RECK","ANK2", "LY96", "ZFPM2", "CSRP2", "EFEMP1", "RARRES2", "PTPRC","PLEKHO1", "RGS2", "F13A1", "JAM2", "CHRDL1", "TUBA1A", "AP1S2","MYLK", "DDR2", "DSE", "SACS", "GLIPR1", "CXCL13", "FLRT2", "PTX3","AKT3", "COL6A2", "DPYSL3", "CDH11", "PDZRN3", "ZEB2", "CCL2","MAFB", "SFRP1", "C14orf139", "MFAP4", "MAF", "UCHL1", "TUBB6","SRGN", "HEG1", "KCNJ8", "AKAP12", "EVI2A", "COL14A1", "AXL","ECM2", "FSTL1", "PLN", "MYL9", "OLFML3", "STON1", "SLIT2", "BICC1","SOBP", "CLIC4", "ENPP2", "SAMSN1", "TPM2", "ASPN", "COL6A1","IGFBP5", "MOXD1", "AKAP2", "SLC2A3", "OLFML2B", "ANGPTL2", "PCOLCE","COLEC12", "CTSK", "TAGLN", "CDH2", "IL10RA", "C1orf54", "CEP170","TNS1", "CLEC2B", "JAM3", "38961", "GREM1", "VCAN", "ZCCHC24","CRYAB", "SFRP4", "RUNX1T1", "FGL2", "MS4A6A", "PTRF", "GIMAP4","TWIST1", "GFPT2", "LHFP", "CXCR4", "SPOCK1", "SPARC", "VSIG4","GPM6B", "TRPC1", "SNAI2", "GUCY1B3", "PLXNC1", "SYT11", "FLI1","MYH10", "CSF2RB", "TNC", "PMP22", "COL5A2", "MMP2", "GNG11","CAV1", "CDK14", "SDC2", "PTGDS", "NR3C1", "SYNM", "FAP", "NUAK1","WWTR1", "FBLN1", "MPDZ", "SYNE1", "EFEMP2", "GIMAP6", "KIAA1462","CCL8", "COL15A1", "CHN1", "CRISPLD2", "PDGFC", "GEM", "ISLR","GZMK", "SPARCL1", "BNC2", "BGN", "MEOX2", "ITM2A", "IFFO1"))
EM_gene_signature_tumor_KS_Mes <- EM_gene_signature_tumor_KS[EM_gene_signature_tumor_KS$...2 == "Mes", ] #EM_gene_signature_tumor_KS_Epi_list<-list(c("KRT19", "AGR2", "RAB25", "CDH1","ERBB3", "FXYD3", "SLC44A4", "S100P", "SCNN1A", "GALNT3", "PRSS8","ELF3", "CEACAM6", "TMPRSS4", "CLDN7", "TACSTD2", "CLDN3", "EPCAM","SPINT1", "TSPAN1", "PLS1", "TMEM30B", "PRR15L", "KRT8", "ST14", "RBM47", "S100A14", "C1orf106", "NQO1", "TOX3", "PTK6", "TFF1","CLDN4", "GPRC5A", "TJP3", "KRT18", "MAP7", "CKMT1A", "ESRP1","MUC1", "SPINT2", "ESRP2", "CDS1", "PPAP2C", "CEACAM7", "TTC39A","OVOL2", "EHF", "AP1M2", "CEACAM5", "LAD1", "ARHGAP8", "TFF3","JUP", "CD24", "TMC5", "MLPH", "ELMO3", "ERBB2", "LLGL2", "DDR1","FA2H", "CBLC", "TMPRSS2", "LSR", "PERP", "POF1B", "MYO5C", "RAB11FIP1","MAPK13", "KRT7", "CEACAM1", "CXADR", "ATP2C2", "RNF128", "MPZL2","EPS8L1", "GALNT7", "CORO2A", "BCAS1", "TPD52", "ARHGAP32", "FUT2","OR7E14P", "GALE", "GRHL2", "BIK", "RAPGEFL1", "STYK1", "F11R","PKP3", "CYB561", "SH3YL1", "GDF15", "PSCA", "EZR", "TJP2", "FGFR3","FUT3", "BSPRY", "TOM1L1", "IRF6", "EPB41L4B", "OCLN", "LRRC1","C19orf21", "ABHD11", "EPS8L2", "MYO6", "TSPAN8", "MST1R", "SLC16A5","GPR56", "AZGP1", "TOB1", "SLC35A3", "TRPM4", "PHLDA2", "VAMP8","SLC22A18", "AKR1B10", "VAV3", "SPAG1", "ABCC3", "SYNGR2", "STAP2","C4orf19", "PPL", "PLLP", "DSG2", "HDHD3", "CD2AP", "MANSC1","DHCR24", "EPN3", "TUFT1", "GMDS", "EXPH5", "DSP", "SDC4", "IL20RA","FAM174B", "PTPRF", "SORD"))

TCGA_SARC_expr <- makeExpressionSet(TCGA_SARC_discrete, state=colnames(TCGA_SARC_discrete)) #Converting to dataexpression class for ssGSEA score calculation

GSVAtumor_SARC_Epi<-gsva(TCGA_SARC_expr, EM_gene_signature_tumor_KS_Epi_list, method=c("ssgsea")) #Calculating ssGSEA score
GSVAtumor_SARC_Mes<-gsva(TCGA_SARC_expr, EM_gene_signature_tumor_KS_Mes_list, method=c("ssgsea")) #Calculating ssGSEA score

scatter_epi_mes_data<-data.frame(rbind(GSVAtumor_SARC_Mes@assayData[["exprs"]],GSVAtumor_SARC_Epi@assayData[["exprs"]])) #Merging the two ssGSEA scores calculated above
row.names(scatter_epi_mes_data)<-c("Mesenchymal","Epithelial") #Giving appropriate RowNames
tscatter_epi_mes_data<-data.frame(t(scatter_epi_mes_data)) #Taking transpose

for ( row in 1:nrow(tscatter_epi_mes_data)){
    row.names(tscatter_epi_mes_data)[row] <-  sub(".01A", ".01", row.names(tscatter_epi_mes_data)[row])
} #Removing A at the end in the rowname
for ( row in 1:nrow(tscatter_epi_mes_data)){
    row.names(tscatter_epi_mes_data)[row] <-  sub(".01B", ".01", row.names(tscatter_epi_mes_data)[row])
} #Removing B at the end in the rowname

sample_subtype <- read_excel("~/Desktop/EMT-Melanoma/DataSets/sample_subtype.xlsx")
for ( row in 1:nrow(sample_subtype)){
    row.names(sample_subtype)[row] <-  sub("-", ".", sample_subtype$`TCGA barcode`[row])
} #Column to rowname with appropriate changes
for ( row in 1:nrow(sample_subtype)){
    row.names(sample_subtype)[row] <-  sub("-", ".", row.names(sample_subtype)[row])
} #making necesarry changes in the rownames before mapping with the subtype
for ( row in 1:nrow(sample_subtype)){
    row.names(sample_subtype)[row] <-  sub("-", ".", row.names(sample_subtype)[row])
} #making necesarry changes in the rownames before mapping with the subtype

sample_subtype$`TCGA barcode`<-row.names(sample_subtype) #rownames to column for mapping
tscatter_epi_mes_data$samples<-row.names(tscatter_epi_mes_data) #rownames to column for mapping
colnames(sample_subtype)<-c("samples","subtypes") #Giving appropriate ColumnNames for mapping
tscatter_epi_mes_data<-data.frame(merge(tscatter_epi_mes_data,sample_subtype,by= "samples")) #Merging the two datasets based on column

colnames(tscatter_epi_mes_data)<-c("Samples","Mesenchymal","Epithelial","Subtypes") #Giving appropriate ColumnNames
ggplot_epi_mes_data<-ggplot(tscatter_epi_mes_data, aes(x=Epithelial, y=Mesenchymal, color=Subtypes)) + geom_point() + geom_rug() + theme_minimal() + labs(title="Subtypes across \nEpithelial and Mesenchymal Axes")
ggplot_epi_mes_data + scale_color_brewer(palette="Dark2")+  theme(text = element_text(size = 20))  + ggpubr::border() #Plotting with appropriate theme


