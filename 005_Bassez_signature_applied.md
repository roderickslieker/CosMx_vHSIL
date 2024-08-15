# Bassez data


-   [Bassez](#bassez)
    -   [Apply signatures](#apply-signatures)
-   [Subset data](#subset-data)
-   [Rename](#rename)
-   [Check](#check)
-   [Change levels](#change-levels)
-   [Frequency across conditions](#frequency-across-conditions)
-   [Add meta info](#add-meta-info)
-   [Add PID](#add-pid)
-   [Main classes](#main-classes)
    -   [Plots](#plots)
-   [Subtypes](#subtypes)

``` r
dir <- "/Users/rcslieker/Documents/ONCO/002_Projects/002_CosMx_signatures/"
knitr::opts_knit$set(root.dir = dir)
setwd(dir)
```

## Bassez

Downloaded from: https://lambrechtslab.sites.vib.be/en/data-access

### Apply signatures

``` r
library(ggplot2)
library(magrittr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rstatix)

source("./001_Scripts/Utils.R")
```

## Subset data

``` r
load("../000_Data/002_Public_Datasets/Bassez_Breast/Cohort1.RData")
```

## Rename

``` r
load("./002_Data/Bassez_Immune_Applied.RData")
Bassez_applied$mainClass <- gsub("plasmablast","plasma cell",Bassez_applied$mainClass)
Bassez_applied$subClass <- gsub("plasmablast","plasma cell",Bassez_applied$subClass)
```

## Check

``` r
t.sub <- table(Bassez_applied$cellType, Bassez_applied$mainClass)  %>% as.data.frame()
t.all <- table(Bassez_applied$mainClass) %>% as.data.frame()
t.sub$Total <- t.all[match(t.sub$Var2, t.all$Var1),2]

t.sub$Percentage <- (t.sub$Freq / t.sub$Total)*100
```

## Change levels

``` r
t.sub <- t.sub[!t.sub$Var2 %in% "Unassigned",]
lvl <- c("B-cell","plasma cell","T CD4","T CD8","Treg","macrophage","neutrophil","mDC","pDC","mast", "fibroblast","endothelial","a","b","c","d","e","f")
t.sub$Var2 <- factor(t.sub$Var2, levels=lvl)
t.sub$Var1 <- factor(t.sub$Var1, levels=rev(c("B_cell","T_cell","Myeloid_cell","pDC","Mast_cell",
                                          "Fibroblast","Endothelial_cell","Cancer_cell")))
```

``` r
ggplot(t.sub, aes(x=Var2, y=Var1, fill=Percentage))+
  geom_tile()+
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle=90, vjust=0.3, hjust=1))+
  scale_fill_gradientn(colours = c("white","blue","black"))
```

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-6-1.png)

## Frequency across conditions

``` r
Bassez_applied$Tumour_ID <- paste0(Bassez_applied$patient_id, "_", Bassez_applied$timepoint)

freqs_main_all <- table(Bassez_applied$mainClass, Bassez_applied$Tumour_ID) %>% reshape2::melt()
freqs_main_overall <- table(Bassez_applied$Tumour_ID) %>% reshape2::melt()

freqs_sub_all <- table(Bassez_applied$subClass, Bassez_applied$Tumour_ID) %>% reshape2::melt()
freqs_sub_overall <- table(Bassez_applied$Tumour_ID) %>% reshape2::melt()
```

## Add meta info

``` r
freqs_main_all$Total <- freqs_main_overall[match(freqs_main_all$Var2, freqs_main_overall$Var1),"value"]
freqs_sub_all$Total <- freqs_sub_overall[match(freqs_sub_all$Var2, freqs_sub_overall$Var1),"value"]


freqs_main_all <- data.frame(freqs_main_all, Bassez_applied[match(freqs_main_all$Var2, Bassez_applied$Tumour_ID),c("BC_type","timepoint")])
freqs_sub_all <- data.frame(freqs_sub_all, Bassez_applied[match(freqs_sub_all$Var2, Bassez_applied$Tumour_ID),c("BC_type","timepoint")])
freqs_main_all$Percentage <- (freqs_main_all$value / freqs_main_all$Total)*100
freqs_sub_all$Percentage <- (freqs_sub_all$value / freqs_sub_all$Total)*100
freqs_main_all$timepoint <- factor(freqs_main_all$timepoint, levels=c("Pre","On"))
freqs_sub_all$timepoint <- factor(freqs_sub_all$timepoint, levels=c("Pre","On"))
```

## Add PID

``` r
freqs_main_all$patient_id <- Bassez_applied[match(freqs_main_all$Var2, Bassez_applied$Tumour_ID),"patient_id"]
freqs_sub_all$patient_id <- Bassez_applied[match(freqs_sub_all$Var2, Bassez_applied$Tumour_ID),"patient_id"]
```

## Main classes

### Plots

``` r
lvl <- levels(freqs_main_all$Var1)
lvl <- lvl[!lvl %in% "Unassigned"]

#i = "T CD4"

pvals.main <- lapply(lvl, function(i){
  cat(i,"\n")
  fma <- freqs_main_all[freqs_main_all$Var1 %in% i,]
  

  if(sum(fma$Percentage) == 0) return()
  fma <- fma[order(fma$patient_id),]
 

  fit <- rstatix::wilcox_test(data = fma, formula = Percentage~timepoint, p.adjust.method = "none", paired = TRUE)
  fit
  fit$CellType <- i
  fit
}) %>% do.call(what=rbind)
```

    a 
    b 
    B-cell 
    c 
    d 
    e 
    endothelial 
    f 
    fibroblast 
    macrophage 
    mast 
    mDC 
    neutrophil 
    pDC 
    plasma cell 
    T CD4 
    T CD8 
    Treg 

``` r
mx <- by(freqs_main_all$Percentage, freqs_main_all$Var1, function(i){quantile(i, probs=0.95)}) %>% as.list() %>% do.call(what=c)
```

``` r
cts <- as.character(unique(freqs_main_all$Var1))
cts <- cts[!cts %in% "Unassigned"]


lapply(cts, function(CT){
  cat(CT)
  fma <- freqs_main_all[freqs_main_all$Var1 %in% CT,]

  # Create a box plot
  bxp <- ggboxplot(outlier.shape=NA,
    fma, x = "timepoint", y = "Percentage",
    fill = "timepoint", palette = c("#96C66A", "#7172B4")
    )
  
  pvals.sub <- pvals.main[pvals.main$CellType %in% CT,]
  
  
  bp <- boxplot(fma$Percentage, plot=F) # box plot is a simple way to identify outliers
  if(length(bp$out) > 0){outlier <- which(fma$Percentage %in% bp$out) # get the position of those outliers in your data
    
   maxv <- by(fma[-outlier,]$Percentage, fma[-outlier,]$timepoint, function(i){quantile(i, probs=0.95)}) %>% as.list() %>% do.call(what=c) %>% max(na.rm=T)
  }else{
    maxv <- by(fma$Percentage, fma$timepoint, function(i){quantile(i, probs=0.95)}) %>% as.list() %>% do.call(what=c) %>% max(na.rm=T)
  }
  
  
  seqx <- (0.1*c(1,3,5,7,9,11))
  
  pvals.sub$y.position <- (maxv)*seqx[1:nrow(pvals.sub)]
  #seqx <- seqx[-(1:nrow(pvals.sub))]
  #pvals.sub$y.position <- (maxv)*seqx[1:nrow(pvals.sub)]

  
  pvals.sub <- add_significance(pvals.sub)

  if(nrow(pvals.sub) >= 1){
      maxx <- max(pvals.sub$y.position)
  }else{
    maxx <- maxv
  }
  pvals.sub$x <- 1.5
  pvals.sub$xmin <- 1
  pvals.sub$xmax <- 2
  pvals.sub <- pvals.sub[pvals.sub$p <= 0.05,]
  bxp + coord_cartesian(ylim = c(0, maxx)) +stat_pvalue_manual(
    pvals.sub,  label = "p.signif", tip.length = 0)+
    theme(legend.position = "bottom")+
    ggtitle(CT)
})
```

    abB-cellcdeendothelialffibroblastmacrophagemastmDCneutrophilpDCplasma cellT CD4T CD8Treg

    [[1]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-1.png)


    [[2]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-2.png)


    [[3]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-3.png)


    [[4]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-4.png)


    [[5]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-5.png)


    [[6]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-6.png)


    [[7]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-7.png)


    [[8]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-8.png)


    [[9]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-9.png)


    [[10]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-10.png)


    [[11]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-11.png)


    [[12]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-12.png)


    [[13]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-13.png)


    [[14]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-14.png)


    [[15]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-15.png)


    [[16]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-16.png)


    [[17]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-17.png)


    [[18]]

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-11-18.png)

## Subtypes

``` r
lvl.sub <- levels(freqs_sub_all$Var1)
lvl.sub <- lvl.sub[grep("_",lvl.sub)]

pvals.sub <- lapply(lvl.sub, function(i){
  cat(i,"\n")
  fma <- freqs_sub_all[freqs_sub_all$Var1 %in% i,]
  
  # Pre / on within E /NE
  fma.sub <- fma#[fma$expansion %in% ep & fma$expansion %in% c("NE","E"),]
  if(sum(fma.sub$Percentage) == 0) return()
  fit <- rstatix::wilcox_test(data = fma.sub, formula = Percentage~timepoint, p.adjust.method = "none", paired = TRUE)
  fit
  

  
  
  
  if(is.null(fit)) return()
  fit$CellType <- i
  fit
}) %>% do.call(what=rbind)
```

    a_0 
    a_1 
    a_10 
    a_2 
    a_3 
    a_4 
    a_5 
    a_6 
    a_7 
    a_8 
    a_9 
    b_0 
    b_1 
    b_10 
    b_11 
    b_12 
    b_13 
    b_2 
    b_3 
    b_4 
    b_5 
    b_6 
    b_7 
    b_8 
    b_9 
    B-cell_0 
    B-cell_1 
    B-cell_2 
    B-cell_3 
    B-cell_4 
    B-cell_5 
    c_0 
    c_1 
    c_10 
    c_11 
    c_12 
    c_2 
    c_3 
    c_4 
    c_5 
    c_6 
    c_7 
    c_8 
    c_9 
    d_0 
    d_1 
    d_10 
    d_11 
    d_12 
    d_2 
    d_3 
    d_4 
    d_6 
    d_9 
    e_0 
    e_1 
    e_2 
    e_3 
    e_6 
    e_7 
    e_8 
    e_9 
    endothelial_0 
    endothelial_1 
    endothelial_10 
    endothelial_2 
    endothelial_3 
    endothelial_4 
    endothelial_6 
    endothelial_7 
    endothelial_8 
    endothelial_9 
    f_6 
    f_7 
    fibroblast_0 
    fibroblast_1 
    fibroblast_2 
    fibroblast_3 
    fibroblast_4 
    fibroblast_5 
    fibroblast_6 
    fibroblast_7 
    fibroblast_8 
    fibroblast_9 
    macrophage_0 
    macrophage_1 
    macrophage_2 
    macrophage_3 
    macrophage_4 
    macrophage_5 
    macrophage_6 
    macrophage_7 
    macrophage_8 
    macrophage_9 
    mast_0 
    mast_1 
    mast_2 
    mast_4 
    mast_5 
    mast_7 
    mDC_0 
    mDC_1 
    mDC_2 
    mDC_3 
    mDC_4 
    mDC_5 
    mDC_6 
    mDC_7 
    mDC_8 
    neutrophil_0 
    neutrophil_1 
    neutrophil_2 
    pDC_0 
    pDC_1 
    pDC_2 
    pDC_4 
    pDC_5 
    plasma cell_0 
    plasma cell_1 
    plasma cell_2 
    plasma cell_3 
    plasma cell_4 
    plasma cell_5 
    plasma cell_6 
    plasma cell_7 
    plasma cell_8 
    T CD4_0 
    T CD4_1 
    T CD4_10 
    T CD4_2 
    T CD4_3 
    T CD4_4 
    T CD4_5 
    T CD4_6 
    T CD4_7 
    T CD4_8 
    T CD4_9 
    T CD8_0 
    T CD8_1 
    T CD8_2 
    T CD8_3 
    T CD8_4 
    T CD8_5 
    T CD8_6 
    T CD8_7 
    Treg_0 
    Treg_1 
    Treg_2 
    Treg_3 
    Treg_4 
    Treg_5 
    Treg_6 
    Treg_7 

``` r
rio::export(pvals.sub, file="./003_Figures/Bassez_noSplit_significant_subsets.xlsx")
```

``` r
library(patchwork)

getPlot <-function(ct.in){
  cts <- freqs_sub_all$Var1[grep(paste0(ct.in,"_"), freqs_sub_all$Var1)] %>% as.character() %>% unique()
  
  cts.o <- paste0(ct.in,"_", 0:20)
  cts.o <- cts.o[cts.o %in% cts]
  
  lapply(cts.o, function(CT){
    cat(CT)
    fma <- freqs_sub_all[freqs_sub_all$Var1 %in% CT,]
  
    # Create a box plot
    bxp <- ggboxplot(outlier.shape=NA,
      fma, x = "timepoint", y = "Percentage",
      fill = "timepoint", palette = c("#96C66A", "#7172B4")
      )
    
    
    
    
    
    pvals.sub.sub <- pvals.sub[pvals.sub$CellType %in% CT,]
  
    
    bp <- boxplot(fma$Percentage, plot=F) # box plot is a simple way to identify outliers
    if(length(bp$out) > 0){outlier <- which(fma$Percentage %in% bp$out) # get the position of those outliers in your data
    
    maxv <- by(fma[-outlier,]$Percentage, fma[-outlier,]$timepoint, function(i){quantile(i, probs=0.99)}) %>% as.list() %>% do.call(what=c) %>% max(na.rm=T)
    }else{
    maxv <- by(fma$Percentage, fma$timepoint, function(i){quantile(i, probs=0.99)}) %>% as.list() %>% do.call(what=c) %>% max(na.rm=T)
     
    }
    

    #maxv <- by(fma$Percentage, fma$timepoint, function(i){quantile(i, probs=0.95)}) %>% as.list() %>% do.call(what=c) %>% max(na.rm=T)
    
    
    seqx <- 1+(0.1*c(1,3,5,7,9,11))
    
    pvals.sub.sub$y.position <- (maxv)*seqx[1:nrow(pvals.sub.sub)]
    seqx <- seqx[-(1:nrow(pvals.sub.sub))]
    pvals.sub.sub$y.position <- (maxv)*seqx[1:nrow(pvals.sub.sub)]
  
    
    pvals.sub.sub <- add_significance(pvals.sub.sub)
  
    if(nrow(pvals.sub.sub) >= 1){
        maxx <- max(pvals.sub.sub$y.position)
    }else{
      maxx <- maxv
    }
    pvals.sub.sub$x <- 1.5
    pvals.sub.sub$xmin <- 1
    pvals.sub.sub$xmax <- 2
    pvals.sub.sub <- pvals.sub.sub[pvals.sub.sub$p <= 0.05,]
    bxp <- bxp + coord_cartesian(ylim = c(0, maxx)) +stat_pvalue_manual(
      pvals.sub.sub,  label = "p.signif", tip.length = 0)+
      theme(legend.position = "bottom")+
      ggtitle(CT)
    bxp
}) %>% wrap_plots(ncol=4)

}
```

``` r
getPlot("T CD4")
```

    T CD4_0T CD4_1T CD4_2T CD4_3T CD4_4T CD4_5T CD4_6T CD4_7T CD4_8T CD4_9T CD4_10

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-14-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_TCD4.pdf", width=15, height=14)
getPlot("T CD4")
```

    T CD4_0T CD4_1T CD4_2T CD4_3T CD4_4T CD4_5T CD4_6T CD4_7T CD4_8T CD4_9T CD4_10

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("T CD8")
```

    T CD8_0T CD8_1T CD8_2T CD8_3T CD8_4T CD8_5T CD8_6T CD8_7

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-15-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_TCD8.pdf", width=15, height=10)
getPlot("T CD8")
```

    T CD8_0T CD8_1T CD8_2T CD8_3T CD8_4T CD8_5T CD8_6T CD8_7

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("Treg")
```

    Treg_0Treg_1Treg_2Treg_3Treg_4Treg_5Treg_6Treg_7

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-16-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_Treg.pdf", width=15, height=9)
getPlot("Treg")
```

    Treg_0Treg_1Treg_2Treg_3Treg_4Treg_5Treg_6Treg_7

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("B-cell")
```

    B-cell_0B-cell_1B-cell_2B-cell_3B-cell_4B-cell_5

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-17-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_Bcell.pdf", width=15, height=9)
getPlot("B-cell")
```

    B-cell_0B-cell_1B-cell_2B-cell_3B-cell_4B-cell_5

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("endothelial")
```

    endothelial_0endothelial_1endothelial_2endothelial_3endothelial_4endothelial_6endothelial_7endothelial_8endothelial_9endothelial_10

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-18-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_endothelial.pdf", width=15, height=13)
getPlot("endothelial")
```

    endothelial_0endothelial_1endothelial_2endothelial_3endothelial_4endothelial_6endothelial_7endothelial_8endothelial_9endothelial_10

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("fibroblast")
```

    fibroblast_0fibroblast_1fibroblast_2fibroblast_3fibroblast_4fibroblast_5fibroblast_6fibroblast_7fibroblast_8fibroblast_9

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-19-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_fibroblast.pdf", width=15, height=13)
getPlot("fibroblast")
```

    fibroblast_0fibroblast_1fibroblast_2fibroblast_3fibroblast_4fibroblast_5fibroblast_6fibroblast_7fibroblast_8fibroblast_9

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("macrophage")
```

    macrophage_0macrophage_1macrophage_2macrophage_3macrophage_4macrophage_5macrophage_6macrophage_7macrophage_8macrophage_9

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-20-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_macrophage.pdf", width=15, height=13)
getPlot("macrophage")
```

    macrophage_0macrophage_1macrophage_2macrophage_3macrophage_4macrophage_5macrophage_6macrophage_7macrophage_8macrophage_9

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("mast")
```

    mast_0mast_1mast_2mast_4mast_5mast_7

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-21-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_mast.pdf", width=15, height=9)
getPlot("mast")
```

    mast_0mast_1mast_2mast_4mast_5mast_7

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("mDC")
```

    mDC_0mDC_1mDC_2mDC_3mDC_4mDC_5mDC_6mDC_7mDC_8

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-22-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_mDC.pdf", width=15, height=14)
getPlot("mDC")
```

    mDC_0mDC_1mDC_2mDC_3mDC_4mDC_5mDC_6mDC_7mDC_8

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("pDC")
```

    pDC_0pDC_1pDC_2pDC_4pDC_5

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-23-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_pDC.pdf", width=15, height=9)
getPlot("pDC")
```

    pDC_0pDC_1pDC_2pDC_4pDC_5

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("plasma cell")
```

    plasma cell_0plasma cell_1plasma cell_2plasma cell_3plasma cell_4plasma cell_5plasma cell_6plasma cell_7plasma cell_8

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-24-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_plasmacell.pdf", width=15, height=13)
getPlot("plasma cell")
```

    plasma cell_0plasma cell_1plasma cell_2plasma cell_3plasma cell_4plasma cell_5plasma cell_6plasma cell_7plasma cell_8

``` r
dev.off()
```

    quartz_off_screen 
                    2 

``` r
getPlot("neutrophil")
```

    neutrophil_0neutrophil_1neutrophil_2

![](005_Bassez_signature_applied.markdown_strict_files/figure-markdown_strict/unnamed-chunk-25-1.png)

``` r
pdf("./004_PaperFigures/Bassez_noSplit_neutrophil.pdf", width=15, height=5)
getPlot("neutrophil")
```

    neutrophil_0neutrophil_1neutrophil_2

``` r
dev.off()
```

    quartz_off_screen 
                    2 
