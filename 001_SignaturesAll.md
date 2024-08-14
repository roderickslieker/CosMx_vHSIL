# Cell transfer


-   [Packages](#packages)
-   [Load data and filter set](#load-data-and-filter-set)
    -   [Conde](#conde)
    -   [Bassez](#bassez)
    -   [Franken](#franken)
    -   [Kong](#kong)

``` r
dir <- "/Users/rcslieker/Documents/ONCO/002_Projects/002_CosMx_signatures//"
knitr::opts_knit$set(root.dir = dir)
setwd(dir)
```

## Packages

``` r
library(SingleR)
library(Seurat)
library(patchwork)
library(magrittr)
library(ggplot2)
source("./001_Scripts/Utils.R")
```

# Load data and filter set

``` r
load("../001_CosMxData/002_Data/CosMx_withPCA.Rdata")
```

## Conde

``` r
load("../000_Data/002_Public_Datasets/Conde_Immune/Conde_immune_CosMx_subset.RData")
Conde_immune_applied <- getSeurat_all(CosMx, Conde_immune)
save(Conde_immune_applied, file="./002_Data/Conde_Immune_Applied.RData")
```

## Bassez

``` r
load("../000_Data/002_Public_Datasets/Bassez_Breast/Cohort1.RData")
Cohort1 <- subset(Cohort1, features = rownames(CosMx))
Bassez_applied <- getSeurat_all(CosMx, Cohort1)
save(Bassez_applied, file="./002_Data/Bassez_Immune_Applied.RData")
```

``` r
load("../000_Data/002_Public_Datasets/Bassez_Breast/Cohort2.RData")
Bassez_applied_2 <- getSeurat_all(CosMx, Cohort2)
save(Bassez_applied_2, file="./002_Data/Bassez_Immune2_Applied.RData")
```

## Franken

``` r
load("../000_Data/002_Public_Datasets/Franken_Immunity/Data/Seurat_Franken.RData")
Franken_applied <- getSeurat_all(CosMx, se)
save(Franken_applied, file="./002_Data/Franken_Applied.RData")
```

## Kong

``` r
load("../000_Data/002_Public_Datasets/Kong_CD/SCP1884/Kong_Data.RData")
Kong_applied <- getSeurat_all(CosMx, cd)
save(Kong_applied, file="./002_Data/Kong_Applied.RData")
```
