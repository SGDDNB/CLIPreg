
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CLIPreg

<!-- badges: start -->

<!-- badges: end -->

The goal of CLIPreg is to discover key RBP regulators in different
datasets. It combines CLIP-seq with RNA- and RIBO-seq to calculate
enrichment of RBP and generate plots for publications.

## Installation

### Check and install required packages

Users may use following codes to check and install all the required
packages.

``` r
list.of.packages <- c("ggplot2","grid","doParallel","foreach","data.table","fastmatch","GGally","ggnet","topGO","ALL","devtools")

## for package "ggplot2", "pheatmap", "grid", "doParallel", "foreach", "data.table", "fastmatch", "GGally"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## for package "topGO", "ALL", "ggnet", "ComplexHeatmap"
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))  BiocManager::install(new.packages)
if("ggnet"%in%new.packages)  devtools::install_github("briatte/ggnet")
install_version("network", version = "1.16.1", repos = "http://cran.us.r-project.org")
if("ComplexHeatmap"%in%new.packages)  devtools::install_github("jokergoo/ComplexHeatmap")
```

### Install CLIPreg

The source code of CLIPreg can be installed from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("SGDDNB/CLIPreg")
```

CLIPreg requires 4 different inputs. A cluster input which is a
dataframe containing geneID and the cluster given by DeltaTE output
(ideally). Curated CLIP-seq data, POSTAR and ENCODE are pre-loaded in
the package. RIBO lfc and TPM.

``` r
## load libraries
library(CLIPreg)
library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(doParallel)
library(ggnet)
library(topGO)
library(GGally)
library(data.table)
library(stringr)
```

### Basic usage

For testing the package with the example data

``` r
#Example=Load_example()
```

For using your own data, you must specify a folder that contains at
least 3 txt files that are named: clusters.txt, ribo\_lfc.txt and
ribo\_tpm.txt. For the format of those 3 files, please refer to Advance
usage step 1.

``` r
#Input_data=Load_input_files(folder = "path/to/folder")
```

Run the analysis with all default parameters.

``` r
#results=run_CLIPreg(Example, is.example=T) # or run_CLIPreg(Input_data, is.example=F)
```

Generate the visual output from the results of the analysis. The
Visualise function will create 2 new files in the folder given: an Rdata
file containing the list object from the run\_CLIPreg function and a pdf
file with the main figures.

``` r
#dir.create("Results_CLIPreg")
#Visualise(results=results,folder="Results_CLIPreg")
```

### Advance usage

## Step 1: Input datasets

Let’s have a look at the cluster file from the example data. It consists
of 2 columns containing the geneID and the cluster for all the DE genes.
This is optional, if you want to provide your own clusters, go to next
step.

``` r
data("clusters")
head(clusters)
#>             geneID      Cluster
#> 1:       LOC728392 forwarded_up
#> 2:      WRB-SH3BGR forwarded_up
#> 3: ENSG00000162408 forwarded_up
#> 4: ENSG00000204859 forwarded_up
#> 5: ENSG00000142583 forwarded_up
#> 6: ENSG00000116649 forwarded_up
```

You can input your own clusters by using the function load\_clusters()
and giving the file location of your file as input.

``` r
#clusters=load_clusters(cluster_file = "cluster_file")
```

Load POSTAR and ENCODE RBP data. Those are 2 public datasets which are
processed in order to have lists of vector. Each vector is named after 1
RBP and contains the geneID of all the targets of that RBP. Combine both
data in a target list.

``` r
data("RBP_ENCODE")
data("RBP_POSTAR")
Targets=combine_targets(RBP_list1=RBP_ENCODE,RBP_list2=RBP_POSTAR,background=clusters$geneID)
```

Load the fold change and identify the RBPs. If you have your own then
provide your own data.

``` r
# load fold change and tpm. optional if you want to use the input data
data("ribo_lfc")
data("tpm_ribo")

head(ribo_lfc)
#>                          geneID IDENTIFIER  FoldChange
#> ENSG00000000003 ENSG00000000003     TSPAN6  0.36690410
#> ENSG00000000419 ENSG00000000419       DPM1  0.10266923
#> ENSG00000000457 ENSG00000000457      SCYL3 -0.42380808
#> ENSG00000000971 ENSG00000000971        CFH -0.21626115
#> ENSG00000001036 ENSG00000001036      FUCA2 -0.04897555
#> ENSG00000001084 ENSG00000001084       GCLC -0.10033803
#head(tpm_ribo)

# If you want to input your own data
#load_ribo_lfc(ribo_lfc_file = "ribo_lfc_file")
#load_ribo_tpm(ribo_tpm_file = "ribo_tpm_file")
```

## Step 2: Data integration and analysis

Run the enrichment analysis using CLIPReg\_V4() function. This takes
several minutes. If you want to have a look at the example results skip
this step.

``` r
# The CLIPreg_V4 function requires a few minutes to run, save the data after running it to be sure not to lose it
# res_Postar=CLIPReg_V4(RBP_data=RBP_POSTAR,clusters=clusters)
# res_Encode=CLIPReg_V4(RBP_data=RBP_ENCODE,clusters=clusters)
# save(res_Encode,file="Res_RBP_Encode.RData")
# save(res_Postar,file="Res_RBP_Postar.RData")
```

If you want to get the results directly you can load it by using the
example data results. The output of CLIPReg\_V4() is a list of
dataframes. One dataframe per cluster containing the RBP and statistical
information calculated during the analysis such as p-value and z-score.

``` r
data("Res_Encode")
data("Res_Postar")
head(res_Encode[[1]])
#>       RBP real_overlap simulated_overlap_mean simulated_overlap_sd          z
#> 1 ZC3H11A          268              282.41101            12.590167 -1.1446242
#> 2    GNL3           14               15.73466             3.306621 -0.5246020
#> 3  HNRNPM          345              396.49065            14.077142 -3.6577488
#> 4   RBM15          461              476.04563            14.839723 -1.0138754
#> 5   DDX24          370              406.79695            14.246241 -2.5829235
#> 6   XRCC6          142              146.71529             9.573704 -0.4925252
#>      pval
#> 1 0.86594
#> 2 0.63920
#> 3 0.99986
#> 4 0.83574
#> 5 0.99474
#> 6 0.66889
```

Then we want to combine POSTAR and ENCODE to work with only one
dataframe and only keep RBPs that are significant in at least one
cluster.

``` r

res=CLIPreg::combine(res1=res_Encode,res2=res_Postar)
```

Extract the RBP LFC from the RIBO\_LFC and keep only detected RBPs in
res

``` r
# Change of RBPs
rbp_lfc=rbp_change(res=res,ribo_lfc=ribo_lfc)

# Cure res data by removing RBPs that are not in the rbp_lfc dataframe
res=cure_res(res=res,rbp_lfc=rbp_lfc)
```

## Step 3: Visualisation

Generate and save heatmap to pdf. The heatmap represents the -logP of
each RBP for each cluster. The blue RBPs are downregulated and the
orange RBPs are upregulated. Only RBPs that are significant in at least
one cluster are shown.

``` r
# Heatmap of RBP scores

HeatmapRBP(res=res,rbp_lfc=rbp_lfc)
```

![](man/figures/README-Generate%20and%20save%20heatmap-1.png)<!-- -->

``` r
# Plot the heatmap with updated colors for the RBPs

# Save the heatmap
# e=HeatmapRBP(res=res,rbp_lfc=rbp_lfc)
# location="Heatmap_fibroblasts.pdf"
# n=length(e$tree_row$order)
# pdf(location,8,3+n*0.15)
# 
# dev.off()
```

A bubble plot can be generated only if the clusters are the output of
DeltaTE program.

``` r
# Bubble plot clusters if clusters are from DeltaTE
BubbleRBPs(res = res,clusters = clusters,rbp_lfc = rbp_lfc)
```

![](man/figures/README-Bubble%20plot-1.png)<!-- -->

From the results, the user can choose a number of RBP to draw the
network for by n. This will pick the n most changing RBPs.

``` r
# Draw network
Draw_network_by_group(rbp_lfc=rbp_lfc,res=res,Targets=Targets,clusters=clusters,n=5,forwarded = F)
```

![](man/figures/README-Network-1.png)<!-- -->

``` r
# plot GO
Plot_GO(rbp_lfc=rbp_lfc,res=res,Targets=Targets,clusters=clusters,n=5,
  tpm_ribo = tpm_ribo,th=200,GO_to_show=3,forwarded = F)
```

![](man/figures/README-GO%20of%20specific%20nodes-1.png)<!-- -->

``` r
Plot_GO_node_name(rbp_lfc=rbp_lfc,res=res,Targets=Targets,clusters=clusters,n=5,
                  tpm_ribo = tpm_ribo,Nodes_to_keep=c(19,15),GO_to_show=3,forwarded = F)
```

![](man/figures/README-GO%20of%20specific%20nodes-2.png)<!-- -->

``` r
Plot_GO_RBP(rbp_of_interest="QKI",tpm_ribo = tpm_ribo,Targets=Targets,clusters=clusters,GO_to_show=3)
```

![](man/figures/README-GO%20of%20specific%20nodes-3.png)<!-- -->
