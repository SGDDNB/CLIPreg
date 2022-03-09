
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CLIPreg

<!-- badges: start -->

<!-- badges: end -->

The goal of CLIPreg is to discover key RBP regulators in different
datasets. It combines CLIP-seq with RNA- and RIBO-seq to calculate
enrichment of RBP and generate plots for publications. Another feature
that can be analyzed by CLIPreg is enrichment of miRNA targets from
TargetScan database.

## Installation

### Check and install required packages

Users may use following codes to check and install all the required
packages.

``` r
list.of.packages <- c("ggplot2","grid","doParallel","foreach","data.table","fastmatch","GGally","ggnet","topGO","ALL","devtools","org.Hs.eg.db")

## for package "ggplot2", "pheatmap", "grid", "doParallel", "foreach", "data.table", "fastmatch", "GGally"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## for package "topGO", "ALL", "ggnet", "ComplexHeatmap","org.Hs.eg.db"
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))  BiocManager::install(new.packages)
if("ggnet"%in%new.packages)  devtools::install_github("briatte/ggnet")
devtools::install_version("network", version = "1.16.1", repos = "http://cran.us.r-project.org")
if("ComplexHeatmap"%in%new.packages)  devtools::install_github("jokergoo/ComplexHeatmap")
```

### Install CLIPreg

The source code of CLIPreg can be installed from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("SGDDNB/CLIPreg")
```

CLIPreg requires 4 different inputs. A gene\_groups input which is a
dataframe containing geneID and the gene\_groups given by DeltaTE output
(ideally). The DeltaTE method can be found in the paper here:
<https://doi.org/10.1002/cpmb.108>. It categorizes transcriptionally
and/or translationally regulated genes into 4 categories: forwarded,
intensified, exclusive and buffered, each category with an up or down
direction. A deltaTE function is included in this package, it can be
found in advance usage step 0. The 3 other required files are :
summarized CLIP-seq data, with POSTAR and ENCODE which are pre-loaded in
the package, RIBO log fold change (lfc) and TPM.

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
library(DESeq2)
```

### Basic usage

For testing the package with the example data. This will load the
summary data of POSTAR and ENCODE as well as the TPM, lfc and
gene\_groups from the fibroblast dataset used in the paper.

``` r
Example=Load_example()
```

For using your own data, you must specify a folder that contains at
least 3 txt files that are named: gene\_groups.txt, ribo\_lfc.txt and
ribo\_tpm.txt. For the format of those 3 files, please refer to Advance
usage step 1.

``` r
Input_data=Load_input_files(folder = "path/to/folder")
```

Run the analysis with all default parameters.

``` r
results=run_CLIPreg(Example, is.example=T) # or run_CLIPreg(Input_data, is.example=F)
```

Generate the visual output from the results of the analysis. The
Visualise function will create 2 new files in the folder given: an Rdata
file containing the list object from the run\_CLIPreg function and a pdf
file with the main figures.

``` r
dir.create("Results_CLIPreg")
Visualise(results=results,folder="Results_CLIPreg")
```

### Advance usage

#### Step 0 : Run deltaTE

This has to be done on your dataset from raw counts of RIBO and RNA as
well as a coldata dataframe that can run through DESeq2. Put batch = 1
if you have a batch column in your design. The counts must have a first
column containing the geneID or gene names, they must be unique. The
other columns are the counts of RNA and RIBO. The column names for the
counts must be the same names as filled in the “SampleID” column of
coldata. The coldata dataframe must have column names “SampleID”,
“SeqType”, “Condition” and “Batch”. (Batch is an optional column).
Condition must contain 1 and 2. For example, in a untreated vs treatment
drug experiment, set the untreated sample condition to 1 and treatment
to 2.

``` r
gene_groups=deltaTE_gene_groups(counts=counts,coldata=coldata,batch=0)
```

#### Step 1: Input datasets

Let’s have a look at the gene\_groups file from the example data. It
consists of 2 columns containing the “geneID” and the “Gene\_group” for
all the DE genes. This is as an example on fibroblasts stimulated by
TGFB, if you start with your own gene groups, go to next step.

``` r
data("gene_groups")
head(gene_groups)
#>             geneID   Gene_group
#> 1:       LOC728392 forwarded_up
#> 2:      WRB-SH3BGR forwarded_up
#> 3: ENSG00000162408 forwarded_up
#> 4: ENSG00000204859 forwarded_up
#> 5: ENSG00000142583 forwarded_up
#> 6: ENSG00000116649 forwarded_up
```

You can input your own gene groups by using the function
load\_gene\_groups() and giving the file location of your file as input.
Each gene group must contain “up” or “down” in the name or you won’t be
able to plot the network and gene ontology.

``` r
gene_groups=load_gene_groups(gene_groups_file = "path/to/gene_groups_file.txt")
```

Load POSTAR and ENCODE RBP data. Those are 2 public datasets which are
processed in order to have lists of vector. Each vector is named after 1
RBP and contains the geneID of all the targets of that RBP. Combine both
data in a target list.

``` r
data("RBP_ENCODE")
data("RBP_POSTAR")
Targets=combine_targets(RBP_list1=RBP_ENCODE,RBP_list2=RBP_POSTAR,background=gene_groups$geneID)
```

Load the fold change and identify the RBPs. If you have your own then
provide your own data.

``` r
# load fold change and tpm. optional if you want to use the input data
data("ribo_lfc")
data("tpm_ribo")

head(ribo_lfc)
#>                          geneID IDENTIFIER FoldChange
#> ENSG00000000003 ENSG00000000003     TSPAN6  0.3669041
#> ENSG00000000419 ENSG00000000419       DPM1  0.1026692
#> ENSG00000000457 ENSG00000000457      SCYL3 -0.4238081
#> ENSG00000000971 ENSG00000000971        CFH -0.2162611
#> ENSG00000001084 ENSG00000001084       GCLC -0.1003380
#> ENSG00000001461 ENSG00000001461     NIPAL3 -0.3866424
#head(tpm_ribo)
```

To load your own data use

``` r
# If you want to input your own data
load_ribo_lfc(ribo_lfc_file = "ribo_lfc_file")
load_ribo_tpm(ribo_tpm_file = "ribo_tpm_file")
```

#### Step 2: Data integration and analysis

Run the enrichment analysis using CLIPreg() function. This takes several
minutes depending on how many gene groups and how many genes per gene
group there are. If you want to have a look at the example results skip
this step.

``` r
# The CLIPreg function requires a few minutes to run, save the data after running it to be sure not to lose it
res_Postar=CLIPreg(RBP_data=RBP_POSTAR,gene_groups=gene_groups)
res_Encode=CLIPreg(RBP_data=RBP_ENCODE,gene_groups=gene_groups)
#res_miR=CLIPreg(RBP_data=miR_data,gene_groups=gene_groups)
save(res_Encode,file="Res_RBP_Encode.RData")
save(res_Postar,file="Res_RBP_Postar.RData")
#save(res_miR,file="Res_miR.RData")
```

If you want to get the results directly you can load it by using the
example data results. The output of CLIPreg() is a list of dataframes.
One dataframe per gene group containing the RBP and statistical
information calculated during the analysis such as p-value and z-score.

``` r
data("res_Encode")
data("res_Postar")
head(res_Encode[[1]])
#>       RBP real_overlap simulated_overlap_mean simulated_overlap_sd          z
#> 1 ZC3H11A          268              282.41101            12.590167 -1.1446242
#> 2    GNL3           14               15.73466             3.306621 -0.5246020
#> 3  HNRNPM          345              396.49065            14.077142 -3.6577488
#> 4   RBM15          461              476.04563            14.839723 -1.0138754
#> 5   DDX24          370              406.79695            14.246241 -2.5829235
#> 6   XRCC6          142              146.71529             9.573704 -0.4925252
#>      pval padj
#> 1 0.86594    1
#> 2 0.63920    1
#> 3 0.99986    1
#> 4 0.83574    1
#> 5 0.99474    1
#> 6 0.66889    1
```

Then we want to combine POSTAR and ENCODE to work with only one
dataframe and only keep RBPs that are significant in at least one gene
group.

``` r

res=CLIPreg::combine(res1=res_Encode,res2=res_Postar,FDR=0.1)
head(res)
#> $forwarded_up
#>             RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 1        ALKBH1           40               40.97943             5.275220
#> 2        ALKBH5          599              646.88781            15.573995
#> 3         ATXN2         1231             1249.49058             9.418429
#> 4      C17orf85          759              809.32365            15.414678
#> 5       CAPRIN1          958             1014.84253            13.853814
#> 6         CELF2          597              673.85913            15.640422
#> 7          CNBP          278              304.99502            12.946082
#> 8         CPSF1          940             1000.85109            14.017409
#> 9         CPSF2          500              568.45342            15.363176
#> 10        CPSF3         1011             1062.38911            13.289958
#> 11        CPSF4          928              992.26129            14.131421
#> 12        CPSF6         1014             1058.46603            13.259708
#> 13        CPSF7          995             1049.55567            13.415458
#> 14        CSTF2         1145             1193.97796            10.879823
#> 15       CSTF2T         1173             1205.54419            10.551160
#> 16        DDX3X         1219             1236.99381             9.771163
#> 17        DGCR8          239              259.88132            12.153228
#> 18       DICER1          395              430.35728            14.432834
#> 19       DIS3L2          110              122.63193             8.843790
#> 20        EIF3A          532              537.61172            15.235891
#> 21        EIF3B          535              598.81263            15.479669
#> 22        EIF3D          621              653.39952            15.604828
#> 23        EIF3G          804              837.26724            15.295357
#> 24       EIF4A3          464              529.01438            15.214827
#> 25       ELAVL1         1257             1289.30913             8.108882
#> 26        EWSR1          812              882.33524            15.011065
#> 27         EZH2            4                5.04836             1.883119
#> 28          FBL          283              319.88571            13.196883
#> 29       FIP1L1         1217             1252.17583             9.343733
#> 30         FMR1          824              911.52147            14.847267
#> 31          FUS         1117             1163.98890            11.480535
#> 32         FXR1          130              199.24487            10.973454
#> 33         FXR2          426              526.85663            15.143310
#> 34      HNRNPA1          992             1046.27125            13.434858
#> 35    HNRNPA2B1           84              112.52593             8.538031
#> 36       HNRNPC         1197             1227.23871            10.037684
#> 37       HNRNPD          725              805.37998            15.422630
#> 38       HNRNPF          183              245.58702            11.929098
#> 39       HNRNPH          190              203.16129            11.085509
#> 40       HNRNPM          206              264.34497            12.219267
#> 41       HNRNPU          423              516.42804            15.120107
#> 42      IGF2BP1          804              892.14414            14.939946
#> 43      IGF2BP2          765              853.23669            15.125983
#> 44      IGF2BP3          879              972.65978            14.245710
#> 45       LIN28A          873              936.40805            14.628477
#> 46       LIN28B         1217             1251.30444             9.362961
#> 47        MBNL2          277              380.43798            13.928886
#> 48      METTL14           80               95.64639             7.914270
#> 49       METTL3          161              177.04546            10.406970
#> 50        MOV10         1005             1062.63564            13.268308
#> 51        NOP56          229              285.07747            12.625234
#> 52        NOP58          268              298.77354            12.835065
#> 53       NUDT21         1191             1235.82224             9.804341
#> 54        PRKRA           60               65.00330             6.599562
#> 55        PTBP1         1193             1219.49851            10.205809
#> 56   PTBP1PTBP2          598              664.66506            15.572072
#> 57         PUM2          334              423.47940            14.381784
#> 58          QKI          317              407.66707            14.191437
#> 59        RBM10          178              203.13212            11.114106
#> 60       RBM15B         1021             1077.21737            13.046817
#> 61        RBM15          811              864.79221            15.151666
#> 62        RBM47          360              473.67510            14.823691
#> 63        RBPMS          245              297.53573            12.841437
#> 64         RTCB          834              889.54657            14.938049
#> 65        SRRM4          674              727.35025            15.592666
#> 66        SRSF1          910              977.39691            14.262009
#> 67          SSB          379              373.00692            13.812845
#> 68        STAU1          500              563.68047            15.332555
#> 69        TAF15          598              693.22360            15.603948
#> 70       TARBP2          375              403.31574            14.202676
#> 71       TARDBP         1266             1288.11378             8.148649
#> 72         TIA1          541              618.04386            15.496954
#> 73        TIAL1          804              871.95967            15.110438
#> 74       TNRC6A           43               49.89956             5.806728
#> 75       TNRC6B           33               40.06955             5.242269
#> 76       TNRC6C           31               34.43893             4.838692
#> 77        U2AF2         1125             1174.62577            11.291837
#> 78         UPF1          978             1039.51535            13.537217
#> 79        WDR33         1003             1057.26291            13.334117
#> 80         WTAP           89               97.13432             7.983991
#> 81       YTHDC1          994             1056.42888            13.347192
#> 82       YTHDC2          269              289.26126            12.697292
#> 83       YTHDF1         1068             1119.38490            12.362137
#> 84       YTHDF2         1129             1166.34334            11.453455
#> 85       YTHDF3          618              670.00857            15.600567
#> 86       ZC3H7B          926              977.99133            14.233938
#> 87      ZC3H11A          268              282.41101            12.590167
#> 510       DDX24          370              406.79695            14.246241
#> 610       XRCC6          142              146.71529             9.573704
#> 710       SF3B4          644              678.87034            15.680244
#> 91        NCBP2          350              377.18147            13.872721
#> 101     FAM120A          458              506.71068            15.046409
#> 111       LARP7           38               36.52782             4.984062
#> 121         FTO           91               95.66366             7.900697
#> 141       RBM22          365              395.88323            14.130708
#> 151       UCHL5          469              513.40313            15.130823
#> 161       DDX42           59               52.56453             5.955049
#> 171        XPO5          250              247.42674            11.926124
#> 181      TROVE2           76               90.86991             7.721426
#> 191        SUB1          264              281.88869            12.579825
#> 201       TBRG4          360              410.40015            14.259464
#> 241      METAP2          121              147.00859             9.605555
#> 251        DKC1          175              167.83299            10.178550
#> 271     KHDRBS1           39               56.13002             6.157265
#> 281        SLTM          176              173.16894            10.311872
#> 291        SND1          275              280.43357            12.550879
#> 311      SERBP1           49               64.15501             6.542334
#> 321        NONO          163              178.79456            10.470182
#> 341       NSUN2           11               13.97398             3.126660
#> 361        PPIG          606              616.81013            15.542345
#> 381      AKAP8L          123              130.33604             9.118056
#> 391      HNRNPL          286              362.30943            13.688775
#> 401       PCBP2          505              559.58582            15.337322
#> 411       SF3A3          413              416.62310            14.310073
#> 431      RBFOX2          359              395.00414            14.095427
#> 451      HNRNPK          290              290.09581            12.691712
#> 471      GEMIN5          367              389.36022            14.075901
#> 481       NOL12           78              100.69574             8.095322
#> 491     SUPV3L1           68               74.23677             7.043703
#> 501       SRSF7          105              100.97112             8.122575
#> 521      GTF2F1          312              346.60956            13.525562
#> 531      EFTUD2          653              703.26779            15.674007
#> 541       PRPF8          896              957.17074            14.438729
#> 551       DDX55          365              374.51490            13.872856
#> 581       BUD13          732              763.82023            15.607458
#> 591       YWHAG            8               14.83608             3.206340
#> 601       FUBP3          167              220.41421            11.384820
#> 611       MATR3          238              284.22742            12.613256
#> 621       TRA2A          136              138.70236             9.419076
#> 631      DROSHA          293              311.83765            13.066934
#> 651       SRSF9           79               80.79503             7.304603
#> 671       KHSRP          376              470.45951            14.835055
#> 681      EXOSC5          113              152.38700             9.756906
#> 691       LARP4          319              339.80580            13.416817
#> 701       CDC40          103              104.56010             8.225064
#> 711       LSM11          222              223.02849            11.464125
#> 721     FASTKD2          183              207.59288            11.145928
#> 741      SMNDC1          169              184.65339            10.631862
#> 751       DHX30          108              115.79241             8.640163
#> 771        RPS3          449              449.09232            14.573602
#> 781        XRN2          360              392.24590            14.070507
#> 791        SBDS            3                2.67244             1.367451
#> 801        RBM5          101              119.11030             8.759311
#> 811        NPM1            7               10.39802             2.688489
#> 821      ZNF622          422              463.27868            14.775255
#> 841       SUGP2          193              236.10828            11.751321
#> 861       NOLC1          163              174.29833            10.349737
#> 871        YBX3          481              517.67104            15.100574
#> 88        SF3B1          675              748.72617            15.522102
#> 89        DDX59          132              152.58389             9.797790
#> 90        GRSF1           74               86.72707             7.548596
#> 93         NKRF          203              209.97086            11.176346
#> 96       ZRANB2          145              161.87870             9.995831
#> 97        MTPAP           74               85.49531             7.517549
#> 98         PUS1            7                7.71937             2.324205
#> 99        SAFB2           53               63.55648             6.501998
#> 103       BCCIP           29               27.30080             4.347423
#> 104    HNRNPUL1           65               59.69178             6.343978
#> 109        AARS          168              173.73775            10.283820
#> 1111       DDX6          370              383.14422            13.933490
#> 112      EIF4G2          163              182.35542            10.579130
#> 115       EIF3H          192              214.13742            11.293253
#> 116       FKBP4           64               63.84357             6.537739
#> 117       GRWD1          567              606.75289            15.487481
#>                 z    pval padj
#> 1    -0.185666193 0.53155    1
#> 2    -3.074857230 0.99897    1
#> 3    -1.963234068 0.97092    1
#> 4    -3.264657887 0.99938    1
#> 5    -4.103023795 0.99997    1
#> 6    -4.914133915 1.00000    1
#> 7    -2.085188448 0.98022    1
#> 8    -4.341108198 0.99999    1
#> 9    -4.455681618 0.99999    1
#> 10   -3.866762342 0.99993    1
#> 11   -4.547404737 1.00000    1
#> 12   -3.353469889 0.99949    1
#> 13   -4.066627612 0.99992    1
#> 14   -4.501723792 1.00000    1
#> 15   -3.084418162 0.99867    1
#> 16   -1.841521833 0.96123    1
#> 17   -1.718170671 0.95268    1
#> 18   -2.449780845 0.99226    1
#> 19   -1.428339027 0.91575    1
#> 20   -0.368322400 0.63002    1
#> 21   -4.122351011 0.99997    1
#> 22   -2.076249727 0.97964    1
#> 23   -2.174989459 0.98421    1
#> 24   -4.273093539 0.99999    1
#> 25   -3.984412435 0.99996    1
#> 26   -4.685559671 1.00000    1
#> 27   -0.556714716 0.59958    1
#> 28   -2.795031911 0.99732    1
#> 29   -3.764644163 0.99991    1
#> 30   -5.894786481 1.00000    1
#> 31   -4.092919155 0.99995    1
#> 32   -6.310216451 1.00000    1
#> 33   -6.660144162 1.00000    1
#> 34   -4.039585061 0.99998    1
#> 35   -3.341043077 0.99960    1
#> 36   -3.012518700 0.99841    1
#> 37   -5.211820578 1.00000    1
#> 38   -5.246584233 1.00000    1
#> 39   -1.187251780 0.87425    1
#> 40   -4.774833990 1.00000    1
#> 41   -6.179059282 1.00000    1
#> 42   -5.899896705 1.00000    1
#> 43   -5.833451438 1.00000    1
#> 44   -6.574595661 1.00000    1
#> 45   -4.334562589 1.00000    1
#> 46   -3.663845350 0.99975    1
#> 47   -7.426149017 1.00000    1
#> 48   -1.976984506 0.97321    1
#> 49   -1.541799448 0.93325    1
#> 50   -4.343857496 1.00000    1
#> 51   -4.441697633 1.00000    1
#> 52   -2.397614589 0.99133    1
#> 53   -4.571672890 0.99999    1
#> 54   -0.758126017 0.74993    1
#> 55   -2.596414429 0.99439    1
#> 56   -4.281065542 0.99999    1
#> 57   -6.221717801 1.00000    1
#> 58   -6.388857631 1.00000    1
#> 59   -2.261281227 0.98731    1
#> 60   -4.308895562 0.99997    1
#> 61   -3.550250618 0.99980    1
#> 62   -7.668474696 1.00000    1
#> 63   -4.091109796 0.99998    1
#> 64   -3.718462040 0.99993    1
#> 65   -3.421496269 0.99974    1
#> 66   -4.725625194 1.00000    1
#> 67    0.433877316 0.31858    1
#> 68   -4.153284967 0.99999    1
#> 69   -6.102532454 1.00000    1
#> 70   -1.993690486 0.97594    1
#> 71   -2.713797146 0.99569    1
#> 72   -4.971548659 1.00000    1
#> 73   -4.497531341 0.99999    1
#> 74   -1.188201016 0.86413    1
#> 75   -1.348566731 0.89571    1
#> 76   -0.710714731 0.72634    1
#> 77   -4.394835946 1.00000    1
#> 78   -4.544165260 0.99999    1
#> 79   -4.069478950 0.99996    1
#> 80   -1.018828777 0.83077    1
#> 81   -4.677304552 0.99999    1
#> 82   -1.595714999 0.94125    1
#> 83   -4.156635851 0.99998    1
#> 84   -3.260443271 0.99934    1
#> 85   -3.333761526 0.99945    1
#> 86   -3.652631389 0.99985    1
#> 87   -1.144624211 0.86594    1
#> 510  -2.582923532 0.99474    1
#> 610  -0.492525154 0.66889    1
#> 710  -2.223839141 0.98603    1
#> 91   -1.959346716 0.97290    1
#> 101  -3.237362525 0.99948    1
#> 111   0.295377543 0.34209    1
#> 121  -0.590284598 0.69912    1
#> 141  -2.185540156 0.98467    1
#> 151  -2.934614324 0.99805    1
#> 161   1.080674538 0.12337    1
#> 171   0.215766664 0.39937    1
#> 181  -1.925798368 0.96981    1
#> 191  -1.422014192 0.91634    1
#> 201  -3.534505139 0.99990    1
#> 241  -2.707661465 0.99631    1
#> 251   0.704128802 0.22577    1
#> 271  -2.782082569 0.99703    1
#> 281   0.274543751 0.37087    1
#> 291  -0.432923465 0.65111    1
#> 311  -2.316453234 0.98882    1
#> 321  -1.508527744 0.92884    1
#> 341  -0.951168438 0.78306    1
#> 361  -0.695527624 0.74617    1
#> 381  -0.804561834 0.77266    1
#> 391  -5.574598957 1.00000    1
#> 401  -3.559018965 0.99976    1
#> 411  -0.253185288 0.58614    1
#> 431  -2.554313490 0.99420    1
#> 451  -0.007549021 0.48795    1
#> 471  -1.588546212 0.94057    1
#> 481  -2.803562291 0.99750    1
#> 491  -0.885439146 0.79136    1
#> 501   0.496010190 0.28625    1
#> 521  -2.558826031 0.99430    1
#> 531  -3.207079722 0.99917    1
#> 541  -4.236573757 1.00000    1
#> 551  -0.685864560 0.74036    1
#> 581  -2.038783623 0.97740    1
#> 591  -2.132051201 0.98032    1
#> 601  -4.691704215 1.00000    1
#> 611  -3.664987010 0.99981    1
#> 621  -0.286902869 0.59043    1
#> 631  -1.441627411 0.92003    1
#> 651  -0.245739556 0.56679    1
#> 671  -6.367317751 1.00000    1
#> 681  -4.036832812 0.99997    1
#> 691  -1.550725447 0.93616    1
#> 701  -0.189676326 0.54892    1
#> 711  -0.089713784 0.51789    1
#> 721  -2.206445218 0.98508    1
#> 741  -1.472309416 0.92284    1
#> 751  -0.901882253 0.80165    1
#> 771  -0.006334741 0.48830    1
#> 781  -2.291736860 0.98856    1
#> 791   0.239540518 0.26317    1
#> 801  -2.067548572 0.97855    1
#> 811  -1.263914479 0.86202    1
#> 821  -2.793771071 0.99740    1
#> 841  -3.668377220 0.99986    1
#> 861  -1.091653804 0.85168    1
#> 871  -2.428453433 0.99142    1
#> 88   -4.749754195 1.00000    1
#> 89   -2.100870777 0.98056    1
#> 90   -1.686018058 0.94847    1
#> 93   -0.623715504 0.71688    1
#> 96   -1.688573898 0.95110    1
#> 97   -1.529130037 0.92930    1
#> 98   -0.309512229 0.52663    1
#> 99   -1.623574885 0.94105    1
#> 103   0.390852219 0.30187    1
#> 104   0.836733608 0.17938    1
#> 109  -0.557939563 0.69543    1
#> 1111 -0.943354468 0.81862    1
#> 112  -1.829585196 0.96375    1
#> 115  -1.960234073 0.97332    1
#> 116   0.023927233 0.45649    1
#> 117  -2.566775721 0.99423    1
#> 
#> $forwarded_down
#>             RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 1        ALKBH1           46               47.15942             5.480677
#> 2        ALKBH5          717              744.09239            16.150812
#> 3         ATXN2         1382             1437.37276             9.772669
#> 4      C17orf85          851              931.01479            16.037915
#> 5       CAPRIN1         1148             1167.45517            14.326591
#> 6         CELF2          781              775.22340            16.187287
#> 7          CNBP          401              350.82544            13.448328
#> 8         CPSF1         1108             1151.28655            14.506101
#> 9         CPSF2          688              653.85825            15.945379
#> 10        CPSF3         1186             1222.02223            13.771951
#> 11        CPSF4         1096             1141.50593            14.620235
#> 12        CPSF6         1157             1217.67478            13.830341
#> 13        CPSF7         1161             1207.40675            13.915379
#> 14        CSTF2          258              198.87608            10.688486
#> 15       CSTF2T         1349             1386.76599            10.978578
#> 16        DDX3X         1356             1423.00811            10.128125
#> 17        DGCR8          249              298.94613            12.624506
#> 18       DICER1          491              495.08169            15.011707
#> 19       DIS3L2          151              141.09509             9.190967
#> 20        EIF3A          597              618.44889            15.765281
#> 21        EIF3B          683              688.80125            16.118545
#> 22        EIF3D          708              751.66512            16.191518
#> 23        EIF3G          918              963.14593            15.914377
#> 24       EIF4A3          570              608.52703            15.719874
#> 25       ELAVL1         1487             1483.12504             8.417869
#> 26        EWSR1         1045             1014.99696            15.592938
#> 27         EZH2            4                5.81143             1.952124
#> 28          FBL          364              367.89055            13.729757
#> 29       FIP1L1         1399             1440.42779             9.703168
#> 30         FMR1         1000             1048.55904            15.368085
#> 31          FUS         1328             1338.91918            11.958578
#> 32         FXR1          246              229.28316            11.366818
#> 33         FXR2          586              606.13158            15.705673
#> 34      HNRNPA1         1191             1203.65424            13.959731
#> 35    HNRNPA2B1          146              129.53738             8.803541
#> 36       HNRNPC         1368             1411.69547            10.403662
#> 37       HNRNPD         1015              926.54498            15.979138
#> 38       HNRNPF          342              282.53048            12.409195
#> 39       HNRNPH          236              233.64716            11.464918
#> 40       HNRNPM          364              304.03039            12.725146
#> 41       HNRNPU          683              594.17494            15.610278
#> 42      IGF2BP1         1037             1026.33820            15.492784
#> 43      IGF2BP2         1015              981.57682            15.754263
#> 44      IGF2BP3         1170             1119.01047            14.855440
#> 45       LIN28A         1052             1077.25660            15.225096
#> 46       LIN28B         1404             1439.40218             9.712718
#> 47        MBNL2          520              437.60129            14.427014
#> 48      METTL14          108              110.03479             8.224063
#> 49       METTL3          190              203.64101            10.843020
#> 50        MOV10         1262             1222.42815            13.706518
#> 51        NOP56          365              328.01803            13.075331
#> 52        NOP58          346              343.71150            13.391052
#> 53       NUDT21         1378             1421.64767            10.172357
#> 54        PRKRA           90               74.83666             6.870624
#> 55        PTBP1          643              608.93825            15.696141
#> 56   PTBP1PTBP2          775              764.60506            16.239725
#> 57         PUM2          550              487.19631            14.924347
#> 58          QKI          562              469.02117            14.731150
#> 59        RBM10          298              233.64311            11.452745
#> 60       RBM15B         1139             1239.20272            13.536737
#> 61        RBM15          911              994.90322            15.719052
#> 62        RBM47          657              545.00449            15.389697
#> 63        RBPMS          341              342.26143            13.353425
#> 64         RTCB          964             1023.25828            15.566335
#> 65        SRRM4          733              836.72271            16.226119
#> 66        SRSF1         1057             1124.37026            14.779233
#> 67          SSB          410              429.16696            14.363496
#> 68        STAU1          587              648.54340            15.902687
#> 69        TAF15          856              797.38566            16.220666
#> 70       TARBP2          485              463.99036            14.752991
#> 71       TARDBP          372              351.95840            13.427179
#> 72         TIA1          731              711.00590            16.122499
#> 73        TIAL1          997             1003.10553            15.748774
#> 74       TNRC6A           62               57.38572             6.040078
#> 75       TNRC6B           48               46.11711             5.439003
#> 76       TNRC6C           44               39.63938             5.057848
#> 77        U2AF2         1309             1351.23172            11.707767
#> 78         UPF1         1164             1195.78753            14.108685
#> 79        WDR33         1243             1216.26961            13.815289
#> 80         WTAP          113              111.71576             8.292334
#> 81       YTHDC1         1119             1215.25222            13.835875
#> 82       YTHDC2          297              332.79791            13.171554
#> 83       YTHDF1         1204             1287.67363            12.852510
#> 84       YTHDF2         1269             1341.72823            11.899702
#> 85       YTHDF3          718              770.77654            16.201723
#> 86       ZC3H7B         1173             1125.08914            14.881164
#> 87      ZC3H11A          289              324.92681            13.015945
#> 510       DDX24          379              468.08897            14.783998
#> 610       XRCC6          153              168.79117             9.984329
#> 710       SF3B4          668              781.13502            16.171532
#> 91        NCBP2          419              433.87516            14.431640
#> 101     FAM120A          534              582.90766            15.545144
#> 111       LARP7           34               42.03857             5.199967
#> 121         FTO           89              109.99457             8.191675
#> 141       RBM22          424              455.43362            14.588606
#> 151       UCHL5          548              590.80377            15.698601
#> 161       DDX42           50               60.47211             6.201481
#> 171        XPO5          254              284.62690            12.369899
#> 181      TROVE2          102              104.53533             7.991312
#> 191        SUB1          298              324.33933            13.062773
#> 201       TBRG4          484              472.18198            14.767771
#> 241      METAP2          131              169.13287             9.989987
#> 251        DKC1          186              193.06539            10.556372
#> 271     KHDRBS1           84               64.53993             6.381554
#> 281        SLTM          164              199.16709            10.718851
#> 291        SND1          264              322.53619            13.028970
#> 311      SERBP1           63               73.79813             6.783179
#> 321        NONO          166              205.72349            10.861140
#> 341       NSUN2           19               16.04510             3.241316
#> 361        PPIG          607              709.67078            16.105227
#> 381      AKAP8L           97              150.01127             9.457502
#> 391      HNRNPL          504              416.84374            14.212624
#> 401       PCBP2          572              643.72831            15.897439
#> 411       SF3A3          371              479.34210            14.795761
#> 431      RBFOX2          395              454.39578            14.589627
#> 451      HNRNPK          277              333.79024            13.205312
#> 471      GEMIN5          388              447.97987            14.528061
#> 481       NOL12           73              115.82013             8.402366
#> 491     SUPV3L1           66               85.45194             7.290766
#> 501       SRSF7           93              116.15368             8.437486
#> 521      GTF2F1          351              398.69558            14.058368
#> 531      EFTUD2          728              809.07935            16.167254
#> 541       PRPF8         1038             1101.22035            14.967237
#> 551       DDX55          329              430.80035            14.380343
#> 581       BUD13          790              878.84391            16.146277
#> 591       YWHAG           16               17.09763             3.344803
#> 601       FUBP3          300              253.53414            11.818857
#> 611       MATR3          380              327.03129            13.124608
#> 621       TRA2A          152              159.56321             9.736560
#> 631      DROSHA          290              358.68479            13.529630
#> 651       SRSF9           83               92.94607             7.570683
#> 671       KHSRP          618              541.29419            15.290466
#> 681      EXOSC5          219              175.29841            10.082944
#> 691       LARP4          311              390.83063            13.923736
#> 701       CDC40          105              120.27244             8.545913
#> 711       LSM11          195              256.56345            11.939002
#> 721     FASTKD2          233              238.81415            11.610674
#> 741      SMNDC1          164              212.53583            11.016688
#> 751       DHX30          115              133.22228             8.944085
#> 771        RPS3          431              516.65169            15.097433
#> 781        XRN2          392              451.32775            14.594255
#> 791        SBDS            3                3.08508             1.419423
#> 801        RBM5          114              136.96972             9.066242
#> 811        NPM1           14               11.95159             2.799308
#> 821      ZNF622          507              533.08197            15.298377
#> 841       SUGP2          337              271.66877            12.189306
#> 861       NOLC1          174              200.56711            10.712009
#> 871        YBX3          559              595.58986            15.642586
#> 88        SF3B1          811              861.38168            16.091019
#> 89        DDX59          150              175.61886            10.172311
#> 90        GRSF1           72               99.78378             7.835884
#> 93         NKRF          219              241.58675            11.644410
#> 96       ZRANB2          178              186.21398            10.426727
#> 97        MTPAP           85               98.38524             7.807655
#> 98         PUS1            5                8.87344             2.416696
#> 99        SAFB2           89               73.09560             6.780617
#> 103       BCCIP           28               31.42298             4.509170
#> 104    HNRNPUL1           71               68.65802             6.564876
#> 109        AARS          155              199.91384            10.730626
#> 1111       DDX6          390              440.75127            14.466473
#> 112      EIF4G2          187              209.80694            10.913050
#> 115       EIF3H          225              246.31384            11.748030
#> 116       FKBP4           48               73.41559             6.798841
#> 117       GRWD1          664              698.08929            16.061792
#>                z    pval         padj
#> 1    -0.21154684 0.54360 1.0000000000
#> 2    -1.67746298 0.95035 1.0000000000
#> 3    -5.66608348 1.00000 1.0000000000
#> 4    -4.98910190 1.00000 1.0000000000
#> 5    -1.35797625 0.90698 1.0000000000
#> 6     0.35686030 0.34869 0.9371043750
#> 7     3.73091424 0.00007 0.0006020000
#> 8    -2.98402378 0.99836 1.0000000000
#> 9     2.14116897 0.01430 0.0723411765
#> 10   -2.61562288 0.99469 1.0000000000
#> 11   -3.11253060 0.99894 1.0000000000
#> 12   -4.38707774 1.00000 1.0000000000
#> 13   -3.33492527 0.99946 1.0000000000
#> 14    5.53155253 0.00001 0.0001950000
#> 15   -3.43997101 0.99972 1.0000000000
#> 16   -6.61604314 1.00000 1.0000000000
#> 17   -3.95628379 0.99997 1.0000000000
#> 18   -0.27190046 0.59476 1.0000000000
#> 19    1.07767878 0.12936 0.4449984000
#> 20   -1.36051426 0.90889 1.0000000000
#> 21   -0.35991151 0.63003 1.0000000000
#> 22   -2.69678975 0.99621 1.0000000000
#> 23   -2.83680153 0.99739 1.0000000000
#> 24   -2.45084851 0.99210 1.0000000000
#> 25    0.46032554 0.30246 0.8670520000
#> 26    1.92414281 0.02446 0.1051780000
#> 27   -0.92792757 0.74445 1.0000000000
#> 28   -0.28336626 0.59649 1.0000000000
#> 29   -4.26951182 0.99996 1.0000000000
#> 30   -3.15973266 0.99907 1.0000000000
#> 31   -0.91308345 0.80694 1.0000000000
#> 32    1.47067015 0.06519 0.2548336364
#> 33   -1.28180311 0.89437 1.0000000000
#> 34   -0.90648168 0.80982 1.0000000000
#> 35    1.86999977 0.02805 0.1148714286
#> 36   -4.20000879 0.99997 1.0000000000
#> 37    5.53565657 0.00001 0.0001075000
#> 38    4.79237548 0.00001 0.0001075000
#> 39    0.20522084 0.40213 1.0000000000
#> 40    4.71268535 0.00001 0.0001075000
#> 41    5.69016523 0.00001 0.0001075000
#> 42    0.68817847 0.23577 0.7241507143
#> 43    2.12153245 0.01556 0.0743422222
#> 44    3.43238092 0.00023 0.0016483333
#> 45   -1.65887957 0.94809 1.0000000000
#> 46   -3.64493022 0.99982 1.0000000000
#> 47    5.71141806 0.00001 0.0001075000
#> 48   -0.24741906 0.57110 1.0000000000
#> 49   -1.25804522 0.88792 1.0000000000
#> 50    2.88708262 0.00178 0.0109342857
#> 51    2.82837740 0.00216 0.0123840000
#> 52    0.17089770 0.41792 1.0000000000
#> 53   -4.29081183 0.99997 1.0000000000
#> 54    2.20698154 0.01165 0.0626187500
#> 55    2.17007162 0.01432 0.1523127273
#> 56    0.64009336 0.25100 0.7443448276
#> 57    4.20813648 0.00003 0.0002866667
#> 58    6.31171577 0.00001 0.0001075000
#> 59    5.61934185 0.00001 0.0001075000
#> 60   -7.40228041 1.00000 1.0000000000
#> 61   -5.33767676 1.00000 1.0000000000
#> 62    7.27730435 0.00001 0.0001075000
#> 63   -0.09446491 0.52123 1.0000000000
#> 64   -3.80682279 0.99996 1.0000000000
#> 65   -6.39232998 1.00000 1.0000000000
#> 66   -4.55844085 1.00000 1.0000000000
#> 67   -1.33442168 0.90230 1.0000000000
#> 68   -3.87000015 0.99996 1.0000000000
#> 69    3.61355945 0.00012 0.0009381818
#> 70    1.42409360 0.07320 0.2737043478
#> 71    1.49261436 0.06364 0.4653675000
#> 72    1.24013651 0.10174 0.3645683333
#> 73   -0.38768288 0.64017 1.0000000000
#> 74    0.76394375 0.19850 0.6322592593
#> 75    0.34618295 0.32865 0.9117387097
#> 76    0.86214923 0.16752 0.5541046154
#> 77   -3.60715427 0.99977 1.0000000000
#> 78   -2.25304694 0.98636 1.0000000000
#> 79    1.93484118 0.02373 0.1051780000
#> 80    0.15487076 0.41243 1.0000000000
#> 81   -6.95671357 1.00000 1.0000000000
#> 82   -2.71781986 0.99623 1.0000000000
#> 83   -6.51029505 1.00000 1.0000000000
#> 84   -6.11176903 1.00000 1.0000000000
#> 85   -3.25746464 0.99946 1.0000000000
#> 86    3.21956410 0.00061 0.0040353846
#> 87   -2.76021523 0.99716 1.0000000000
#> 510  -6.02604041 1.00000 1.0000000000
#> 610  -1.58159559 0.93792 1.0000000000
#> 710  -6.99593709 1.00000 1.0000000000
#> 91   -1.03073246 0.84044 1.0000000000
#> 101  -3.14616958 0.99896 1.0000000000
#> 111  -1.54588877 0.92799 1.0000000000
#> 121  -2.56291555 0.99449 1.0000000000
#> 141  -2.15466921 0.98321 1.0000000000
#> 151  -2.72659766 0.99661 1.0000000000
#> 161  -1.68864663 0.94735 1.0000000000
#> 171  -2.47592155 0.99303 1.0000000000
#> 181  -0.31726079 0.59794 1.0000000000
#> 191  -2.01636582 0.97642 1.0000000000
#> 201   0.80025753 0.20243 0.9473724000
#> 241  -3.81710920 0.99995 1.0000000000
#> 251  -0.66930097 0.73200 1.0000000000
#> 271   3.04942490 0.00104 0.0135200000
#> 281  -3.28086388 0.99947 1.0000000000
#> 291  -4.49277203 0.99999 1.0000000000
#> 311  -1.59189811 0.93701 1.0000000000
#> 321  -3.65739585 0.99989 1.0000000000
#> 341   0.91163582 0.14460 0.8459100000
#> 361  -6.37499749 1.00000 1.0000000000
#> 381  -5.60520864 1.00000 1.0000000000
#> 391   6.13231305 0.00001 0.0001950000
#> 401  -4.51194107 1.00000 1.0000000000
#> 411  -7.32250958 1.00000 1.0000000000
#> 431  -4.07109644 0.99999 1.0000000000
#> 451  -4.30056027 0.99999 1.0000000000
#> 471  -4.12855302 0.99997 1.0000000000
#> 481  -5.09619879 1.00000 1.0000000000
#> 491  -2.66802436 0.99580 1.0000000000
#> 501  -2.74414432 0.99669 1.0000000000
#> 521  -3.39268262 0.99964 1.0000000000
#> 531  -5.01503528 1.00000 1.0000000000
#> 541  -4.22391595 1.00000 1.0000000000
#> 551  -7.07913203 1.00000 1.0000000000
#> 581  -5.50243944 1.00000 1.0000000000
#> 591  -0.32815979 0.56521 1.0000000000
#> 601   3.93150186 0.00005 0.0007312500
#> 611   4.03583179 0.00003 0.0005014286
#> 621  -0.77678465 0.76634 1.0000000000
#> 631  -5.07662005 1.00000 1.0000000000
#> 651  -1.31376133 0.89481 1.0000000000
#> 671   5.01657755 0.00001 0.0001950000
#> 681   4.33420937 0.00001 0.0001950000
#> 691  -5.73342028 1.00000 1.0000000000
#> 701  -1.78710458 0.95993 1.0000000000
#> 711  -5.15649865 1.00000 1.0000000000
#> 721  -0.50075903 0.67561 1.0000000000
#> 741  -4.40566447 1.00000 1.0000000000
#> 751  -2.03735548 0.97732 1.0000000000
#> 771  -5.67326173 1.00000 1.0000000000
#> 781  -4.06514412 0.99993 1.0000000000
#> 791  -0.05993985 0.37263 1.0000000000
#> 801  -2.53354365 0.99372 1.0000000000
#> 811   0.73175587 0.17990 0.8770125000
#> 821  -1.70488473 0.95282 1.0000000000
#> 841   5.35971692 0.00001 0.0001950000
#> 861  -2.48012387 0.99295 1.0000000000
#> 871  -2.33911838 0.98983 1.0000000000
#> 88   -3.13104349 0.99894 1.0000000000
#> 89   -2.51848970 0.99383 1.0000000000
#> 90   -3.54571085 0.99979 1.0000000000
#> 93   -1.93970751 0.97149 1.0000000000
#> 96   -0.78778123 0.77040 1.0000000000
#> 97   -1.71437386 0.95190 1.0000000000
#> 98   -1.60278311 0.92252 1.0000000000
#> 99    2.34556842 0.00838 0.0980460000
#> 103  -0.75911537 0.73913 1.0000000000
#> 104   0.35674397 0.32904 1.0000000000
#> 109  -4.18557512 0.99999 1.0000000000
#> 1111 -3.50819931 0.99974 1.0000000000
#> 112  -2.08987773 0.98017 1.0000000000
#> 115  -1.81424803 0.96234 1.0000000000
#> 116  -3.73822398 0.99993 1.0000000000
#> 117  -2.12238400 0.98216 1.0000000000
#> 
#> $exclusive_up
#>             RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 1        ALKBH1            6                5.37530            2.2403482
#> 2        ALKBH5          113               84.71910            6.5891327
#> 3         ATXN2          180              163.59850            4.0062100
#> 4      C17orf85          122              105.98602            6.5137784
#> 5       CAPRIN1          170              132.86828            5.8843025
#> 6         CELF2          126               88.24509            6.6218834
#> 7          CNBP           41               39.95718            5.4678410
#> 8         CPSF1          156              131.05240            5.9574566
#> 9         CPSF2          106               74.43014            6.5107775
#> 10        CPSF3          165              139.10598            5.6239083
#> 11        CPSF4          167              129.92143            5.9920794
#> 12        CPSF6          176              138.58517            5.6538470
#> 13        CPSF7          177              137.42694            5.7001603
#> 14        CSTF2          180              156.33098            4.6193880
#> 15       CSTF2T          179              157.83387            4.4990125
#> 16        DDX3X          179              161.95023            4.1500488
#> 17        DGCR8           34               34.02309            5.1757246
#> 18       DICER1           55               56.37151            6.1333031
#> 19       DIS3L2           21               16.06965            3.7504319
#> 20        EIF3A           79               70.37887            6.4435986
#> 21        EIF3B           84               78.39929            6.5715728
#> 22        EIF3D           88               85.54499            6.6024247
#> 23        EIF3G          119              109.63739            6.4645094
#> 24       EIF4A3           76               69.24117            6.4299580
#> 25       ELAVL1          182              168.80798            3.4304148
#> 26        EWSR1          164              115.55597            6.3733016
#> 27         EZH2            3                0.66151            0.7952615
#> 28          FBL           45               41.89964            5.5763912
#> 29       FIP1L1          182              163.95899            3.9631787
#> 30         FMR1          145              119.33305            6.2884309
#> 31          FUS          178              152.38571            4.8795795
#> 32         FXR1           55               26.11463            4.6348295
#> 33         FXR2          107               68.97834            6.4088206
#> 34      HNRNPA1          176              136.98222            5.7173605
#> 35    HNRNPA2B1           24               14.73707            3.6136171
#> 36       HNRNPC          182              160.69826            4.2612715
#> 37       HNRNPD          162              105.46488            6.5257837
#> 38       HNRNPF           57               32.16059            5.0494194
#> 39       HNRNPH           25               26.61086            4.6510586
#> 40       HNRNPM           54               34.60214            5.1817966
#> 41       HNRNPU           99               67.62095            6.3806472
#> 42      IGF2BP1          163              116.83672            6.3690663
#> 43      IGF2BP2          167              111.73327            6.4489705
#> 44      IGF2BP3          175              127.35236            6.0897663
#> 45       LIN28A          149              122.60024            6.1930134
#> 46       LIN28B          180              163.84451            3.9808657
#> 47        MBNL2           67               49.82545            5.9068444
#> 48      METTL14           18               12.52510            3.3459620
#> 49       METTL3           29               23.18325            4.4245593
#> 50        MOV10          174              139.14097            5.6278179
#> 51        NOP56           51               37.34404            5.3389195
#> 52        NOP58           65               39.12373            5.4626275
#> 53       NUDT21          179              161.80582            4.1687777
#> 54        PRKRA            5                8.51499            2.8127521
#> 55        PTBP1          180              159.67048            4.3448758
#> 56   PTBP1PTBP2          105               87.03934            6.6077446
#> 57         PUM2           86               55.45614            6.0755399
#> 58          QKI           88               53.41721            6.0283803
#> 59        RBM10           26               26.61460            4.6744310
#> 60       RBM15B          172              141.03569            5.5223891
#> 61        RBM15          135              113.21712            6.4094282
#> 62        RBM47          100               62.04289            6.2593005
#> 63        RBPMS           61               38.96723            5.4220983
#> 64         RTCB          147              116.45607            6.3863995
#> 65        SRRM4          130               95.21804            6.6135676
#> 66        SRSF1          157              127.95291            6.0583067
#> 67          SSB           53               48.83378            5.8727145
#> 68        STAU1          100               73.83200            6.5020473
#> 69        TAF15          132               90.80896            6.5895400
#> 70       TARBP2           46               52.82026            6.0233609
#> 71       TARDBP          180              168.65969            3.4501476
#> 72         TIA1          119               80.90386            6.5794050
#> 73        TIAL1          153              114.15064            6.3653258
#> 74       TNRC6A            7                6.54626            2.4619059
#> 75       TNRC6B            5                5.24066            2.2070866
#> 76       TNRC6C            7                4.51571            2.0540339
#> 77        U2AF2          177              153.78678            4.7935485
#> 78         UPF1          170              136.08271            5.7515563
#> 79        WDR33          178              138.43263            5.6609789
#> 80         WTAP           18               12.73928            3.3855103
#> 81       YTHDC1          169              138.31190            5.6434136
#> 82       YTHDC2           46               37.85955            5.3648403
#> 83       YTHDF1          176              146.54885            5.2442052
#> 84       YTHDF2          172              152.70817            4.8629766
#> 85       YTHDF3           90               87.72570            6.6164974
#> 86       ZC3H7B          174              128.06706            6.0618620
#> 87      ZC3H11A           36               37.00524            5.3297220
#> 510       DDX24           54               53.30842            6.0300896
#> 610       XRCC6           12               19.23606            4.0710811
#> 710       SF3B4           80               88.90964            6.6234125
#> 91        NCBP2           51               49.37618            5.8889588
#> 101     FAM120A           77               66.38116            6.3930451
#> 111       LARP7            4                4.78706            2.1194201
#> 121         FTO            9               12.51361            3.3404590
#> 141       RBM22           42               51.84665            5.9814690
#> 151       UCHL5           78               67.26804            6.3839832
#> 161       DDX42            8                6.87948            2.5276785
#> 171        XPO5           25               32.42602            5.0669245
#> 181      TROVE2            8               11.90820            3.2600152
#> 191        SUB1           43               36.91274            5.3031657
#> 201       TBRG4           47               53.78239            6.0308274
#> 241      METAP2            7               19.25905            4.0710968
#> 251        DKC1           11               21.97705            4.3189177
#> 271     KHDRBS1            7                7.36004            2.6031095
#> 281        SLTM           10               22.67154            4.3701058
#> 291        SND1           25               36.71909            5.3084707
#> 311      SERBP1           10                8.39074            2.7779920
#> 321        NONO           15               23.41061            4.4333312
#> 341       NSUN2            1                1.83105            1.3161699
#> 361        PPIG           61               80.79242            6.5622908
#> 381      AKAP8L            6               17.08318            3.8598251
#> 391      HNRNPL           60               47.45669            5.7977634
#> 401       PCBP2           62               73.26724            6.4806344
#> 411       SF3A3           30               54.59775            6.0563002
#> 431      RBFOX2           39               51.74552            5.9788977
#> 451      HNRNPK           19               38.02223            5.3841811
#> 471      GEMIN5           55               50.98269            5.9545499
#> 481       NOL12            5               13.21219            3.4318542
#> 491     SUPV3L1            2                9.73388            2.9881214
#> 501       SRSF7           11               13.22244            3.4275353
#> 521      GTF2F1           23               45.40975            5.7307716
#> 531      EFTUD2           80               92.15475            6.6123248
#> 541       PRPF8          135              125.35744            6.1096637
#> 551       DDX55           35               49.07950            5.8720716
#> 581       BUD13           90              100.05423            6.5807661
#> 591       YWHAG            0                1.94047            1.3552655
#> 601       FUBP3           49               28.84618            4.8327769
#> 611       MATR3           45               37.23776            5.3374296
#> 621       TRA2A           23               18.17279            3.9717794
#> 631      DROSHA           35               40.83926            5.5369711
#> 651       SRSF9            6               10.57914            3.0764641
#> 671       KHSRP          107               61.60886            6.2371098
#> 681      EXOSC5           18               19.94338            4.1188984
#> 691       LARP4           32               44.53398            5.6565577
#> 701       CDC40            5               13.68347            3.4845086
#> 711       LSM11           14               29.24427            4.8655646
#> 721     FASTKD2           16               27.19096            4.7154021
#> 741      SMNDC1           13               24.20415            4.5077285
#> 751       DHX30            4               15.17076            3.6659508
#> 771        RPS3           41               58.80168            6.1877275
#> 781        XRN2           41               51.40128            5.9878822
#> 791        SBDS            0                0.35098            0.5790824
#> 801        RBM5            5               15.58681            3.7052667
#> 811        NPM1            5                1.35747            1.1328539
#> 821      ZNF622           72               60.69684            6.2255282
#> 841       SUGP2           27               30.92755            4.9660978
#> 861       NOLC1           14               22.82992            4.3853626
#> 871        YBX3           83               67.79643            6.4228982
#> 88        SF3B1          109               98.07803            6.6067751
#> 89        DDX59            9               19.98095            4.1275801
#> 90        GRSF1            9               11.35489            3.2039703
#> 93         NKRF           18               27.51722            4.7646879
#> 96       ZRANB2           17               21.21564            4.2643946
#> 97        MTPAP            8               11.19070            3.1780551
#> 98         PUS1            0                1.01346            0.9840572
#> 99        SAFB2            5                8.31952            2.7655205
#> 103       BCCIP            4                3.56698            1.8327540
#> 104    HNRNPUL1            4                7.82141            2.6866611
#> 109        AARS           18               22.76861            4.3783490
#> 1111       DDX6           38               50.20626            5.9338157
#> 112      EIF4G2           23               23.88773            4.4835105
#> 115       EIF3H           18               28.04853            4.7738227
#> 116       FKBP4            3                8.36430            2.7624775
#> 117       GRWD1           94               79.45923            6.5709611
#>                 z    pval         padj
#> 1     0.278840583 0.29304 3.360192e-01
#> 2     4.292051991 0.00001 1.622642e-05
#> 3     4.094019076 0.00001 1.622642e-05
#> 4     2.458477868 0.00548 7.725902e-03
#> 5     6.310300992 0.00001 1.622642e-05
#> 6     5.701536539 0.00001 1.622642e-05
#> 7     0.190718786 0.38229 4.161638e-01
#> 8     4.187625974 0.00001 1.622642e-05
#> 9     4.848861763 0.00001 1.622642e-05
#> 10    4.604274934 0.00001 1.622642e-05
#> 11    6.187930330 0.00001 1.622642e-05
#> 12    6.617588007 0.00001 1.622642e-05
#> 13    6.942446895 0.00001 1.622642e-05
#> 14    5.123843229 0.00001 1.622642e-05
#> 15    4.704616874 0.00001 1.622642e-05
#> 16    4.108329987 0.00001 1.622642e-05
#> 17   -0.004461211 0.45637 4.845410e-01
#> 18   -0.223616864 0.55243 5.723973e-01
#> 19    1.314608588 0.07724 9.914388e-02
#> 20    1.337937161 0.07900 9.966171e-02
#> 21    0.852263248 0.17545 2.095653e-01
#> 22    0.371834608 0.32714 3.695766e-01
#> 23    1.448309440 0.06362 8.289879e-02
#> 24    1.051146827 0.13008 1.575617e-01
#> 25    3.845604862 0.00001 1.622642e-05
#> 26    7.601088554 0.00001 1.622642e-05
#> 27    2.940529646 0.00363 5.291186e-03
#> 28    0.555979645 0.25622 2.977692e-01
#> 29    4.552156604 0.00001 1.622642e-05
#> 30    4.081614384 0.00001 1.622642e-05
#> 31    5.249282260 0.00001 1.622642e-05
#> 32    6.232240001 0.00001 1.622642e-05
#> 33    5.932707799 0.00001 1.622642e-05
#> 34    6.824439389 0.00001 1.622642e-05
#> 35    2.563340220 0.00592 8.211613e-03
#> 36    4.998916436 0.00001 1.622642e-05
#> 37    8.663345722 0.00001 1.622642e-05
#> 38    4.919260663 0.00001 1.622642e-05
#> 39   -0.346342656 0.58644 6.004029e-01
#> 40    3.743462263 0.00012 1.876364e-04
#> 41    4.917847547 0.00001 1.622642e-05
#> 42    7.248045155 0.00001 1.622642e-05
#> 43    8.569853076 0.00001 1.622642e-05
#> 44    7.824214896 0.00001 1.622642e-05
#> 45    4.262829437 0.00001 1.622642e-05
#> 46    4.058285663 0.00001 1.622642e-05
#> 47    2.907567684 0.00172 2.595088e-03
#> 48    1.636270807 0.04306 5.878032e-02
#> 49    1.314650697 0.08058 9.966171e-02
#> 50    6.194057927 0.00001 1.622642e-05
#> 51    2.557813426 0.00516 7.396000e-03
#> 52    4.736964027 0.00001 1.622642e-05
#> 53    4.124513476 0.00001 1.622642e-05
#> 54   -1.249662205 0.86064 8.606400e-01
#> 55    4.678964632 0.00001 1.622642e-05
#> 56    2.718122601 0.00254 3.766207e-03
#> 57    5.027349069 0.00001 1.622642e-05
#> 58    5.736663613 0.00001 1.622642e-05
#> 59   -0.131481243 0.50009 5.244846e-01
#> 60    5.607049690 0.00001 1.622642e-05
#> 61    3.398568383 0.00018 2.764286e-04
#> 62    6.064113753 0.00001 1.622642e-05
#> 63    4.063513556 0.00001 1.622642e-05
#> 64    4.782652606 0.00001 1.622642e-05
#> 65    5.259182681 0.00001 1.622642e-05
#> 66    4.794588924 0.00001 1.622642e-05
#> 67    0.709419807 0.21077 2.483044e-01
#> 68    4.024578551 0.00003 4.777778e-05
#> 69    6.250973475 0.00001 1.622642e-05
#> 70   -1.132301406 0.85476 8.606400e-01
#> 71    3.286905799 0.00001 1.622642e-05
#> 72    5.790210524 0.00001 1.622642e-05
#> 73    6.103279096 0.00001 1.622642e-05
#> 74    0.184304364 0.33090 3.695766e-01
#> 75   -0.109039670 0.42761 4.596807e-01
#> 76    1.209468823 0.08112 9.966171e-02
#> 77    4.842596274 0.00001 1.622642e-05
#> 78    5.897063027 0.00001 1.622642e-05
#> 79    6.989492616 0.00001 1.622642e-05
#> 80    1.553892796 0.04973 6.682469e-02
#> 81    5.437861243 0.00001 1.622642e-05
#> 82    1.517370419 0.05546 7.337785e-02
#> 83    5.615941520 0.00001 1.622642e-05
#> 84    3.967082603 0.00001 1.622642e-05
#> 85    0.343731717 0.33683 3.713767e-01
#> 86    7.577364893 0.00001 1.622642e-05
#> 87   -0.188610213 0.53025 9.999800e-01
#> 510   0.114688180 0.41740 9.999800e-01
#> 610  -1.777429583 0.95686 9.999800e-01
#> 710  -1.345173645 0.89750 9.999800e-01
#> 91    0.275739747 0.35528 9.999800e-01
#> 101   1.660998760 0.04180 2.574000e-01
#> 111  -0.371356300 0.52665 9.999800e-01
#> 121  -1.051834492 0.81501 9.999800e-01
#> 141  -1.646192596 0.94399 9.999800e-01
#> 151   1.681075840 0.04007 2.574000e-01
#> 161   0.443300040 0.24879 9.096384e-01
#> 171  -1.465587267 0.91773 9.999800e-01
#> 181  -1.198828770 0.85492 9.999800e-01
#> 191   1.147853997 0.10824 5.065632e-01
#> 201  -1.124620140 0.85074 9.999800e-01
#> 241  -3.011240145 0.99939 9.999800e-01
#> 251  -2.541620573 0.99569 9.999800e-01
#> 271  -0.138311509 0.45707 9.999800e-01
#> 281  -2.899595680 0.99883 9.999800e-01
#> 291  -2.207620729 0.98508 9.999800e-01
#> 311   0.579288930 0.21726 8.679060e-01
#> 321  -1.897130967 0.96817 9.999800e-01
#> 341  -0.631415440 0.55317 9.999800e-01
#> 361  -3.016083934 0.99855 9.999800e-01
#> 381  -2.871420236 0.99891 9.999800e-01
#> 391   2.163473937 0.01340 1.406700e-01
#> 401  -1.738601382 0.95209 9.999800e-01
#> 411  -4.061514343 0.99998 9.999800e-01
#> 431  -2.131750817 0.98158 9.999800e-01
#> 451  -3.532984827 0.99995 9.999800e-01
#> 471   0.674662240 0.22254 8.679060e-01
#> 481  -2.392930919 0.99317 9.999800e-01
#> 491  -2.588208106 0.99746 9.999800e-01
#> 501  -0.648407627 0.68273 9.999800e-01
#> 521  -3.910424545 0.99995 9.999800e-01
#> 531  -1.838196143 0.96036 9.999800e-01
#> 541   1.578247274 0.04563 2.669355e-01
#> 551  -2.397705788 0.99121 9.999800e-01
#> 581  -1.527820603 0.92726 9.999800e-01
#> 591  -1.431800644 0.86261 9.999800e-01
#> 601   4.170235921 0.00004 1.170000e-03
#> 611   1.454303014 0.06237 3.316950e-01
#> 621   1.215377177 0.09314 4.737991e-01
#> 631  -1.054594627 0.83377 9.999800e-01
#> 651  -1.488442516 0.91503 9.999800e-01
#> 671   7.277591952 0.00001 3.900000e-04
#> 681  -0.471820332 0.62565 9.999800e-01
#> 691  -2.215831713 0.98577 9.999800e-01
#> 701  -2.492021401 0.99530 9.999800e-01
#> 711  -3.133093740 0.99939 9.999800e-01
#> 721  -2.373278004 0.99170 9.999800e-01
#> 741  -2.485542346 0.99450 9.999800e-01
#> 751  -3.047165810 0.99950 9.999800e-01
#> 771  -2.876933397 0.99775 9.999800e-01
#> 781  -1.737054889 0.95318 9.999800e-01
#> 791  -0.606096848 0.30152 9.999800e-01
#> 801  -2.857232922 0.99880 9.999800e-01
#> 811   3.215357157 0.00178 2.603250e-02
#> 821   1.815614612 0.02937 2.315820e-01
#> 841  -0.790872459 0.75259 9.999800e-01
#> 861  -2.013498269 0.97625 9.999800e-01
#> 871   2.367088724 0.00773 1.004900e-01
#> 88    1.653146923 0.04076 2.574000e-01
#> 89   -2.660384472 0.99718 9.999800e-01
#> 90   -0.734991206 0.70862 9.999800e-01
#> 93   -1.997448786 0.97507 9.999800e-01
#> 96   -0.988567054 0.80606 9.999800e-01
#> 97   -1.003978803 0.79860 9.999800e-01
#> 98   -1.029879190 0.64667 9.999800e-01
#> 99   -1.200323787 0.84735 9.999800e-01
#> 103   0.236267386 0.28401 9.773285e-01
#> 104  -1.422363996 0.89940 9.999800e-01
#> 109  -1.089134272 0.83406 9.999800e-01
#> 1111 -2.057067592 0.97823 9.999800e-01
#> 112  -0.197998867 0.52337 9.999800e-01
#> 115  -2.104923188 0.98133 9.999800e-01
#> 116  -1.941843885 0.97248 9.999800e-01
#> 117   2.212883292 0.01137 1.330290e-01
#> 
#> $exclusive_down
#>             RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 1        ALKBH1            7                8.33111            2.7625735
#> 2        ALKBH5          116              131.22288            8.1199288
#> 3         ATXN2          261              253.47956            4.9162632
#> 4      C17orf85          186              164.15552            8.0227587
#> 5       CAPRIN1          179              205.89282            7.2651800
#> 6         CELF2          112              136.72033            8.1271013
#> 7          CNBP           38               61.94491            6.7676978
#> 8         CPSF1          202              203.06069            7.3071404
#> 9         CPSF2           94              115.34442            7.9683656
#> 10        CPSF3          202              215.53211            6.9314883
#> 11        CPSF4          189              201.33701            7.3809646
#> 12        CPSF6           17               13.67988            3.4873292
#> 13        CPSF7          193              212.94062            7.0200817
#> 14        CSTF2          235              242.20847            5.6724221
#> 15       CSTF2T          114               75.61940            7.1829938
#> 16        DDX3X          220              182.22444            7.7613935
#> 17        DGCR8          128               89.79479            7.5545807
#> 18       DICER1           95               87.33473            7.5102177
#> 19       DIS3L2           22               24.89057            4.6133489
#> 20        EIF3A          114              109.05015            7.9252081
#> 21        EIF3B          118              121.43967            8.0566202
#> 22        EIF3D          171              132.52318            8.1329370
#> 23        EIF3G          180              169.86184            8.0149457
#> 24       EIF4A3          130              107.30673            7.9208051
#> 25       ELAVL1          251              261.57757            4.2227102
#> 26        EWSR1           57               25.45226            4.6502341
#> 27         EZH2            1                1.02084            0.9795689
#> 28          FBL           64               64.92047            6.8681174
#> 29       FIP1L1          256              254.03266            4.8705103
#> 30         FMR1           74               56.97499            6.5406049
#> 31          FUS          222              236.15512            6.0296236
#> 32         FXR1           17               11.43492            3.2198988
#> 33         FXR2           65               45.74684            5.9989391
#> 34      HNRNPA1            1                3.73562            1.8603166
#> 35    HNRNPA2B1           15               22.83487            4.4473925
#> 36       HNRNPC           51               48.34004            6.1213779
#> 37       HNRNPD          104              163.42362            8.0588284
#> 38       HNRNPF           22               49.85540            6.2142978
#> 39       HNRNPH           40               41.24649            5.7497916
#> 40       HNRNPM           77               80.43285            7.3230135
#> 41       HNRNPU           10                9.46408            2.9426444
#> 42      IGF2BP1           82               88.94059            7.5197916
#> 43      IGF2BP2           52               51.79202            6.3134636
#> 44      IGF2BP3           26               22.94837            4.4514764
#> 45       LIN28A          183              189.97008            7.6315298
#> 46       LIN28B          119               98.21578            7.7607384
#> 47        MBNL2           66               77.23079            7.2689940
#> 48      METTL14           19               19.40925            4.1295442
#> 49       METTL3           40               35.90435            5.4374200
#> 50        MOV10          184              215.60573            6.9327932
#> 51        NOP56           53               57.84435            6.5937484
#> 52        NOP58           38               60.63831            6.6960525
#> 53       NUDT21          262              250.71616            5.1258890
#> 54        PRKRA           17               13.19908            3.4552462
#> 55        PTBP1          115              107.40036            7.8605617
#> 56   PTBP1PTBP2          121              134.83430            8.1310654
#> 57         PUM2           62               85.93277            7.4688920
#> 58          QKI           91               95.39361            7.7018851
#> 59        RBM10           30               41.19457            5.7609656
#> 60       RBM15B          224              218.51532            6.8509601
#> 61        RBM15          112               96.60761            7.7072559
#> 62        RBM47           38               96.11678            7.7174114
#> 63        RBPMS           52               60.41977            6.6740699
#> 64         RTCB          180              180.46063            7.8350395
#> 65        SRRM4          156              147.52691            8.1241674
#> 66        SRSF1           68               53.10104            6.3492956
#> 67          SSB           87               75.68203            7.2102708
#> 68        STAU1          123              114.35618            8.0340240
#> 69        TAF15           18               15.76416            3.7513998
#> 70       TARBP2           69               81.77966            7.3835828
#> 71       TARDBP          264              261.32311            4.2499212
#> 72         TIA1           80               85.28731            7.4839390
#> 73        TIAL1          140              176.89857            7.8873953
#> 74       TNRC6A           13               10.12716            3.0367190
#> 75       TNRC6B           12                8.14371            2.7267841
#> 76       TNRC6C            9                7.00273            2.5356748
#> 77        U2AF2          241              238.31387            5.8931166
#> 78         UPF1           98               78.83323            7.2802505
#> 79        WDR33          187              214.53616            6.9295969
#> 80         WTAP           20               19.72253            4.1412162
#> 81       YTHDC1          230              214.31582            6.9934774
#> 82       YTHDC2           67               58.70538            6.6294629
#> 83       YTHDF1          238              227.09690            6.4628189
#> 84       YTHDF2          251              236.63017            5.9952677
#> 85       YTHDF3          164              135.92766            8.1123465
#> 86       ZC3H7B          141              198.43642            7.4762956
#> 87      ZC3H11A           81               57.33332            6.5298012
#> 510       DDX24          121               82.58047            7.4072217
#> 610       XRCC6           46               29.75390            4.9714064
#> 710       SF3B4          203              137.77681            8.1236209
#> 91        NCBP2           73               76.56581            7.2421153
#> 101     FAM120A          121              102.78765            7.8188637
#> 111       LARP7           13                7.42051            2.6124451
#> 121         FTO           32               19.41221            4.1374224
#> 141       RBM22          115               80.33876            7.3611564
#> 151       UCHL5          114              104.16633            7.8300062
#> 161       DDX42           14               10.66942            3.1023013
#> 171        XPO5           81               50.19892            6.2268177
#> 181      TROVE2           30               18.43307            4.0269619
#> 191        SUB1           76               57.19697            6.5263327
#> 201       TBRG4          103               83.23211            7.4212496
#> 241      METAP2           67               29.83893            5.0235226
#> 251        DKC1           50               34.03692            5.2936865
#> 271     KHDRBS1            9               11.39091            3.1888087
#> 281        SLTM           55               35.12678            5.3618312
#> 291        SND1           94               56.88892            6.5421426
#> 311      SERBP1           31               13.02977            3.4048935
#> 321        NONO           60               36.28751            5.4184353
#> 341       NSUN2            4                2.83650            1.6229400
#> 361        PPIG          181              125.14984            8.0928032
#> 381      AKAP8L           55               26.47226            4.7394172
#> 391      HNRNPL           61               73.51455            7.1395769
#> 401       PCBP2          167              113.52356            7.9952739
#> 411       SF3A3          136               84.53449            7.4593878
#> 431      RBFOX2          138               80.10785            7.3027674
#> 451      HNRNPK          110               58.88217            6.6136165
#> 471      GEMIN5           99               78.99584            7.2878154
#> 481       NOL12           49               20.43373            4.2252677
#> 491     SUPV3L1           35               15.06406            3.6647470
#> 501       SRSF7           38               20.48714            4.1998989
#> 521      GTF2F1          117               70.32837            7.0327034
#> 531      EFTUD2          188              142.68451            8.1235593
#> 541       PRPF8          228              194.18893            7.5118466
#> 551       DDX55          126               75.97415            7.2101180
#> 581       BUD13          195              154.95469            8.0771845
#> 591       YWHAG           11                3.00633            1.6752486
#> 601       FUBP3           23               44.69288            5.9640081
#> 611       MATR3           48               57.65964            6.5606284
#> 621       TRA2A           28               28.09618            4.8851559
#> 631      DROSHA          113               63.27024            6.7828600
#> 651       SRSF9           29               16.38688            3.8216083
#> 671       KHSRP           52               95.44124            7.7003169
#> 681      EXOSC5           21               30.89772            5.0822138
#> 691       LARP4          124               68.91859            6.9865893
#> 701       CDC40           41               21.19824            4.2942572
#> 711       LSM11           81               45.22624            5.9775156
#> 721     FASTKD2           70               42.09784            5.7925783
#> 741      SMNDC1           65               37.48859            5.5121224
#> 751       DHX30           42               23.50187            4.5062600
#> 771        RPS3          144               91.10023            7.6184634
#> 781        XRN2          112               79.58863            7.3051583
#> 791        SBDS            2                0.54249            0.7148564
#> 801        RBM5           44               24.14394            4.5555448
#> 811        NPM1            2                2.10318            1.4035860
#> 821      ZNF622           85               93.97778            7.6536378
#> 841       SUGP2           30               47.88570            6.1059945
#> 861       NOLC1           63               35.37240            5.4083834
#> 871        YBX3          124              105.01701            7.8567766
#> 88        SF3B1          176              151.93262            8.1008022
#> 89        DDX59           69               30.97425            5.0868208
#> 90        GRSF1           34               17.57802            3.9353663
#> 93         NKRF           68               42.58728            5.8074142
#> 96       ZRANB2           55               32.84997            5.2092411
#> 97        MTPAP           29               17.34499            3.8988106
#> 98         PUS1            5                1.56824            1.2111969
#> 99        SAFB2           10               12.89409            3.3976887
#> 103       BCCIP           10                5.53045            2.2554963
#> 104    HNRNPUL1           12               12.10939            3.3020710
#> 109        AARS           52               35.27258            5.3772790
#> 1111       DDX6          101               77.72861            7.2707693
#> 112      EIF4G2           53               37.01162            5.4900115
#> 115       EIF3H           76               43.43814            5.8645594
#> 116       FKBP4           31               12.95467            3.4079218
#> 117       GRWD1          129              123.05038            8.0498838
#>                z    pval         padj
#> 1    -0.48183695 0.59950 1.000000e+00
#> 2    -1.87475535 0.96556 1.000000e+00
#> 3     1.52970654 0.04587 3.944820e-01
#> 4     2.72281404 0.00276 4.747200e-02
#> 5    -3.70160410 0.99984 1.000000e+00
#> 6    -3.04171551 0.99875 1.000000e+00
#> 7    -3.53811751 0.99987 1.000000e+00
#> 8    -0.14515802 0.53335 1.000000e+00
#> 9    -2.67864468 0.99600 1.000000e+00
#> 10   -1.95226617 0.96745 1.000000e+00
#> 11   -1.67146309 0.94387 1.000000e+00
#> 12    0.95205237 0.13865 2.053424e-01
#> 13   -2.84051111 0.99646 1.000000e+00
#> 14   -1.27079225 0.88043 1.000000e+00
#> 15    5.34325954 0.00001 3.078947e-05
#> 16    4.86711053 0.00001 3.078947e-05
#> 17    5.05722442 0.00001 3.078947e-05
#> 18    1.02064551 0.13841 5.668219e-01
#> 19   -0.62656653 0.69233 1.000000e+00
#> 20    0.62457035 0.24446 8.085985e-01
#> 21   -0.42693709 0.64171 1.000000e+00
#> 22    4.73098710 0.00001 8.600000e-04
#> 23    1.26490689 0.09160 5.668219e-01
#> 24    2.86502061 0.00191 4.106500e-02
#> 25   -2.50492446 0.98874 1.000000e+00
#> 26    6.78411868 0.00001 3.078947e-05
#> 27   -0.02127466 0.27157 8.650007e-01
#> 28   -0.13402071 0.51969 1.000000e+00
#> 29    0.40392893 0.31194 9.250634e-01
#> 30    2.60297179 0.00442 8.080313e-03
#> 31   -2.34759594 0.98596 1.000000e+00
#> 32    1.72834003 0.03603 5.854875e-02
#> 33    3.20942750 0.00076 1.778400e-03
#> 34   -1.47051314 0.89635 9.710458e-01
#> 35   -1.76167722 0.95627 1.000000e+00
#> 36    0.43453615 0.29911 3.976803e-01
#> 37   -7.37372939 1.00000 1.000000e+00
#> 38   -4.48246945 1.00000 1.000000e+00
#> 39   -0.21678873 0.54476 1.000000e+00
#> 40   -0.46877559 0.65267 7.814891e-01
#> 41    0.18212190 0.34786 4.472486e-01
#> 42   -0.92297637 0.80301 8.947826e-01
#> 43    0.03294230 0.44879 5.527203e-01
#> 44    0.68553210 0.20938 2.952300e-01
#> 45   -0.91332671 0.80192 1.000000e+00
#> 46    2.67812401 0.00371 7.115902e-03
#> 47   -1.54502672 0.93120 1.000000e+00
#> 48   -0.09910295 0.47634 1.000000e+00
#> 49    0.75323406 0.19710 7.062750e-01
#> 50   -4.55887389 0.99999 1.000000e+00
#> 51   -0.73468833 0.74226 1.000000e+00
#> 52   -3.38084418 0.99975 1.000000e+00
#> 53    2.20134304 0.00788 9.040750e-02
#> 54    1.10004316 0.10870 5.668219e-01
#> 55    0.96680623 0.15130 2.185444e-01
#> 56   -1.70141295 0.94934 1.000000e+00
#> 57   -3.20432669 0.99935 1.000000e+00
#> 58   -0.57045904 0.69138 8.089146e-01
#> 59   -1.94317599 0.97160 1.000000e+00
#> 60    0.80057100 0.19143 7.062750e-01
#> 61    1.99712973 0.02041 3.411386e-02
#> 62   -7.53060545 1.00000 1.000000e+00
#> 63   -1.26156456 0.88286 1.000000e+00
#> 64   -0.05879102 0.50189 1.000000e+00
#> 65    1.04294872 0.13488 5.668219e-01
#> 66    2.34655323 0.00866 1.534970e-02
#> 67    1.56970110 0.05219 4.080309e-01
#> 68    1.07590169 0.12724 5.668219e-01
#> 69    0.59600153 0.22815 3.103901e-01
#> 70   -1.73082097 0.95362 1.000000e+00
#> 71    0.62986814 0.23325 8.023800e-01
#> 72   -0.70648759 0.73815 8.384811e-01
#> 73   -4.67816925 1.00000 1.000000e+00
#> 74    0.94603418 0.13418 5.668219e-01
#> 75    1.41422635 0.06211 4.451217e-01
#> 76    0.78766803 0.16055 6.276045e-01
#> 77    0.45580805 0.29755 9.139036e-01
#> 78    2.63270748 0.00402 7.465714e-03
#> 79   -3.97370302 0.99990 1.000000e+00
#> 80    0.06700206 0.41394 1.000000e+00
#> 81    2.24268686 0.00841 9.040750e-02
#> 82    1.25117526 0.09336 5.668219e-01
#> 83    1.68705021 0.03583 3.423756e-01
#> 84    2.39686210 0.00456 6.536000e-02
#> 85    3.46044639 0.00025 1.075000e-02
#> 86   -7.68247040 1.00000 1.000000e+00
#> 87    3.62441048 0.00023 5.850000e-04
#> 510   5.18676659 0.00001 3.078947e-05
#> 610   3.26790823 0.00061 1.486875e-03
#> 710   8.02883234 0.00001 3.078947e-05
#> 91   -0.49237134 0.66126 7.814891e-01
#> 101   2.32928349 0.00879 1.534970e-02
#> 111   2.13573480 0.01602 2.716435e-02
#> 121   3.04242323 0.00151 3.333396e-03
#> 141   4.70866779 0.00001 3.078947e-05
#> 151   1.25589555 0.09437 1.452801e-01
#> 161   1.07358368 0.11184 1.699387e-01
#> 171   4.94652028 0.00001 3.078947e-05
#> 181   2.87237132 0.00231 4.659828e-03
#> 191   2.88110196 0.00215 4.413158e-03
#> 201   2.66368753 0.00325 6.337500e-03
#> 241   7.39741280 0.00001 3.078947e-05
#> 251   3.01549400 0.00161 3.488333e-03
#> 271  -0.74978157 0.71353 8.184609e-01
#> 281   3.70642400 0.00018 4.680000e-04
#> 291   5.67261867 0.00001 3.078947e-05
#> 311   5.27776568 0.00001 3.078947e-05
#> 321   4.37626150 0.00002 5.707317e-05
#> 341   0.71690884 0.15115 2.185444e-01
#> 361   6.90121320 0.00001 3.078947e-05
#> 381   6.01925067 0.00001 3.078947e-05
#> 391  -1.75284197 0.95547 1.000000e+00
#> 401   6.68850635 0.00001 3.078947e-05
#> 411   6.89942810 0.00001 3.078947e-05
#> 431   7.92742625 0.00001 3.078947e-05
#> 451   7.72917962 0.00001 3.078947e-05
#> 471   2.74487742 0.00307 6.087966e-03
#> 481   6.76081903 0.00001 3.078947e-05
#> 491   5.43992258 0.00001 3.078947e-05
#> 501   4.16982892 0.00003 8.357143e-05
#> 521   6.63637112 0.00001 3.078947e-05
#> 531   5.57828019 0.00001 3.078947e-05
#> 541   4.50103306 0.00001 3.078947e-05
#> 551   6.93828450 0.00001 3.078947e-05
#> 581   4.95783029 0.00001 3.078947e-05
#> 591   4.77163200 0.00002 5.707317e-05
#> 601  -3.63729886 0.99991 1.000000e+00
#> 611  -1.47236505 0.92084 9.884246e-01
#> 621  -0.01968822 0.45741 5.574684e-01
#> 631   7.33168010 0.00001 3.078947e-05
#> 651   3.30047432 0.00076 1.778400e-03
#> 671  -5.64148731 1.00000 1.000000e+00
#> 681  -1.94752137 0.97240 1.000000e+00
#> 691   7.88387686 0.00001 3.078947e-05
#> 701   4.61121886 0.00001 3.078947e-05
#> 711   5.98472047 0.00001 3.078947e-05
#> 721   4.81688097 0.00001 3.078947e-05
#> 741   4.99107383 0.00001 3.078947e-05
#> 751   4.10498503 0.00005 1.329545e-04
#> 771   6.94362725 0.00001 3.078947e-05
#> 781   4.43677857 0.00001 3.078947e-05
#> 791   2.03888494 0.01396 2.401941e-02
#> 801   4.35865757 0.00001 3.078947e-05
#> 811  -0.07351171 0.35192 4.475504e-01
#> 821  -1.17300821 0.86663 9.476235e-01
#> 841  -2.92920344 0.99859 1.000000e+00
#> 861   5.10829172 0.00001 3.078947e-05
#> 871   2.41612955 0.00668 1.202400e-02
#> 88    2.97098725 0.00110 2.475000e-03
#> 89    7.47534690 0.00001 3.078947e-05
#> 90    4.17292288 0.00005 1.329545e-04
#> 93    4.37590969 0.00002 5.707317e-05
#> 96    4.25206470 0.00001 3.078947e-05
#> 97    2.98937580 0.00196 4.095000e-03
#> 98    2.83336256 0.00385 7.265323e-03
#> 99   -0.85178198 0.75318 8.473275e-01
#> 103   1.98162593 0.02103 3.465507e-02
#> 104  -0.03312769 0.43808 5.452698e-01
#> 109   3.11075919 0.00109 2.475000e-03
#> 1111  3.20067782 0.00057 1.418936e-03
#> 112   2.91226712 0.00187 3.978000e-03
#> 115   5.55231137 0.00001 3.078947e-05
#> 116   5.29511262 0.00001 3.078947e-05
#> 117   0.73909389 0.21124 2.952300e-01
#> 
#> $buffered_up
#>             RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 1        ALKBH1            8               14.10532            3.5000243
#> 2        ALKBH5          318              222.87349           10.3492723
#> 3         ATXN2          463              430.58151            6.2625672
#> 4      C17orf85          331              278.89626           10.2577030
#> 5       CAPRIN1          445              349.66914            9.2204675
#> 6         CELF2          340              232.17709           10.3703329
#> 7          CNBP          143              105.13656            8.5987750
#> 8         CPSF1          421              344.89387            9.3235763
#> 9         CPSF2          278              195.86613           10.2257369
#> 10        CPSF3          443              366.05327            8.8176907
#> 11        CPSF4          424              341.93605            9.3896029
#> 12        CPSF6          446              364.73170            8.8491660
#> 13        CPSF7          450              361.63420            8.9350315
#> 14        CSTF2          461              411.39549            7.1989552
#> 15       CSTF2T          463              415.43166            7.0319602
#> 16        DDX3X          463              426.26435            6.4956086
#> 17        DGCR8           96               89.54289            8.0729952
#> 18       DICER1          181              148.31067            9.5800664
#> 19       DIS3L2           45               42.27877            5.8859599
#> 20        EIF3A          191              185.20036           10.1004137
#> 21        EIF3B          255              206.33974           10.2620100
#> 22        EIF3D          232              225.10011           10.3390191
#> 23        EIF3G          338              288.48903           10.1677733
#> 24       EIF4A3          185              182.28914           10.0861646
#> 25       ELAVL1          471              444.26731            5.3966403
#> 26        EWSR1          435              304.06293           10.0081043
#> 27         EZH2            3                1.74410            1.2506682
#> 28          FBL          134              110.21281            8.7340795
#> 29       FIP1L1          463              431.47267            6.1985851
#> 30         FMR1          411              314.05759            9.9096062
#> 31          FUS          463              401.08502            7.6658711
#> 32         FXR1          168               68.66940            7.2745112
#> 33         FXR2          296              181.51268           10.0834139
#> 34      HNRNPA1          449              360.56082            8.9561511
#> 35    HNRNPA2B1           63               38.77045            5.6418166
#> 36       HNRNPC          463              422.86423            6.6946490
#> 37       HNRNPD          418              277.50034           10.2500805
#> 38       HNRNPF          158               84.64139            7.9468698
#> 39       HNRNPH           81               70.02215            7.3409576
#> 40       HNRNPM          152               91.03472            8.1677158
#> 41       HNRNPU          295              177.88711           10.0304602
#> 42      IGF2BP1          420              307.43971            9.9762558
#> 43      IGF2BP2          429              294.01620           10.1256744
#> 44      IGF2BP3          446              335.16443            9.5269019
#> 45       LIN28A          405              322.68610            9.7455742
#> 46       LIN28B          466              431.19708            6.2307967
#> 47        MBNL2          216              131.08386            9.2392305
#> 48      METTL14           42               32.94085            5.2748298
#> 49       METTL3           62               60.95483            6.9181868
#> 50        MOV10          453              366.16674            8.8260057
#> 51        NOP56          160               98.23697            8.3759786
#> 52        NOP58          144              102.91566            8.5041196
#> 53       NUDT21          464              425.85009            6.5404224
#> 54        PRKRA           18               22.40552            4.3877222
#> 55        PTBP1          458              420.24166            6.8165684
#> 56   PTBP1PTBP2          296              228.99031           10.3786460
#> 57         PUM2          259              145.89908            9.5520494
#> 58          QKI          229              140.48793            9.4353116
#> 59        RBM10           72               70.00670            7.3149320
#> 60       RBM15B          443              371.17101            8.6850366
#> 61        RBM15          342              297.97768           10.0901685
#> 62        RBM47          300              163.22249            9.8234268
#> 63        RBPMS          179              102.58011            8.5285468
#> 64         RTCB          392              306.49170            9.9489698
#> 65        SRRM4          313              250.62617           10.3617158
#> 66        SRSF1          425              336.77660            9.4985080
#> 67          SSB          109              128.55295            9.1717445
#> 68        STAU1          266              194.24164           10.2009858
#> 69        TAF15          366              238.84331           10.3860155
#> 70       TARBP2          164              138.98804            9.4192624
#> 71       TARDBP          468              443.87558            5.4169007
#> 72         TIA1          303              212.92048           10.3000436
#> 73        TIAL1          398              300.43488           10.0175368
#> 74       TNRC6A           24               17.18610            3.8549183
#> 75       TNRC6B           19               13.81028            3.4747355
#> 76       TNRC6C           13               11.88234            3.2308018
#> 77        U2AF2          453              404.76978            7.5293709
#> 78         UPF1          444              358.16134            9.0009921
#> 79        WDR33          449              364.29893            8.8563410
#> 80         WTAP           30               33.46585            5.3004353
#> 81       YTHDC1          439              363.99954            8.8951307
#> 82       YTHDC2          119               99.67591            8.3803149
#> 83       YTHDF1          444              385.70228            8.2280982
#> 84       YTHDF2          441              401.88500            7.6403585
#> 85       YTHDF3          262              230.91877           10.3641336
#> 86       ZC3H7B          450              336.98766            9.4656518
#> 87      ZC3H11A          105               97.32813            8.3501256
#> 510       DDX24          112              140.22452            9.4034597
#> 610       XRCC6           48               50.53254            6.3720615
#> 710       SF3B4          229              233.96397           10.3823364
#> 91        NCBP2          172              129.95402            9.2493471
#> 101     FAM120A          229              174.55645            9.9616498
#> 111       LARP7            5               12.57745            3.3132207
#> 121         FTO           33               32.98588            5.2572804
#> 141       RBM22          131              136.43467            9.3435991
#> 151       UCHL5          190              176.95460           10.0191997
#> 161       DDX42           12               18.11761            3.9584435
#> 171        XPO5           65               85.23791            7.9212471
#> 181      TROVE2           17               31.29945            5.1237274
#> 191        SUB1           94               97.12964            8.3017903
#> 201       TBRG4          141              141.41669            9.4484492
#> 241      METAP2           26               50.64882            6.3592717
#> 251        DKC1           37               57.76200            6.7533852
#> 271     KHDRBS1           28               19.35650            4.0722713
#> 281        SLTM           38               59.61607            6.8696623
#> 291        SND1           57               96.62843            8.3214108
#> 311      SERBP1            9               22.09365            4.3403396
#> 321        NONO           61               61.59444            6.9527566
#> 341       NSUN2            1                4.80463            2.0658033
#> 361        PPIG          177              212.56427           10.3022585
#> 381      AKAP8L           27               44.91586            6.0559563
#> 391      HNRNPL          174              124.83559            9.0889606
#> 401       PCBP2          209              192.80094           10.1969885
#> 411       SF3A3           88              143.58775            9.4924276
#> 431      RBFOX2          112              136.09384            9.3472203
#> 451      HNRNPK           66               99.96962            8.4095413
#> 471      GEMIN5          141              134.19041            9.3317139
#> 481       NOL12            8               34.69999            5.3865178
#> 491     SUPV3L1           14               25.59212            4.6419769
#> 501       SRSF7           31               34.79332            5.3761224
#> 521      GTF2F1          109              119.40971            8.9340755
#> 531      EFTUD2          228              242.30207           10.3423988
#> 541       PRPF8          342              329.81939            9.5792029
#> 551       DDX55           94              129.05584            9.1773615
#> 581       BUD13          248              263.23318           10.3101615
#> 591       YWHAG            1                5.11895            2.1333932
#> 601       FUBP3          125               75.91748            7.5932027
#> 611       MATR3          135               97.95541            8.3479470
#> 621       TRA2A           43               47.81325            6.2251716
#> 631      DROSHA           95              107.43889            8.6386047
#> 651       SRSF9           20               27.82805            4.8448280
#> 671       KHSRP          278              162.10403            9.8045117
#> 681      EXOSC5           86               52.53892            6.5035873
#> 691       LARP4           75              117.08972            8.8914476
#> 701       CDC40           17               36.03708            5.4784784
#> 711       LSM11           56               76.88667            7.6001121
#> 721     FASTKD2           58               71.52327            7.3976176
#> 741      SMNDC1           33               63.67961            7.0290906
#> 751       DHX30           28               39.88496            5.7344707
#> 771        RPS3          100              154.73475            9.7013686
#> 781        XRN2          129              135.14173            9.3276092
#> 791        SBDS            0                0.91884            0.9083949
#> 801        RBM5           13               41.01635            5.8142464
#> 811        NPM1            4                3.58836            1.7917881
#> 821      ZNF622          202              159.68602            9.7797164
#> 841       SUGP2          134               81.38204            7.7751444
#> 861       NOLC1           37               60.07557            6.8506385
#> 871        YBX3          195              178.34693           10.0268696
#> 88        SF3B1          306              258.04346           10.3428575
#> 89        DDX59           32               52.57565            6.4595800
#> 90        GRSF1           20               29.87557            4.9991956
#> 93         NKRF           63               72.32616            7.4260562
#> 96       ZRANB2           45               55.80061            6.6461504
#> 97        MTPAP           19               29.45772            5.0034411
#> 98         PUS1            6                2.65958            1.5411937
#> 99        SAFB2           29               21.90376            4.3373247
#> 103       BCCIP            9                9.41321            2.8846127
#> 104    HNRNPUL1            6               20.54104            4.1929598
#> 109        AARS           53               59.83428            6.8784744
#> 1111       DDX6          114              132.04228            9.2266648
#> 112      EIF4G2           60               62.82293            7.0094656
#> 115       EIF3H           48               73.79983            7.4843077
#> 116       FKBP4            7               21.99365            4.3347915
#> 117       GRWD1          241              209.07895           10.2731394
#>                 z    pval         padj
#> 1    -1.744365025 0.95358 9.647986e-01
#> 2     9.191613441 0.00001 1.365079e-05
#> 3     5.176549617 0.00001 1.365079e-05
#> 4     5.079474435 0.00001 1.365079e-05
#> 5    10.339048371 0.00001 1.365079e-05
#> 6    10.397246741 0.00001 1.365079e-05
#> 7     4.403352825 0.00002 2.687500e-05
#> 8     8.162761504 0.00001 1.365079e-05
#> 9     8.032073489 0.00001 1.365079e-05
#> 10    8.726403811 0.00001 1.365079e-05
#> 11    8.739874427 0.00001 1.365079e-05
#> 12    9.183724248 0.00001 1.365079e-05
#> 13    9.889813993 0.00001 1.365079e-05
#> 14    6.890515172 0.00001 1.365079e-05
#> 15    6.764591765 0.00001 1.365079e-05
#> 16    5.655459295 0.00001 1.365079e-05
#> 17    0.799840679 0.19359 2.219832e-01
#> 18    3.412223745 0.00036 4.763077e-04
#> 19    0.462325611 0.29002 3.197656e-01
#> 20    0.574198265 0.26568 2.967335e-01
#> 21    4.741786440 0.00001 1.365079e-05
#> 22    0.667364080 0.23728 2.685011e-01
#> 23    4.869401454 0.00001 1.365079e-05
#> 24    0.268770154 0.37355 3.966086e-01
#> 25    4.953580094 0.00001 1.365079e-05
#> 26   13.083104076 0.00001 1.365079e-05
#> 27    1.004183241 0.08890 1.033162e-01
#> 28    2.723491353 0.00323 4.145970e-03
#> 29    5.086213945 0.00001 1.365079e-05
#> 30    9.782670256 0.00001 1.365079e-05
#> 31    8.076705105 0.00001 1.365079e-05
#> 32   13.654608209 0.00001 1.365079e-05
#> 33   11.354023655 0.00001 1.365079e-05
#> 34    9.874685981 0.00001 1.365079e-05
#> 35    4.294636198 0.00001 1.365079e-05
#> 36    5.995201577 0.00001 1.365079e-05
#> 37   13.707176233 0.00001 1.365079e-05
#> 38    9.231132684 0.00001 1.365079e-05
#> 39    1.495424796 0.06208 7.313534e-02
#> 40    7.464177425 0.00001 1.365079e-05
#> 41   11.675724502 0.00001 1.365079e-05
#> 42   11.282819118 0.00001 1.365079e-05
#> 43   13.330845348 0.00001 1.365079e-05
#> 44   11.633957273 0.00001 1.365079e-05
#> 45    8.446285286 0.00001 1.365079e-05
#> 46    5.585629189 0.00001 1.365079e-05
#> 47    9.190823798 0.00001 1.365079e-05
#> 48    1.717429818 0.03821 4.628254e-02
#> 49    0.151075713 0.40333 4.230046e-01
#> 50    9.838341700 0.00001 1.365079e-05
#> 51    7.373828572 0.00001 1.365079e-05
#> 52    4.831110331 0.00001 1.365079e-05
#> 53    5.832942866 0.00001 1.365079e-05
#> 54   -1.004056275 0.81094 8.302481e-01
#> 55    5.539200613 0.00001 1.365079e-05
#> 56    6.456496334 0.00001 1.365079e-05
#> 57   11.840487352 0.00001 1.365079e-05
#> 58    9.380937691 0.00001 1.365079e-05
#> 59    0.272497407 0.36264 3.898380e-01
#> 60    8.270430364 0.00001 1.365079e-05
#> 61    4.362892462 0.00001 1.365079e-05
#> 62   13.923604606 0.00001 1.365079e-05
#> 63    8.960482037 0.00001 1.365079e-05
#> 64    8.594688832 0.00001 1.365079e-05
#> 65    6.019642989 0.00001 1.365079e-05
#> 66    9.288132388 0.00001 1.365079e-05
#> 67   -2.131868148 0.98231 9.823100e-01
#> 68    7.034453494 0.00001 1.365079e-05
#> 69   12.243067672 0.00001 1.365079e-05
#> 70    2.655405374 0.00384 4.856471e-03
#> 71    4.453546664 0.00001 1.365079e-05
#> 72    8.745547474 0.00001 1.365079e-05
#> 73    9.739432184 0.00001 1.365079e-05
#> 74    1.767586082 0.03355 4.121857e-02
#> 75    1.493558270 0.05659 6.759361e-02
#> 76    0.345938890 0.29940 3.259291e-01
#> 77    6.405610887 0.00001 1.365079e-05
#> 78    9.536577593 0.00001 1.365079e-05
#> 79    9.563889916 0.00001 1.365079e-05
#> 80   -0.653880254 0.70597 7.314870e-01
#> 81    8.431630982 0.00001 1.365079e-05
#> 82    2.305890675 0.00957 1.192783e-02
#> 83    7.085199843 0.00001 1.365079e-05
#> 84    5.119524136 0.00001 1.365079e-05
#> 85    2.998922160 0.00112 1.459394e-03
#> 86   11.939203192 0.00001 1.365079e-05
#> 87    0.918773004 0.16438 5.197962e-01
#> 510  -3.001503786 0.99860 1.000000e+00
#> 610  -0.397444378 0.62005 1.000000e+00
#> 710  -0.478116851 0.66704 1.000000e+00
#> 91    4.545832225 0.00001 6.882353e-05
#> 101   5.465314618 0.00001 6.882353e-05
#> 111  -2.287034500 0.98958 1.000000e+00
#> 121   0.002685799 0.45443 1.000000e+00
#> 141  -0.581646314 0.70003 1.000000e+00
#> 151   1.302040122 0.08943 3.264300e-01
#> 161  -1.545458479 0.92709 1.000000e+00
#> 171  -2.554889355 0.99473 1.000000e+00
#> 181  -2.790829595 0.99783 1.000000e+00
#> 191  -0.376983744 0.62137 1.000000e+00
#> 201  -0.044101417 0.49323 1.000000e+00
#> 241  -3.876044463 0.99995 1.000000e+00
#> 251  -3.074310037 0.99917 1.000000e+00
#> 271   2.122525577 0.01564 7.624500e-02
#> 281  -3.146598647 0.99935 1.000000e+00
#> 291  -4.762224917 1.00000 1.000000e+00
#> 311  -3.016733973 0.99916 1.000000e+00
#> 321  -0.085497024 0.50031 1.000000e+00
#> 341  -1.841719401 0.95969 1.000000e+00
#> 361  -3.452084794 0.99962 1.000000e+00
#> 381  -2.958386585 0.99894 1.000000e+00
#> 391   5.409244464 0.00001 6.882353e-05
#> 401   1.588612164 0.05133 2.070900e-01
#> 411  -5.856009924 1.00000 1.000000e+00
#> 431  -2.577647597 0.99487 1.000000e+00
#> 451  -4.039414153 0.99999 1.000000e+00
#> 471   0.729725542 0.21613 6.465712e-01
#> 481  -4.956818288 1.00000 1.000000e+00
#> 491  -2.497237770 0.99425 1.000000e+00
#> 501  -0.705586610 0.72599 1.000000e+00
#> 521  -1.165169238 0.86673 1.000000e+00
#> 531  -1.382858102 0.90844 1.000000e+00
#> 541   1.271568226 0.09207 3.264300e-01
#> 551  -3.819816839 0.99996 1.000000e+00
#> 581  -1.477491890 0.92322 1.000000e+00
#> 591  -1.930703661 0.97105 1.000000e+00
#> 601   6.464007618 0.00001 6.882353e-05
#> 611   4.437568918 0.00001 6.882353e-05
#> 621  -0.773191532 0.75250 1.000000e+00
#> 631  -1.439918873 0.91943 1.000000e+00
#> 651  -1.615753969 0.93972 1.000000e+00
#> 671  11.820677452 0.00001 6.882353e-05
#> 681   5.145018939 0.00001 6.882353e-05
#> 691  -4.733730858 1.00000 1.000000e+00
#> 701  -3.474884581 0.99991 1.000000e+00
#> 711  -2.748205517 0.99713 1.000000e+00
#> 721  -1.828057461 0.96312 1.000000e+00
#> 741  -4.364662778 0.99999 1.000000e+00
#> 751  -2.072546978 0.97976 1.000000e+00
#> 771  -5.641961665 1.00000 1.000000e+00
#> 781  -0.658446328 0.72678 1.000000e+00
#> 791  -1.011498406 0.61986 1.000000e+00
#> 801  -4.818569470 1.00000 1.000000e+00
#> 811   0.229736987 0.28829 8.226812e-01
#> 821   4.326708285 0.00001 6.882353e-05
#> 841   6.767457619 0.00001 6.882353e-05
#> 861  -3.368382355 0.99965 1.000000e+00
#> 871   1.660844369 0.04372 1.857375e-01
#> 88    4.636681890 0.00001 6.882353e-05
#> 89   -3.185292211 0.99947 1.000000e+00
#> 90   -1.975431790 0.97422 1.000000e+00
#> 93   -1.255869834 0.88295 1.000000e+00
#> 96   -1.625092621 0.94205 1.000000e+00
#> 97   -2.090105552 0.98165 1.000000e+00
#> 98    2.167423881 0.01332 6.969130e-02
#> 99    1.636086853 0.04445 1.857375e-01
#> 103  -0.143246264 0.47026 1.000000e+00
#> 104  -3.467965552 0.99992 1.000000e+00
#> 109  -0.993574971 0.82110 1.000000e+00
#> 1111 -1.955449813 0.97246 1.000000e+00
#> 112  -0.402731132 0.62486 1.000000e+00
#> 115  -3.447189900 0.99975 1.000000e+00
#> 116  -3.458909133 0.99992 1.000000e+00
#> 117   3.107234204 0.00070 4.310526e-03
#> 
#> $buffered_down
#>             RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 1        ALKBH1           26               16.76099            3.7878817
#> 2        ALKBH5          238              264.29409           11.1612019
#> 3         ATXN2          524              510.60061            6.7350064
#> 4      C17orf85          376              330.71134           11.0343297
#> 5       CAPRIN1          384              414.61692            9.9003247
#> 6         CELF2          230              275.37189           11.1370377
#> 7          CNBP           87              124.66612            9.2652978
#> 8         CPSF1          408              409.01211            9.9967466
#> 9         CPSF2          187              232.25944           10.9551555
#> 10        CPSF3          430              434.09006            9.5003427
#> 11        CPSF4          408              405.47406           10.0818264
#> 12        CPSF6           25               27.51397            4.7839433
#> 13        CPSF7          418              428.83500            9.6012175
#> 14        CSTF2          489              487.88709            7.7800904
#> 15       CSTF2T          195              152.32679            9.9097982
#> 16        DDX3X          420              367.01828           10.7021928
#> 17        DGCR8          154              106.17463            8.7268915
#> 18       DICER1          179              175.92309           10.3087302
#> 19       DIS3L2           52               50.07618            6.3536084
#> 20        EIF3A          225              219.66132           10.9141554
#> 21        EIF3B          264              244.61561           11.0707010
#> 22        EIF3D           53               35.04589            5.3728570
#> 23        EIF3G          356              342.07248           10.9603005
#> 24       EIF4A3          288              216.13247           10.7823534
#> 25       ELAVL1          521              526.85308            5.7839579
#> 26        EWSR1          101               51.19384            6.4034972
#> 27         EZH2            2                2.06549            1.3371234
#> 28          FBL          147              130.68168            9.3886013
#> 29       FIP1L1          530              511.66931            6.6827779
#> 30         FMR1          155              114.82521            8.9819566
#> 31          FUS          458              475.59383            8.2200117
#> 32         FXR1           42               23.03509            4.4174239
#> 33         FXR2          150               92.07904            8.2419507
#> 34      HNRNPA1           10                7.51459            2.5528753
#> 35    HNRNPA2B1           33               46.00606            6.0894593
#> 36       HNRNPC          517              501.47555            7.2031938
#> 37       HNRNPD          196              329.08036           11.0345965
#> 38       HNRNPF           39              100.33219            8.5278477
#> 39       HNRNPH           86               83.02043            7.8481659
#> 40       HNRNPM          166              162.01694           10.0724479
#> 41       HNRNPU           26               19.06871            4.0198645
#> 42      IGF2BP1          305              364.49864           10.7192820
#> 43      IGF2BP2          123              104.25928            8.6396771
#> 44      IGF2BP3           58               46.25077            6.0928151
#> 45       LIN28A          371              382.63603           10.4784757
#> 46       LIN28B          225              197.81056           10.6402521
#> 47        MBNL2           85              155.47089            9.9572026
#> 48      METTL14           42               39.07614            5.6580282
#> 49       METTL3           95               72.30757            7.4486714
#> 50        MOV10          360              434.15968            9.4797448
#> 51        NOP56           74              116.48943            9.0339872
#> 52        NOP58          112              122.07207            9.1810217
#> 53       NUDT21          521              504.98768            6.9842864
#> 54        PRKRA           22               26.56319            4.7089031
#> 55        PTBP1          511              498.30207            7.3151363
#> 56   PTBP1PTBP2          255              271.57739           11.1758230
#> 57         PUM2           80              173.07125           10.3090230
#> 58          QKI           75              166.61606           10.1898238
#> 59        RBM10           62               83.00724            7.8787745
#> 60       RBM15B          479              440.17636            9.3134067
#> 61        RBM15          252              194.52098           10.5678870
#> 62        RBM47           74              193.55307           10.5802147
#> 63        RBPMS           90              121.57266            9.1574843
#> 64         RTCB          365              363.46662           10.7188943
#> 65        SRRM4          346              297.22531           11.0946517
#> 66        SRSF1          156              106.92610            8.7391145
#> 67          SSB          176              152.36540            9.9258737
#> 68        STAU1          246              230.29610           10.9715463
#> 69        TAF15           37               31.82890            5.1332980
#> 70       TARBP2          165              164.84291           10.1592620
#> 71       TARDBP          536              526.36013            5.8172515
#> 72         TIA1          166              171.72040           10.2526189
#> 73        TIAL1          321              356.26016           10.8207000
#> 74       TNRC6A           14               20.39064            4.1472488
#> 75       TNRC6B           14               16.40103            3.7472637
#> 76       TNRC6C            9               14.09449            3.4750140
#> 77        U2AF2          498              479.97742            8.0804878
#> 78         UPF1          165              158.77730           10.0392038
#> 79        WDR33          363              432.00526            9.5147989
#> 80         WTAP           46               39.68468            5.6977415
#> 81       YTHDC1          469              431.68272            9.5468427
#> 82       YTHDC2          139              118.19189            9.0738741
#> 83       YTHDF1          485              457.40982            8.8081737
#> 84       YTHDF2          503              476.59744            8.1912622
#> 85       YTHDF3          304              273.80615           11.1696329
#> 86       ZC3H7B          304              399.60191           10.1934358
#> 87      ZC3H11A          138              115.43217            9.0168682
#> 510       DDX24          286              166.22980           10.2061292
#> 610       XRCC6           73               59.94602            6.8637612
#> 710       SF3B4          374              277.38540           11.2013420
#> 91        NCBP2          158              154.06702            9.9120538
#> 101     FAM120A          220              207.09524           10.7716921
#> 111       LARP7           27               14.90794            3.5783366
#> 121         FTO           56               39.06991            5.6394397
#> 141       RBM22          205              161.74891           10.0966976
#> 151       UCHL5          268              209.87087           10.8111648
#> 161       DDX42           31               21.45576            4.2641839
#> 171        XPO5          120              101.10699            8.5067836
#> 181      TROVE2           58               37.12277            5.5307778
#> 191        SUB1          149              115.14156            8.9906590
#> 201       TBRG4          194              167.73783           10.1560272
#> 241      METAP2          121               60.07804            6.9001526
#> 251        DKC1           87               68.56953            7.2687202
#> 271     KHDRBS1           13               22.91891            4.3973910
#> 281        SLTM          108               70.74760            7.3700283
#> 291        SND1          190              114.54747            8.9873853
#> 311      SERBP1           51               26.19087            4.6964241
#> 321        NONO          113               73.03954            7.4859787
#> 341       NSUN2           11                5.72277            2.2244332
#> 361        PPIG          374              252.09711           11.1276933
#> 381      AKAP8L          112               53.26130            6.5184498
#> 391      HNRNPL           94              148.05911            9.8221678
#> 401       PCBP2          285              228.63631           10.9429803
#> 411       SF3A3          302              170.23357           10.2604019
#> 431      RBFOX2          239              161.43552           10.0632070
#> 451      HNRNPK          176              118.58712            9.0731204
#> 471      GEMIN5          204              159.10338           10.0386304
#> 481       NOL12          113               41.13218            5.7899520
#> 491     SUPV3L1           52               30.33546            5.0241954
#> 501       SRSF7           49               41.28290            5.7806766
#> 521      GTF2F1          202              141.69659            9.7019490
#> 531      EFTUD2          400              287.35784           11.1698755
#> 541       PRPF8          465              391.11108           10.3410498
#> 551       DDX55          263              153.06908            9.9313752
#> 581       BUD13          419              312.15691           11.1509162
#> 591       YWHAG           13                6.05141            2.3016733
#> 601       FUBP3           50               90.07148            8.1621319
#> 611       MATR3           76              116.14189            8.9932700
#> 621       TRA2A           67               56.65426            6.6772786
#> 631      DROSHA          182              127.40097            9.3188755
#> 651       SRSF9           45               33.04037            5.2363665
#> 671       KHSRP          110              192.18659           10.5733350
#> 681      EXOSC5           44               62.22377            7.0134520
#> 691       LARP4          242              138.82322            9.6280702
#> 701       CDC40           68               42.72346            5.8713704
#> 711       LSM11          153               91.13664            8.2036932
#> 721     FASTKD2          107               84.79094            7.9797062
#> 741      SMNDC1          155               75.47400            7.6256413
#> 751       DHX30           80               47.30909            6.1725177
#> 771        RPS3          284              183.50911           10.4302054
#> 781        XRN2          223              160.32658           10.0437181
#> 791        SBDS            1                1.08875            0.9786741
#> 801        RBM5          109               48.66012            6.2550806
#> 811        NPM1            1                4.24990            1.9201268
#> 821      ZNF622          219              189.27374           10.5315989
#> 841       SUGP2           47               96.50396            8.4182441
#> 861       NOLC1          113               71.23063            7.4328146
#> 871        YBX3          244              211.52658           10.8447715
#> 88        SF3B1          349              305.92142           11.1658977
#> 89        DDX59          101               62.38326            6.9634989
#> 90        GRSF1           74               35.42806            5.4204408
#> 93         NKRF          101               85.79544            7.9459856
#> 96       ZRANB2           91               66.09523            7.1766159
#> 97        MTPAP           59               34.93228            5.3855756
#> 98         PUS1            2                3.15278            1.6592908
#> 99        SAFB2           20               25.94594            4.6657063
#> 103       BCCIP           11               11.13920            3.1026472
#> 104    HNRNPUL1           38               24.39263            4.5261857
#> 109        AARS          112               70.98194            7.3782273
#> 1111       DDX6          219              156.54738           10.0054063
#> 112      EIF4G2          104               74.50297            7.5327776
#> 115       EIF3H          132               87.50674            8.0704898
#> 116       FKBP4           50               26.08355            4.6773035
#> 117       GRWD1          272              247.90130           11.1009860
#>                 z    pval         padj
#> 1      2.43909677 0.00700 3.728353e-02
#> 2     -2.35584755 0.98970 1.000000e+00
#> 3      1.98951407 0.01686 6.590727e-02
#> 4      4.10434173 0.00002 2.457143e-04
#> 5     -3.09251674 0.99851 1.000000e+00
#> 6     -4.07396395 0.99997 1.000000e+00
#> 7     -4.06528975 0.99999 1.000000e+00
#> 8     -0.10124394 0.52261 1.000000e+00
#> 9     -4.13133708 0.99995 1.000000e+00
#> 10    -0.43051710 0.65038 1.000000e+00
#> 11     0.25054389 0.38533 8.720626e-01
#> 12    -0.52550163 0.65508 7.901480e-01
#> 13    -1.12850272 0.85937 1.000000e+00
#> 14     0.14304590 0.42135 9.139435e-01
#> 15     4.30616337 0.00001 2.600000e-05
#> 16     4.95054807 0.00001 2.600000e-05
#> 17     5.48022968 0.00001 1.433333e-04
#> 18     0.29847614 0.36028 8.374076e-01
#> 19     0.30279172 0.34688 8.286578e-01
#> 20     0.48915192 0.29585 7.710030e-01
#> 21     1.75096320 0.03690 1.217696e-01
#> 22     3.34163181 0.00058 1.233818e-03
#> 23     1.27072428 0.09404 2.788772e-01
#> 24     6.66529166 0.00001 1.433333e-04
#> 25    -1.01195066 0.82294 1.000000e+00
#> 26     7.77796236 0.00001 2.600000e-05
#> 27    -0.04897828 0.34185 8.286578e-01
#> 28     1.73809916 0.03823 1.217696e-01
#> 29     2.74297459 0.00153 1.096500e-02
#> 30     4.47283277 0.00001 2.600000e-05
#> 31    -2.14036558 0.97958 1.000000e+00
#> 32     4.29320583 0.00001 2.600000e-05
#> 33     7.02757901 0.00001 2.600000e-05
#> 34     0.97357282 0.12417 1.793567e-01
#> 35    -2.13583168 0.98345 1.000000e+00
#> 36     2.15521759 0.01080 4.422857e-02
#> 37   -12.06028331 1.00000 1.000000e+00
#> 38    -7.19198934 1.00000 1.000000e+00
#> 39     0.37965176 0.32640 8.256000e-01
#> 40     0.39544111 0.32667 4.246710e-01
#> 41     1.72425961 0.03652 5.622158e-02
#> 42    -5.55061802 1.00000 1.000000e+00
#> 43     2.16914588 0.01408 2.288000e-02
#> 44     1.92837462 0.02466 3.952356e-02
#> 45    -1.11046972 0.85702 1.000000e+00
#> 46     2.55533795 0.00483 8.829844e-03
#> 47    -7.07737832 1.00000 1.000000e+00
#> 48     0.51676307 0.26770 7.194437e-01
#> 49     3.04650706 0.00120 9.381818e-03
#> 50    -7.82296167 1.00000 1.000000e+00
#> 51    -4.70328651 1.00000 1.000000e+00
#> 52    -1.09705329 0.85206 1.000000e+00
#> 53     2.29262077 0.00737 3.728353e-02
#> 54    -0.96905582 0.80447 1.000000e+00
#> 55     1.73584326 0.03261 1.168525e-01
#> 56    -1.48332610 0.92429 1.000000e+00
#> 57    -9.02813483 1.00000 1.000000e+00
#> 58    -8.99093665 1.00000 1.000000e+00
#> 59    -2.66630805 0.99622 1.000000e+00
#> 60     4.16857560 0.00001 1.433333e-04
#> 61     5.43902672 0.00001 2.600000e-05
#> 62   -11.29968280 1.00000 1.000000e+00
#> 63    -3.44774382 0.99979 1.000000e+00
#> 64     0.14305393 0.42509 9.139435e-01
#> 65     4.39623445 0.00001 1.433333e-04
#> 66     5.61543162 0.00001 2.600000e-05
#> 67     2.38111030 0.00797 3.807889e-02
#> 68     1.43132969 0.07022 2.156757e-01
#> 69     1.00736408 0.13525 1.929787e-01
#> 70     0.01546274 0.47247 9.910346e-01
#> 71     1.65711763 0.03697 1.217696e-01
#> 72    -0.55794525 0.69458 8.208673e-01
#> 73    -3.25858402 0.99938 1.000000e+00
#> 74    -1.54093481 0.92757 1.000000e+00
#> 75    -0.64074220 0.68343 1.000000e+00
#> 76    -1.46603438 0.91293 1.000000e+00
#> 77     2.23038268 0.00935 4.232105e-02
#> 78     0.61983999 0.24982 3.521559e-01
#> 79    -7.25241394 1.00000 1.000000e+00
#> 80     1.10839006 0.11619 3.330780e-01
#> 81     3.90886090 0.00001 1.433333e-04
#> 82     2.29318918 0.01048 4.422857e-02
#> 83     3.13233833 0.00047 4.042000e-03
#> 84     3.22325905 0.00034 3.248889e-03
#> 85     2.70320880 0.00294 1.944923e-02
#> 86    -9.37877197 1.00000 1.000000e+00
#> 87     2.50284572 0.00635 1.125682e-02
#> 510   11.73512479 0.00001 2.600000e-05
#> 610    1.90186979 0.02663 4.154280e-02
#> 710    8.62527012 0.00001 2.600000e-05
#> 91     0.39678760 0.32473 4.246710e-01
#> 101    1.19802534 0.10711 1.566484e-01
#> 111    3.37924049 0.00056 1.213333e-03
#> 121    3.00208722 0.00150 3.025862e-03
#> 141    4.28368676 0.00003 7.468085e-05
#> 151    5.37676847 0.00001 2.600000e-05
#> 161    2.23823365 0.01177 1.995783e-02
#> 171    2.22093460 0.01220 2.039143e-02
#> 181    3.77473669 0.00018 4.050000e-04
#> 191    3.76595754 0.00011 2.523529e-04
#> 201    2.58587039 0.00460 8.542857e-03
#> 241    8.82907430 0.00001 2.600000e-05
#> 251    2.53558665 0.00557 1.002600e-02
#> 271   -2.25563521 0.98850 1.000000e+00
#> 281    5.05458031 0.00001 2.600000e-05
#> 291    8.39538170 0.00001 2.600000e-05
#> 311    5.28255744 0.00001 2.600000e-05
#> 321    5.33804084 0.00001 2.600000e-05
#> 341    2.37239312 0.00848 1.480836e-02
#> 361   10.95491104 0.00001 2.600000e-05
#> 381    9.01114563 0.00001 2.600000e-05
#> 391   -5.50378602 1.00000 1.000000e+00
#> 401    5.15067090 0.00001 2.600000e-05
#> 411   12.84222888 0.00001 2.600000e-05
#> 431    7.70772976 0.00001 2.600000e-05
#> 451    6.32779878 0.00001 2.600000e-05
#> 471    4.47238500 0.00001 2.600000e-05
#> 481   12.41250714 0.00001 2.600000e-05
#> 491    4.31204172 0.00002 5.086957e-05
#> 501    1.33498214 0.08042 1.206300e-01
#> 521    6.21559753 0.00001 2.600000e-05
#> 531   10.08445978 0.00001 2.600000e-05
#> 541    7.14520494 0.00001 2.600000e-05
#> 551   11.06905310 0.00001 2.600000e-05
#> 581    9.58155257 0.00001 2.600000e-05
#> 591    3.01892975 0.00163 3.178500e-03
#> 601   -4.90943795 1.00000 1.000000e+00
#> 611   -4.46354772 1.00000 1.000000e+00
#> 621    1.54939469 0.05476 8.320675e-02
#> 631    5.85897190 0.00001 2.600000e-05
#> 651    2.28395586 0.01083 1.863397e-02
#> 671   -7.77300543 1.00000 1.000000e+00
#> 681   -2.59840233 0.99511 1.000000e+00
#> 691   10.71624716 0.00001 2.600000e-05
#> 701    4.30504947 0.00004 9.750000e-05
#> 711    7.54091582 0.00001 2.600000e-05
#> 721    2.78319270 0.00283 5.340484e-03
#> 741   10.42876221 0.00001 2.600000e-05
#> 751    5.29620355 0.00001 2.600000e-05
#> 771    9.63460316 0.00001 2.600000e-05
#> 781    6.24006161 0.00001 2.600000e-05
#> 791   -0.09068392 0.29950 4.074593e-01
#> 801    9.64653922 0.00001 2.600000e-05
#> 811   -1.69254449 0.93817 1.000000e+00
#> 821    2.82257806 0.00229 4.392295e-03
#> 841   -5.88055652 1.00000 1.000000e+00
#> 861    5.61958995 0.00001 2.600000e-05
#> 871    2.99438489 0.00126 2.586316e-03
#> 88     3.85804895 0.00006 1.432653e-04
#> 89     5.54559434 0.00001 2.600000e-05
#> 90     7.11601536 0.00001 2.600000e-05
#> 93     1.91348951 0.02571 4.064959e-02
#> 96     3.47026652 0.00030 6.622642e-04
#> 97     4.46892254 0.00001 2.600000e-05
#> 98    -0.69474262 0.62676 7.638637e-01
#> 99    -1.27439226 0.88068 1.000000e+00
#> 103   -0.04486491 0.43849 5.492596e-01
#> 104    3.00636586 0.00162 3.178500e-03
#> 109    5.55933807 0.00001 2.600000e-05
#> 1111   6.24188742 0.00001 2.600000e-05
#> 112    3.91582382 0.00011 2.523529e-04
#> 115    5.51308049 0.00001 2.600000e-05
#> 116    5.11329870 0.00001 2.600000e-05
#> 117    2.17086121 0.01403 2.288000e-02
```

Extract the RBP LFC from the RIBO\_LFC and keep only detected RBPs in
res

``` r
# Change of RBPs
rbp_lfc=rbp_change(res=res,ribo_lfc=ribo_lfc)
head(rbp_lfc)
#>                  geneID IDENTIFIER  FoldChange
#> FKBP4   ENSG00000004478      FKBP4 -0.14794320
#> AKAP8L  ENSG00000011243     AKAP8L -0.25094281
#> PTBP1   ENSG00000011304      PTBP1 -0.06274011
#> CELF2   ENSG00000048740      CELF2 -0.27227128
#> FAM120A ENSG00000048828    FAM120A  0.07316692
#> PUM2    ENSG00000055917       PUM2 -0.09197738
# Cure res data by removing RBPs that are not in the rbp_lfc dataframe
res=cure_res(res=res,rbp_lfc=rbp_lfc)
head(res)
#> $forwarded_up
#>         RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 6     CELF2          597              673.85913            15.640422
#> 8     CPSF1          940             1000.85109            14.017409
#> 11    CPSF4          928              992.26129            14.131421
#> 15   CSTF2T         1173             1205.54419            10.551160
#> 22    EIF3D          621              653.39952            15.604828
#> 26    EWSR1          812              882.33524            15.011065
#> 31      FUS         1117             1163.98890            11.480535
#> 37   HNRNPD          725              805.37998            15.422630
#> 38   HNRNPF          183              245.58702            11.929098
#> 41   HNRNPU          423              516.42804            15.120107
#> 42  IGF2BP1          804              892.14414            14.939946
#> 43  IGF2BP2          765              853.23669            15.125983
#> 47    MBNL2          277              380.43798            13.928886
#> 49   METTL3          161              177.04546            10.406970
#> 50    MOV10         1005             1062.63564            13.268308
#> 55    PTBP1         1193             1219.49851            10.205809
#> 57     PUM2          334              423.47940            14.381784
#> 58      QKI          317              407.66707            14.191437
#> 62    RBM47          360              473.67510            14.823691
#> 63    RBPMS          245              297.53573            12.841437
#> 71   TARDBP         1266             1288.11378             8.148649
#> 72     TIA1          541              618.04386            15.496954
#> 74   TNRC6A           43               49.89956             5.806728
#> 81   YTHDC1          994             1056.42888            13.347192
#> 83   YTHDF1         1068             1119.38490            12.362137
#> 86   ZC3H7B          926              977.99133            14.233938
#> 510   DDX24          370              406.79695            14.246241
#> 91    NCBP2          350              377.18147            13.872721
#> 101 FAM120A          458              506.71068            15.046409
#> 141   RBM22          365              395.88323            14.130708
#> 191    SUB1          264              281.88869            12.579825
#> 201   TBRG4          360              410.40015            14.259464
#> 241  METAP2          121              147.00859             9.605555
#> 281    SLTM          176              173.16894            10.311872
#> 381  AKAP8L          123              130.33604             9.118056
#> 401   PCBP2          505              559.58582            15.337322
#> 491 SUPV3L1           68               74.23677             7.043703
#> 501   SRSF7          105              100.97112             8.122575
#> 521  GTF2F1          312              346.60956            13.525562
#> 621   TRA2A          136              138.70236             9.419076
#> 671   KHSRP          376              470.45951            14.835055
#> 691   LARP4          319              339.80580            13.416817
#> 711   LSM11          222              223.02849            11.464125
#> 821  ZNF622          422              463.27868            14.775255
#> 871    YBX3          481              517.67104            15.100574
#> 96   ZRANB2          145              161.87870             9.995831
#> 99    SAFB2           53               63.55648             6.501998
#> 109    AARS          168              173.73775            10.283820
#> 112  EIF4G2          163              182.35542            10.579130
#> 116   FKBP4           64               63.84357             6.537739
#>               z    pval padj
#> 6   -4.91413391 1.00000    1
#> 8   -4.34110820 0.99999    1
#> 11  -4.54740474 1.00000    1
#> 15  -3.08441816 0.99867    1
#> 22  -2.07624973 0.97964    1
#> 26  -4.68555967 1.00000    1
#> 31  -4.09291916 0.99995    1
#> 37  -5.21182058 1.00000    1
#> 38  -5.24658423 1.00000    1
#> 41  -6.17905928 1.00000    1
#> 42  -5.89989671 1.00000    1
#> 43  -5.83345144 1.00000    1
#> 47  -7.42614902 1.00000    1
#> 49  -1.54179945 0.93325    1
#> 50  -4.34385750 1.00000    1
#> 55  -2.59641443 0.99439    1
#> 57  -6.22171780 1.00000    1
#> 58  -6.38885763 1.00000    1
#> 62  -7.66847470 1.00000    1
#> 63  -4.09110980 0.99998    1
#> 71  -2.71379715 0.99569    1
#> 72  -4.97154866 1.00000    1
#> 74  -1.18820102 0.86413    1
#> 81  -4.67730455 0.99999    1
#> 83  -4.15663585 0.99998    1
#> 86  -3.65263139 0.99985    1
#> 510 -2.58292353 0.99474    1
#> 91  -1.95934672 0.97290    1
#> 101 -3.23736253 0.99948    1
#> 141 -2.18554016 0.98467    1
#> 191 -1.42201419 0.91634    1
#> 201 -3.53450514 0.99990    1
#> 241 -2.70766147 0.99631    1
#> 281  0.27454375 0.37087    1
#> 381 -0.80456183 0.77266    1
#> 401 -3.55901896 0.99976    1
#> 491 -0.88543915 0.79136    1
#> 501  0.49601019 0.28625    1
#> 521 -2.55882603 0.99430    1
#> 621 -0.28690287 0.59043    1
#> 671 -6.36731775 1.00000    1
#> 691 -1.55072545 0.93616    1
#> 711 -0.08971378 0.51789    1
#> 821 -2.79377107 0.99740    1
#> 871 -2.42845343 0.99142    1
#> 96  -1.68857390 0.95110    1
#> 99  -1.62357488 0.94105    1
#> 109 -0.55793956 0.69543    1
#> 112 -1.82958520 0.96375    1
#> 116  0.02392723 0.45649    1
#> 
#> $forwarded_down
#>         RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 6     CELF2          781              775.22340            16.187287
#> 8     CPSF1         1108             1151.28655            14.506101
#> 11    CPSF4         1096             1141.50593            14.620235
#> 15   CSTF2T         1349             1386.76599            10.978578
#> 22    EIF3D          708              751.66512            16.191518
#> 26    EWSR1         1045             1014.99696            15.592938
#> 31      FUS         1328             1338.91918            11.958578
#> 37   HNRNPD         1015              926.54498            15.979138
#> 38   HNRNPF          342              282.53048            12.409195
#> 41   HNRNPU          683              594.17494            15.610278
#> 42  IGF2BP1         1037             1026.33820            15.492784
#> 43  IGF2BP2         1015              981.57682            15.754263
#> 47    MBNL2          520              437.60129            14.427014
#> 49   METTL3          190              203.64101            10.843020
#> 50    MOV10         1262             1222.42815            13.706518
#> 55    PTBP1          643              608.93825            15.696141
#> 57     PUM2          550              487.19631            14.924347
#> 58      QKI          562              469.02117            14.731150
#> 62    RBM47          657              545.00449            15.389697
#> 63    RBPMS          341              342.26143            13.353425
#> 71   TARDBP          372              351.95840            13.427179
#> 72     TIA1          731              711.00590            16.122499
#> 74   TNRC6A           62               57.38572             6.040078
#> 81   YTHDC1         1119             1215.25222            13.835875
#> 83   YTHDF1         1204             1287.67363            12.852510
#> 86   ZC3H7B         1173             1125.08914            14.881164
#> 510   DDX24          379              468.08897            14.783998
#> 91    NCBP2          419              433.87516            14.431640
#> 101 FAM120A          534              582.90766            15.545144
#> 141   RBM22          424              455.43362            14.588606
#> 191    SUB1          298              324.33933            13.062773
#> 201   TBRG4          484              472.18198            14.767771
#> 241  METAP2          131              169.13287             9.989987
#> 281    SLTM          164              199.16709            10.718851
#> 381  AKAP8L           97              150.01127             9.457502
#> 401   PCBP2          572              643.72831            15.897439
#> 491 SUPV3L1           66               85.45194             7.290766
#> 501   SRSF7           93              116.15368             8.437486
#> 521  GTF2F1          351              398.69558            14.058368
#> 621   TRA2A          152              159.56321             9.736560
#> 671   KHSRP          618              541.29419            15.290466
#> 691   LARP4          311              390.83063            13.923736
#> 711   LSM11          195              256.56345            11.939002
#> 821  ZNF622          507              533.08197            15.298377
#> 871    YBX3          559              595.58986            15.642586
#> 96   ZRANB2          178              186.21398            10.426727
#> 99    SAFB2           89               73.09560             6.780617
#> 109    AARS          155              199.91384            10.730626
#> 112  EIF4G2          187              209.80694            10.913050
#> 116   FKBP4           48               73.41559             6.798841
#>               z    pval         padj
#> 6    0.35686030 0.34869 0.9371043750
#> 8   -2.98402378 0.99836 1.0000000000
#> 11  -3.11253060 0.99894 1.0000000000
#> 15  -3.43997101 0.99972 1.0000000000
#> 22  -2.69678975 0.99621 1.0000000000
#> 26   1.92414281 0.02446 0.1051780000
#> 31  -0.91308345 0.80694 1.0000000000
#> 37   5.53565657 0.00001 0.0001075000
#> 38   4.79237548 0.00001 0.0001075000
#> 41   5.69016523 0.00001 0.0001075000
#> 42   0.68817847 0.23577 0.7241507143
#> 43   2.12153245 0.01556 0.0743422222
#> 47   5.71141806 0.00001 0.0001075000
#> 49  -1.25804522 0.88792 1.0000000000
#> 50   2.88708262 0.00178 0.0109342857
#> 55   2.17007162 0.01432 0.1523127273
#> 57   4.20813648 0.00003 0.0002866667
#> 58   6.31171577 0.00001 0.0001075000
#> 62   7.27730435 0.00001 0.0001075000
#> 63  -0.09446491 0.52123 1.0000000000
#> 71   1.49261436 0.06364 0.4653675000
#> 72   1.24013651 0.10174 0.3645683333
#> 74   0.76394375 0.19850 0.6322592593
#> 81  -6.95671357 1.00000 1.0000000000
#> 83  -6.51029505 1.00000 1.0000000000
#> 86   3.21956410 0.00061 0.0040353846
#> 510 -6.02604041 1.00000 1.0000000000
#> 91  -1.03073246 0.84044 1.0000000000
#> 101 -3.14616958 0.99896 1.0000000000
#> 141 -2.15466921 0.98321 1.0000000000
#> 191 -2.01636582 0.97642 1.0000000000
#> 201  0.80025753 0.20243 0.9473724000
#> 241 -3.81710920 0.99995 1.0000000000
#> 281 -3.28086388 0.99947 1.0000000000
#> 381 -5.60520864 1.00000 1.0000000000
#> 401 -4.51194107 1.00000 1.0000000000
#> 491 -2.66802436 0.99580 1.0000000000
#> 501 -2.74414432 0.99669 1.0000000000
#> 521 -3.39268262 0.99964 1.0000000000
#> 621 -0.77678465 0.76634 1.0000000000
#> 671  5.01657755 0.00001 0.0001950000
#> 691 -5.73342028 1.00000 1.0000000000
#> 711 -5.15649865 1.00000 1.0000000000
#> 821 -1.70488473 0.95282 1.0000000000
#> 871 -2.33911838 0.98983 1.0000000000
#> 96  -0.78778123 0.77040 1.0000000000
#> 99   2.34556842 0.00838 0.0980460000
#> 109 -4.18557512 0.99999 1.0000000000
#> 112 -2.08987773 0.98017 1.0000000000
#> 116 -3.73822398 0.99993 1.0000000000
#> 
#> $exclusive_up
#>         RBP real_overlap simulated_overlap_mean simulated_overlap_sd          z
#> 6     CELF2          126               88.24509             6.621883  5.7015365
#> 8     CPSF1          156              131.05240             5.957457  4.1876260
#> 11    CPSF4          167              129.92143             5.992079  6.1879303
#> 15   CSTF2T          179              157.83387             4.499012  4.7046169
#> 22    EIF3D           88               85.54499             6.602425  0.3718346
#> 26    EWSR1          164              115.55597             6.373302  7.6010886
#> 31      FUS          178              152.38571             4.879579  5.2492823
#> 37   HNRNPD          162              105.46488             6.525784  8.6633457
#> 38   HNRNPF           57               32.16059             5.049419  4.9192607
#> 41   HNRNPU           99               67.62095             6.380647  4.9178475
#> 42  IGF2BP1          163              116.83672             6.369066  7.2480452
#> 43  IGF2BP2          167              111.73327             6.448971  8.5698531
#> 47    MBNL2           67               49.82545             5.906844  2.9075677
#> 49   METTL3           29               23.18325             4.424559  1.3146507
#> 50    MOV10          174              139.14097             5.627818  6.1940579
#> 55    PTBP1          180              159.67048             4.344876  4.6789646
#> 57     PUM2           86               55.45614             6.075540  5.0273491
#> 58      QKI           88               53.41721             6.028380  5.7366636
#> 62    RBM47          100               62.04289             6.259300  6.0641138
#> 63    RBPMS           61               38.96723             5.422098  4.0635136
#> 71   TARDBP          180              168.65969             3.450148  3.2869058
#> 72     TIA1          119               80.90386             6.579405  5.7902105
#> 74   TNRC6A            7                6.54626             2.461906  0.1843044
#> 81   YTHDC1          169              138.31190             5.643414  5.4378612
#> 83   YTHDF1          176              146.54885             5.244205  5.6159415
#> 86   ZC3H7B          174              128.06706             6.061862  7.5773649
#> 510   DDX24           54               53.30842             6.030090  0.1146882
#> 91    NCBP2           51               49.37618             5.888959  0.2757397
#> 101 FAM120A           77               66.38116             6.393045  1.6609988
#> 141   RBM22           42               51.84665             5.981469 -1.6461926
#> 191    SUB1           43               36.91274             5.303166  1.1478540
#> 201   TBRG4           47               53.78239             6.030827 -1.1246201
#> 241  METAP2            7               19.25905             4.071097 -3.0112401
#> 281    SLTM           10               22.67154             4.370106 -2.8995957
#> 381  AKAP8L            6               17.08318             3.859825 -2.8714202
#> 401   PCBP2           62               73.26724             6.480634 -1.7386014
#> 491 SUPV3L1            2                9.73388             2.988121 -2.5882081
#> 501   SRSF7           11               13.22244             3.427535 -0.6484076
#> 521  GTF2F1           23               45.40975             5.730772 -3.9104245
#> 621   TRA2A           23               18.17279             3.971779  1.2153772
#> 671   KHSRP          107               61.60886             6.237110  7.2775920
#> 691   LARP4           32               44.53398             5.656558 -2.2158317
#> 711   LSM11           14               29.24427             4.865565 -3.1330937
#> 821  ZNF622           72               60.69684             6.225528  1.8156146
#> 871    YBX3           83               67.79643             6.422898  2.3670887
#> 96   ZRANB2           17               21.21564             4.264395 -0.9885671
#> 99    SAFB2            5                8.31952             2.765520 -1.2003238
#> 109    AARS           18               22.76861             4.378349 -1.0891343
#> 112  EIF4G2           23               23.88773             4.483511 -0.1979989
#> 116   FKBP4            3                8.36430             2.762477 -1.9418439
#>        pval         padj
#> 6   0.00001 1.622642e-05
#> 8   0.00001 1.622642e-05
#> 11  0.00001 1.622642e-05
#> 15  0.00001 1.622642e-05
#> 22  0.32714 3.695766e-01
#> 26  0.00001 1.622642e-05
#> 31  0.00001 1.622642e-05
#> 37  0.00001 1.622642e-05
#> 38  0.00001 1.622642e-05
#> 41  0.00001 1.622642e-05
#> 42  0.00001 1.622642e-05
#> 43  0.00001 1.622642e-05
#> 47  0.00172 2.595088e-03
#> 49  0.08058 9.966171e-02
#> 50  0.00001 1.622642e-05
#> 55  0.00001 1.622642e-05
#> 57  0.00001 1.622642e-05
#> 58  0.00001 1.622642e-05
#> 62  0.00001 1.622642e-05
#> 63  0.00001 1.622642e-05
#> 71  0.00001 1.622642e-05
#> 72  0.00001 1.622642e-05
#> 74  0.33090 3.695766e-01
#> 81  0.00001 1.622642e-05
#> 83  0.00001 1.622642e-05
#> 86  0.00001 1.622642e-05
#> 510 0.41740 9.999800e-01
#> 91  0.35528 9.999800e-01
#> 101 0.04180 2.574000e-01
#> 141 0.94399 9.999800e-01
#> 191 0.10824 5.065632e-01
#> 201 0.85074 9.999800e-01
#> 241 0.99939 9.999800e-01
#> 281 0.99883 9.999800e-01
#> 381 0.99891 9.999800e-01
#> 401 0.95209 9.999800e-01
#> 491 0.99746 9.999800e-01
#> 501 0.68273 9.999800e-01
#> 521 0.99995 9.999800e-01
#> 621 0.09314 4.737991e-01
#> 671 0.00001 3.900000e-04
#> 691 0.98577 9.999800e-01
#> 711 0.99939 9.999800e-01
#> 821 0.02937 2.315820e-01
#> 871 0.00773 1.004900e-01
#> 96  0.80606 9.999800e-01
#> 99  0.84735 9.999800e-01
#> 109 0.83406 9.999800e-01
#> 112 0.52337 9.999800e-01
#> 116 0.97248 9.999800e-01
#> 
#> $exclusive_down
#>         RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 6     CELF2          112              136.72033             8.127101
#> 8     CPSF1          202              203.06069             7.307140
#> 11    CPSF4          189              201.33701             7.380965
#> 15   CSTF2T          114               75.61940             7.182994
#> 22    EIF3D          171              132.52318             8.132937
#> 26    EWSR1           57               25.45226             4.650234
#> 31      FUS          222              236.15512             6.029624
#> 37   HNRNPD          104              163.42362             8.058828
#> 38   HNRNPF           22               49.85540             6.214298
#> 41   HNRNPU           10                9.46408             2.942644
#> 42  IGF2BP1           82               88.94059             7.519792
#> 43  IGF2BP2           52               51.79202             6.313464
#> 47    MBNL2           66               77.23079             7.268994
#> 49   METTL3           40               35.90435             5.437420
#> 50    MOV10          184              215.60573             6.932793
#> 55    PTBP1          115              107.40036             7.860562
#> 57     PUM2           62               85.93277             7.468892
#> 58      QKI           91               95.39361             7.701885
#> 62    RBM47           38               96.11678             7.717411
#> 63    RBPMS           52               60.41977             6.674070
#> 71   TARDBP          264              261.32311             4.249921
#> 72     TIA1           80               85.28731             7.483939
#> 74   TNRC6A           13               10.12716             3.036719
#> 81   YTHDC1          230              214.31582             6.993477
#> 83   YTHDF1          238              227.09690             6.462819
#> 86   ZC3H7B          141              198.43642             7.476296
#> 510   DDX24          121               82.58047             7.407222
#> 91    NCBP2           73               76.56581             7.242115
#> 101 FAM120A          121              102.78765             7.818864
#> 141   RBM22          115               80.33876             7.361156
#> 191    SUB1           76               57.19697             6.526333
#> 201   TBRG4          103               83.23211             7.421250
#> 241  METAP2           67               29.83893             5.023523
#> 281    SLTM           55               35.12678             5.361831
#> 381  AKAP8L           55               26.47226             4.739417
#> 401   PCBP2          167              113.52356             7.995274
#> 491 SUPV3L1           35               15.06406             3.664747
#> 501   SRSF7           38               20.48714             4.199899
#> 521  GTF2F1          117               70.32837             7.032703
#> 621   TRA2A           28               28.09618             4.885156
#> 671   KHSRP           52               95.44124             7.700317
#> 691   LARP4          124               68.91859             6.986589
#> 711   LSM11           81               45.22624             5.977516
#> 821  ZNF622           85               93.97778             7.653638
#> 871    YBX3          124              105.01701             7.856777
#> 96   ZRANB2           55               32.84997             5.209241
#> 99    SAFB2           10               12.89409             3.397689
#> 109    AARS           52               35.27258             5.377279
#> 112  EIF4G2           53               37.01162             5.490012
#> 116   FKBP4           31               12.95467             3.407922
#>               z    pval         padj
#> 6   -3.04171551 0.99875 1.000000e+00
#> 8   -0.14515802 0.53335 1.000000e+00
#> 11  -1.67146309 0.94387 1.000000e+00
#> 15   5.34325954 0.00001 3.078947e-05
#> 22   4.73098710 0.00001 8.600000e-04
#> 26   6.78411868 0.00001 3.078947e-05
#> 31  -2.34759594 0.98596 1.000000e+00
#> 37  -7.37372939 1.00000 1.000000e+00
#> 38  -4.48246945 1.00000 1.000000e+00
#> 41   0.18212190 0.34786 4.472486e-01
#> 42  -0.92297637 0.80301 8.947826e-01
#> 43   0.03294230 0.44879 5.527203e-01
#> 47  -1.54502672 0.93120 1.000000e+00
#> 49   0.75323406 0.19710 7.062750e-01
#> 50  -4.55887389 0.99999 1.000000e+00
#> 55   0.96680623 0.15130 2.185444e-01
#> 57  -3.20432669 0.99935 1.000000e+00
#> 58  -0.57045904 0.69138 8.089146e-01
#> 62  -7.53060545 1.00000 1.000000e+00
#> 63  -1.26156456 0.88286 1.000000e+00
#> 71   0.62986814 0.23325 8.023800e-01
#> 72  -0.70648759 0.73815 8.384811e-01
#> 74   0.94603418 0.13418 5.668219e-01
#> 81   2.24268686 0.00841 9.040750e-02
#> 83   1.68705021 0.03583 3.423756e-01
#> 86  -7.68247040 1.00000 1.000000e+00
#> 510  5.18676659 0.00001 3.078947e-05
#> 91  -0.49237134 0.66126 7.814891e-01
#> 101  2.32928349 0.00879 1.534970e-02
#> 141  4.70866779 0.00001 3.078947e-05
#> 191  2.88110196 0.00215 4.413158e-03
#> 201  2.66368753 0.00325 6.337500e-03
#> 241  7.39741280 0.00001 3.078947e-05
#> 281  3.70642400 0.00018 4.680000e-04
#> 381  6.01925067 0.00001 3.078947e-05
#> 401  6.68850635 0.00001 3.078947e-05
#> 491  5.43992258 0.00001 3.078947e-05
#> 501  4.16982892 0.00003 8.357143e-05
#> 521  6.63637112 0.00001 3.078947e-05
#> 621 -0.01968822 0.45741 5.574684e-01
#> 671 -5.64148731 1.00000 1.000000e+00
#> 691  7.88387686 0.00001 3.078947e-05
#> 711  5.98472047 0.00001 3.078947e-05
#> 821 -1.17300821 0.86663 9.476235e-01
#> 871  2.41612955 0.00668 1.202400e-02
#> 96   4.25206470 0.00001 3.078947e-05
#> 99  -0.85178198 0.75318 8.473275e-01
#> 109  3.11075919 0.00109 2.475000e-03
#> 112  2.91226712 0.00187 3.978000e-03
#> 116  5.29511262 0.00001 3.078947e-05
#> 
#> $buffered_up
#>         RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 6     CELF2          340              232.17709            10.370333
#> 8     CPSF1          421              344.89387             9.323576
#> 11    CPSF4          424              341.93605             9.389603
#> 15   CSTF2T          463              415.43166             7.031960
#> 22    EIF3D          232              225.10011            10.339019
#> 26    EWSR1          435              304.06293            10.008104
#> 31      FUS          463              401.08502             7.665871
#> 37   HNRNPD          418              277.50034            10.250081
#> 38   HNRNPF          158               84.64139             7.946870
#> 41   HNRNPU          295              177.88711            10.030460
#> 42  IGF2BP1          420              307.43971             9.976256
#> 43  IGF2BP2          429              294.01620            10.125674
#> 47    MBNL2          216              131.08386             9.239231
#> 49   METTL3           62               60.95483             6.918187
#> 50    MOV10          453              366.16674             8.826006
#> 55    PTBP1          458              420.24166             6.816568
#> 57     PUM2          259              145.89908             9.552049
#> 58      QKI          229              140.48793             9.435312
#> 62    RBM47          300              163.22249             9.823427
#> 63    RBPMS          179              102.58011             8.528547
#> 71   TARDBP          468              443.87558             5.416901
#> 72     TIA1          303              212.92048            10.300044
#> 74   TNRC6A           24               17.18610             3.854918
#> 81   YTHDC1          439              363.99954             8.895131
#> 83   YTHDF1          444              385.70228             8.228098
#> 86   ZC3H7B          450              336.98766             9.465652
#> 510   DDX24          112              140.22452             9.403460
#> 91    NCBP2          172              129.95402             9.249347
#> 101 FAM120A          229              174.55645             9.961650
#> 141   RBM22          131              136.43467             9.343599
#> 191    SUB1           94               97.12964             8.301790
#> 201   TBRG4          141              141.41669             9.448449
#> 241  METAP2           26               50.64882             6.359272
#> 281    SLTM           38               59.61607             6.869662
#> 381  AKAP8L           27               44.91586             6.055956
#> 401   PCBP2          209              192.80094            10.196989
#> 491 SUPV3L1           14               25.59212             4.641977
#> 501   SRSF7           31               34.79332             5.376122
#> 521  GTF2F1          109              119.40971             8.934076
#> 621   TRA2A           43               47.81325             6.225172
#> 671   KHSRP          278              162.10403             9.804512
#> 691   LARP4           75              117.08972             8.891448
#> 711   LSM11           56               76.88667             7.600112
#> 821  ZNF622          202              159.68602             9.779716
#> 871    YBX3          195              178.34693            10.026870
#> 96   ZRANB2           45               55.80061             6.646150
#> 99    SAFB2           29               21.90376             4.337325
#> 109    AARS           53               59.83428             6.878474
#> 112  EIF4G2           60               62.82293             7.009466
#> 116   FKBP4            7               21.99365             4.334792
#>               z    pval         padj
#> 6   10.39724674 0.00001 1.365079e-05
#> 8    8.16276150 0.00001 1.365079e-05
#> 11   8.73987443 0.00001 1.365079e-05
#> 15   6.76459176 0.00001 1.365079e-05
#> 22   0.66736408 0.23728 2.685011e-01
#> 26  13.08310408 0.00001 1.365079e-05
#> 31   8.07670510 0.00001 1.365079e-05
#> 37  13.70717623 0.00001 1.365079e-05
#> 38   9.23113268 0.00001 1.365079e-05
#> 41  11.67572450 0.00001 1.365079e-05
#> 42  11.28281912 0.00001 1.365079e-05
#> 43  13.33084535 0.00001 1.365079e-05
#> 47   9.19082380 0.00001 1.365079e-05
#> 49   0.15107571 0.40333 4.230046e-01
#> 50   9.83834170 0.00001 1.365079e-05
#> 55   5.53920061 0.00001 1.365079e-05
#> 57  11.84048735 0.00001 1.365079e-05
#> 58   9.38093769 0.00001 1.365079e-05
#> 62  13.92360461 0.00001 1.365079e-05
#> 63   8.96048204 0.00001 1.365079e-05
#> 71   4.45354666 0.00001 1.365079e-05
#> 72   8.74554747 0.00001 1.365079e-05
#> 74   1.76758608 0.03355 4.121857e-02
#> 81   8.43163098 0.00001 1.365079e-05
#> 83   7.08519984 0.00001 1.365079e-05
#> 86  11.93920319 0.00001 1.365079e-05
#> 510 -3.00150379 0.99860 1.000000e+00
#> 91   4.54583222 0.00001 6.882353e-05
#> 101  5.46531462 0.00001 6.882353e-05
#> 141 -0.58164631 0.70003 1.000000e+00
#> 191 -0.37698374 0.62137 1.000000e+00
#> 201 -0.04410142 0.49323 1.000000e+00
#> 241 -3.87604446 0.99995 1.000000e+00
#> 281 -3.14659865 0.99935 1.000000e+00
#> 381 -2.95838659 0.99894 1.000000e+00
#> 401  1.58861216 0.05133 2.070900e-01
#> 491 -2.49723777 0.99425 1.000000e+00
#> 501 -0.70558661 0.72599 1.000000e+00
#> 521 -1.16516924 0.86673 1.000000e+00
#> 621 -0.77319153 0.75250 1.000000e+00
#> 671 11.82067745 0.00001 6.882353e-05
#> 691 -4.73373086 1.00000 1.000000e+00
#> 711 -2.74820552 0.99713 1.000000e+00
#> 821  4.32670828 0.00001 6.882353e-05
#> 871  1.66084437 0.04372 1.857375e-01
#> 96  -1.62509262 0.94205 1.000000e+00
#> 99   1.63608685 0.04445 1.857375e-01
#> 109 -0.99357497 0.82110 1.000000e+00
#> 112 -0.40273113 0.62486 1.000000e+00
#> 116 -3.45890913 0.99992 1.000000e+00
#> 
#> $buffered_down
#>         RBP real_overlap simulated_overlap_mean simulated_overlap_sd
#> 6     CELF2          230              275.37189            11.137038
#> 8     CPSF1          408              409.01211             9.996747
#> 11    CPSF4          408              405.47406            10.081826
#> 15   CSTF2T          195              152.32679             9.909798
#> 22    EIF3D           53               35.04589             5.372857
#> 26    EWSR1          101               51.19384             6.403497
#> 31      FUS          458              475.59383             8.220012
#> 37   HNRNPD          196              329.08036            11.034596
#> 38   HNRNPF           39              100.33219             8.527848
#> 41   HNRNPU           26               19.06871             4.019864
#> 42  IGF2BP1          305              364.49864            10.719282
#> 43  IGF2BP2          123              104.25928             8.639677
#> 47    MBNL2           85              155.47089             9.957203
#> 49   METTL3           95               72.30757             7.448671
#> 50    MOV10          360              434.15968             9.479745
#> 55    PTBP1          511              498.30207             7.315136
#> 57     PUM2           80              173.07125            10.309023
#> 58      QKI           75              166.61606            10.189824
#> 62    RBM47           74              193.55307            10.580215
#> 63    RBPMS           90              121.57266             9.157484
#> 71   TARDBP          536              526.36013             5.817251
#> 72     TIA1          166              171.72040            10.252619
#> 74   TNRC6A           14               20.39064             4.147249
#> 81   YTHDC1          469              431.68272             9.546843
#> 83   YTHDF1          485              457.40982             8.808174
#> 86   ZC3H7B          304              399.60191            10.193436
#> 510   DDX24          286              166.22980            10.206129
#> 91    NCBP2          158              154.06702             9.912054
#> 101 FAM120A          220              207.09524            10.771692
#> 141   RBM22          205              161.74891            10.096698
#> 191    SUB1          149              115.14156             8.990659
#> 201   TBRG4          194              167.73783            10.156027
#> 241  METAP2          121               60.07804             6.900153
#> 281    SLTM          108               70.74760             7.370028
#> 381  AKAP8L          112               53.26130             6.518450
#> 401   PCBP2          285              228.63631            10.942980
#> 491 SUPV3L1           52               30.33546             5.024195
#> 501   SRSF7           49               41.28290             5.780677
#> 521  GTF2F1          202              141.69659             9.701949
#> 621   TRA2A           67               56.65426             6.677279
#> 671   KHSRP          110              192.18659            10.573335
#> 691   LARP4          242              138.82322             9.628070
#> 711   LSM11          153               91.13664             8.203693
#> 821  ZNF622          219              189.27374            10.531599
#> 871    YBX3          244              211.52658            10.844772
#> 96   ZRANB2           91               66.09523             7.176616
#> 99    SAFB2           20               25.94594             4.665706
#> 109    AARS          112               70.98194             7.378227
#> 112  EIF4G2          104               74.50297             7.532778
#> 116   FKBP4           50               26.08355             4.677304
#>               z    pval         padj
#> 6    -4.0739640 0.99997 1.000000e+00
#> 8    -0.1012439 0.52261 1.000000e+00
#> 11    0.2505439 0.38533 8.720626e-01
#> 15    4.3061634 0.00001 2.600000e-05
#> 22    3.3416318 0.00058 1.233818e-03
#> 26    7.7779624 0.00001 2.600000e-05
#> 31   -2.1403656 0.97958 1.000000e+00
#> 37  -12.0602833 1.00000 1.000000e+00
#> 38   -7.1919893 1.00000 1.000000e+00
#> 41    1.7242596 0.03652 5.622158e-02
#> 42   -5.5506180 1.00000 1.000000e+00
#> 43    2.1691459 0.01408 2.288000e-02
#> 47   -7.0773783 1.00000 1.000000e+00
#> 49    3.0465071 0.00120 9.381818e-03
#> 50   -7.8229617 1.00000 1.000000e+00
#> 55    1.7358433 0.03261 1.168525e-01
#> 57   -9.0281348 1.00000 1.000000e+00
#> 58   -8.9909367 1.00000 1.000000e+00
#> 62  -11.2996828 1.00000 1.000000e+00
#> 63   -3.4477438 0.99979 1.000000e+00
#> 71    1.6571176 0.03697 1.217696e-01
#> 72   -0.5579452 0.69458 8.208673e-01
#> 74   -1.5409348 0.92757 1.000000e+00
#> 81    3.9088609 0.00001 1.433333e-04
#> 83    3.1323383 0.00047 4.042000e-03
#> 86   -9.3787720 1.00000 1.000000e+00
#> 510  11.7351248 0.00001 2.600000e-05
#> 91    0.3967876 0.32473 4.246710e-01
#> 101   1.1980253 0.10711 1.566484e-01
#> 141   4.2836868 0.00003 7.468085e-05
#> 191   3.7659575 0.00011 2.523529e-04
#> 201   2.5858704 0.00460 8.542857e-03
#> 241   8.8290743 0.00001 2.600000e-05
#> 281   5.0545803 0.00001 2.600000e-05
#> 381   9.0111456 0.00001 2.600000e-05
#> 401   5.1506709 0.00001 2.600000e-05
#> 491   4.3120417 0.00002 5.086957e-05
#> 501   1.3349821 0.08042 1.206300e-01
#> 521   6.2155975 0.00001 2.600000e-05
#> 621   1.5493947 0.05476 8.320675e-02
#> 671  -7.7730054 1.00000 1.000000e+00
#> 691  10.7162472 0.00001 2.600000e-05
#> 711   7.5409158 0.00001 2.600000e-05
#> 821   2.8225781 0.00229 4.392295e-03
#> 871   2.9943849 0.00126 2.586316e-03
#> 96    3.4702665 0.00030 6.622642e-04
#> 99   -1.2743923 0.88068 1.000000e+00
#> 109   5.5593381 0.00001 2.600000e-05
#> 112   3.9158238 0.00011 2.523529e-04
#> 116   5.1132987 0.00001 2.600000e-05
```

#### Step 3: Visualisation

Generate and save heatmap to pdf. The heatmap represents the -logP of
each RBP for each gene group. The blue RBPs are downregulated and the
orange RBPs are upregulated. Only RBPs that are significant in at least
one gene group are shown.

``` r
# Heatmap of RBP scores

HeatmapRBP(res=res,rbp_lfc=rbp_lfc)
```

![](man/figures/README-Generate%20heatmap-1.png)<!-- -->

``` r

# If there is not at least 1 positive and 1 negative RBP lfc then use the heatmap for miRNA
# HeatmapmiRNA(res=res)
```

Suggestion for saving your heatmap

``` r
# Save the heatmap
e=HeatmapRBP(res=res,rbp_lfc=rbp_lfc)
location="Heatmap_fibroblasts.pdf"
n=length(e$tree_row$order)
pdf(location,length(names(res)),3+n*0.15)
e
dev.off()
```

A bubble plot can be generated only if the gene groups are the output of
DeltaTE program.

``` r
# Bubble plot gene_groups if gene_groups are from DeltaTE. FDR has to be lower or equal to the FDR put in CLIPreg::combine() step
BubbleRBPs(res = res,gene_groups = gene_groups,rbp_lfc = rbp_lfc,FDR=0.1)
```

![](man/figures/README-Bubble%20plot-1.png)<!-- -->

From the results, the user can choose a number of RBP to draw the
network for by n. This will pick the n most changing RBPs.

``` r
# Draw network of the RBP that are most changing or choose specific RBPs
Draw_network_by_group(rbp_lfc=rbp_lfc,res=res,Targets=Targets,gene_groups=gene_groups,n=5,forwarded = F)
```

![](man/figures/README-Network-1.png)<!-- -->

``` r

# You can subset your RBP_lfc to keep only your RBPs of interest
```

``` r
# plot GO
Plot_GO(rbp_lfc=rbp_lfc,res=res,Targets=Targets,gene_groups=gene_groups,n=5,
  tpm_ribo = tpm_ribo,th=200,GO_to_show=3,forwarded = F)
```

![](man/figures/README-GO%20of%20specific%20nodes-1.png)<!-- -->

``` r
Plot_GO_node_name(rbp_lfc=rbp_lfc,res=res,Targets=Targets,gene_groups=gene_groups,n=5,
                  tpm_ribo = tpm_ribo,Nodes_to_keep=c(19,15),GO_to_show=3,forwarded = F)
```

![](man/figures/README-GO%20of%20specific%20nodes-2.png)<!-- -->

``` r
Plot_GO_RBP(rbp_of_interest="QKI",tpm_ribo = tpm_ribo,Targets=Targets,gene_groups=gene_groups,GO_to_show=3)
```

![](man/figures/README-GO%20of%20specific%20nodes-3.png)<!-- -->

### Analysis of miRNA target enrichment

The same analysis can be applied to run miRNA target enrichment. Here is
an example of the code that is very similar to the steps followed for
the RBPs.

``` r
data("miR_data") # preparation can be found in Example/data_processing.R
data("gene_groups") # for example on fibroblasts
Targets=GetTarget(RBP_data=miR_data,background=gene_groups$geneID)
data("ribo_lfc")
data("tpm_ribo")

res_miR=CLIPreg(RBP_data=miR_data,gene_groups=gene_groups) # Takes several minutes
save(res_miR,file="Res_miR.RData")
res=CLIPreg::combine(res1=res_miR,res2 = res_miR,FDR=0.1)
HeatmapmiRNA(res=res)


# Save the heatmap
e=HeatmapmiRNA(res=res)
location="Heatmap_HeLa_EGF_miRNA.pdf"
n=nrow(e@matrix)
pdf(location,length(names(res))+3,3+n*0.15)
e
dev.off()


# Network and GO
Draw_network_by_group(rbp_lfc=c("hsa-miR-133a-3p.1","hsa-miR-590-5p","hsa-miR-29a-3p","hsa-miR-30a-5p","hsa-miR-155-5p"),
                      res=res,Targets=Targets,gene_groups=gene_groups,n=5,forwarded = F)
Plot_GO_RBP(rbp_of_interest="hsa-miR-133a-3p.1",tpm_ribo = tpm_ribo,Targets=Targets,gene_groups=gene_groups,GO_to_show=3)
```

### Another example

HeLa stimulated by EGF. groups were obtained by categorizing into TE up
and down from data given in the paper. Data preparation can be found in
Example/data\_processing.R on github.

``` r
data("gene_groups_HeLa_EGF")
data("ribo_lfc_HeLa_EGF")
data("tpm_ribo")

data("RBP_ENCODE")
data("RBP_POSTAR")
Targets=combine_targets(RBP_list1=RBP_ENCODE,RBP_list2=RBP_POSTAR,background=gene_groups$geneID)

res=CLIPreg::combine(res1=res_Encode,res2=res_Postar,FDR=0.1)
rbp_lfc=rbp_change(res=res,ribo_lfc=ribo_lfc)
res=cure_res(res=res,rbp_lfc=rbp_lfc)

HeatmapmiRNA(res=res)
Draw_network_by_group(rbp_lfc=rbp_lfc[c("GRWD1","GRSF1","SUGP2","DHX30","LSM11"),],res=res,Targets=Targets,gene_groups=gene_groups,n=5,forwarded = F)
Plot_GO_RBP(rbp_of_interest="GRWD1",tpm_ribo = tpm_ribo,Targets=Targets,gene_groups=gene_groups,GO_to_show=3)
```
