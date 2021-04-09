---
title: "analyze_go_results"
output:
  html_document:
    keep_md: yes
---



## Analyze go results

```r
go1 = read.csv('data/GOlists/GOresults_plastic_expasref.csv')
go2 = read.csv('data/GOlists/GOresults_evolved_expasref.csv')

##merge
goall = dplyr::left_join(go1, go2, by = "GO.biological.process.complete")


##plastic

plastictable = data.frame(bio= goall$GO.biological.process.complete, enrich = goall$plastic_ref_forGO.txt..fold.Enrichment., observed= goall$plastic_ref_forGO.txt..7679., expected = goall$plastic_ref_forGO.txt..expected., pval = goall$plastic_ref_forGO.txt..raw.P.value., fdr=goall$plastic_ref_forGO.txt..FDR.)
dplyr::filter(plastictable, pval < 0.01, enrich>1) ##only want enrichments not depletions. ##lowering p value cutoff because nothing sig at 0.01
```

```
##                                                       bio enrich observed
## 1                             photosynthesis (GO:0015979)   1.41      117
## 2         defense response to other organism (GO:0098542)    1.2      337
## 3 interspecies interaction between organisms (GO:0044419)   1.19      455
## 4                response to biotic stimulus (GO:0009607)   1.18      447
## 5       response to external biotic stimulus (GO:0043207)   1.18      446
## 6                 response to other organism (GO:0051707)   1.18      446
##   expected    pval fdr
## 1    83.12 0.00633   1
## 2   279.87 0.00833   1
## 3   383.92 0.00476   1
## 4   377.94 0.00545   1
## 5   377.34 0.00595   1
## 6   377.34 0.00595   1
```

```r
##plastic32
plastic32table = data.frame(bio= goall$GO.biological.process.complete, enrich = goall$plastic32_ref_forGO.txt..fold.Enrichment., observed= goall$plastic32_ref_forGO.txt..7279., expected = goall$plastic32_ref_forGO.txt..expected., pval = goall$plastic32_ref_forGO.txt..raw.P.value., fdr = goall$plastic32_ref_forGO.txt..FDR.)
dplyr::filter(plastic32table, pval<0.001, enrich>1)
```

```
##                                 bio enrich observed expected     pval   fdr
## 1 response to chemical (GO:0042221)   1.13     1029   907.54 0.000802 0.286
```

```r
##plastic8

plastic8table = data.frame(bio= goall$GO.biological.process.complete, enrich = goall$plastic8_ref_forGO.txt..fold.Enrichment., observed= goall$plastic8_ref_forGO.txt..4605., expected = goall$plastic8_ref_forGO.txt..expected., pval = goall$plastic8_ref_forGO.txt..raw.P.value., fdr=goall$plastic8_ref_forGO.txt..FDR.)
dplyr::filter(plastic8table, pval<0.001, enrich>1)
```

```
##                                                            bio enrich observed
## 1                                  photosynthesis (GO:0015979)   1.97       98
## 2                  photosynthesis, light reaction (GO:0019684)      2       68
## 3                              response to chitin (GO:0010200)   1.93       52
## 4                            plastid organization (GO:0009657)   1.46      121
## 5  generation of precursor metabolites and energy (GO:0006091)   1.49      130
## 6              defense response to other organism (GO:0098542)   1.33      224
## 7      interspecies interaction between organisms (GO:0044419)   1.29      298
## 8                     response to biotic stimulus (GO:0009607)   1.29      293
## 9            response to external biotic stimulus (GO:0043207)   1.29      293
## 10                     response to other organism (GO:0051707)   1.29      293
## 11               small molecule metabolic process (GO:0044281)   1.22      418
## 12         cellular response to chemical stimulus (GO:0070887)   1.29      302
## 13                   response to abiotic stimulus (GO:0009628)   1.18      569
## 14                           response to chemical (GO:0042221)   1.18      677
## 15                           response to stimulus (GO:0050896)    1.1     1256
##    expected     pval     fdr
## 1     49.85 5.17e-07 0.00295
## 2     34.07 2.23e-05 0.02550
## 3     26.90 3.62e-04 0.10300
## 4     82.84 9.38e-04 0.20600
## 5     87.50 2.91e-04 0.09220
## 6    167.83 3.59e-04 0.10800
## 7    230.23 1.91e-04 0.09070
## 8    226.65 2.31e-04 0.08240
## 9    226.29 2.28e-04 0.09280
## 10   226.29 2.28e-04 0.08670
## 11   343.56 5.75e-04 0.14900
## 12   233.46 1.80e-04 0.09350
## 13   481.98 5.30e-04 0.14400
## 14   574.15 1.31e-04 0.07490
## 15  1139.68 7.97e-04 0.18200
```

```r
##different dir

difftable = data.frame(bio= goall$GO.biological.process.complete, 
                       enrich = goall$different_dir_832_genes_ref_forGO.txt..fold.Enrichment.,
                       observed= goall$different_dir_832_genes_ref_forGO.txt..129., 
                       expected = goall$different_dir_832_genes_ref_forGO.txt..expected., 
                       pval = goall$different_dir_832_genes_ref_forGO.txt..raw.P.value.,
                       fdr = goall$different_dir_832_genes_ref_forGO.txt..FDR.)
dplyr::filter(difftable, pval < 0.001, enrich>1)
```

```
##                                           bio enrich observed expected     pval
## 1         cell redox homeostasis (GO:0045454)   9.39        5     0.53 2.71e-04
## 2 proton transmembrane transport (GO:1902600)   9.76        5     0.51 2.29e-04
## 3           sulfate assimilation (GO:0000103)  33.18        4     0.12 1.55e-05
## 4              sulfate reduction (GO:0019419)  66.36        2     0.03 9.63e-04
##      fdr
## 1 0.5150
## 2 0.6540
## 3 0.0885
## 4 1.0000
```

```r
##evolved
evolvedtable = data.frame(bio= goall$GO.biological.process.complete, enrich = goall$evolved_ref_forGO.txt..fold.Enrichment., observed= goall$evolved_ref_forGO.txt..166., expected = goall$evolved_ref_forGO.txt..expected., pval = goall$evolved_ref_forGO.txt..raw.P.value., fdr=goall$evolved_ref_forGO.txt..FDR.)
dplyr::filter(evolvedtable, pval<0.01, enrich>1)
```

```
##                                                            bio enrich observed
## 1       positive regulation of meiotic cell cycle (GO:0051446)  30.94        2
## 2 positive regulation of meiotic nuclear division (GO:0045836)  51.57        2
## 3         positive regulation of nuclear division (GO:0051785)  15.47        2
## 4                         protein phosphorylation (GO:0006468)   2.14       16
## 5                                 phosphorylation (GO:0016310)   2.02       18
## 6          regulation of meiotic nuclear division (GO:0040020)  15.47        2
## 7       cellular response to gibberellin stimulus (GO:0071370)    9.1        4
##   expected    pval fdr
## 1     0.06 0.00326   1
## 2     0.04 0.00158   1
## 3     0.13 0.00982   1
## 4     7.47 0.00426   1
## 5     8.89 0.00489   1
## 6     0.13 0.00982   1
## 7     0.44 0.00135   1
```

```r
##adaptive plasticity
adaptive_plasticitytable = data.frame(bio= goall$GO.biological.process.complete, enrich = goall$adaptive_plasticity_ref_forGO.txt..fold.Enrichment., observed= goall$adaptive_plasticity_ref_forGO.txt..41., expected = goall$adaptive_plasticity_ref_forGO.txt..expected., pval = goall$adaptive_plasticity_ref_forGO.txt..raw.P.value., fdr=goall$adaptive_plasticity_ref_forGO.txt..FDR.)
dplyr::filter(adaptive_plasticitytable, pval < 0.01, enrich>1)
```

```
##                                                bio enrich observed expected
## 1             protein phosphorylation (GO:0006468)   3.25        6     1.85
## 2                     phosphorylation (GO:0016310)   3.19        7     2.20
## 3           pollen-pistil interaction (GO:0009875)  15.28        2     0.13
## 4               recognition of pollen (GO:0048544)  18.98        2     0.11
## 5                    cell recognition (GO:0008037)   17.9        2     0.11
## 6  regulation of root meristem growth (GO:0010082)  32.97        2     0.06
## 7                meristem maintenance (GO:0010073)  10.33        3     0.29
## 8       regulation of meristem growth (GO:0010075)  16.06        2     0.12
## 9                         sporulation (GO:0043934)  22.37        2     0.09
## 10                 sexual sporulation (GO:0034293)  22.37        2     0.09
## 11            plant-type sporogenesis (GO:0048236)  22.37        2     0.09
##       pval fdr
## 1  0.00984   1
## 2  0.00583   1
## 3  0.00822   1
## 4  0.00550   1
## 5  0.00613   1
## 6  0.00200   1
## 7  0.00328   1
## 8  0.00749   1
## 9  0.00406   1
## 10 0.00406   1
## 11 0.00406   1
```

```r
##maladaptive plasticity
maladaptive_plasticitytable = data.frame(bio= goall$GO.biological.process.complete, enrich = goall$maladaptive_plasticity_ref_forGO.txt..fold.Enrichment., observed= goall$maladaptive_plasticity_ref_forGO.txt..15., expected = goall$maladaptive_plasticity_ref_forGO.txt..expected., pval = goall$maladaptive_plasticity_ref_forGO.txt..raw.P.value., fdr=goall$maladaptive_plasticity_ref_forGO.txt..FDR.)
dplyr::filter(maladaptive_plasticitytable, pval < 0.01, enrich>1)
```

```
##                                                            bio enrich observed
## 1                         protein phosphorylation (GO:0006468)   5.92        4
## 2                                 phosphorylation (GO:0016310)   4.98        4
## 3                    phosphorus metabolic process (GO:0006793)   4.52        6
## 4 phosphate-containing compound metabolic process (GO:0006796)   4.63        6
## 5                    protein modification process (GO:0036211)   3.45        6
## 6           cellular protein modification process (GO:0006464)   3.45        6
##   expected    pval fdr
## 1     0.68 0.00381   1
## 2     0.80 0.00705   1
## 3     1.33 0.00121   1
## 4     1.30 0.00106   1
## 5     1.74 0.00480   1
## 6     1.74 0.00480   1
```

```r
###looking at up and down seperately 
dir1 = read.csv('data/GOlists/plasticresults1.csv')
dir2 = read.csv('data/GOlists/plasticresults2.csv')

##plastic up
plasticuptable = data.frame(bio= dir1$GO.biological.process.complete, enrich = dir1$plasticup.txt..fold.Enrichment., observed= dir1$plasticup.txt..4158., expected = dir1$plasticup.txt..expected., pval = dir1$plasticup.txt..raw.P.value., fdr=dir1$plasticup.txt..FDR.)
#dplyr::filter(plasticuptable, fdr < 0.05, enrich>1)

dplyr::filter(plasticuptable, fdr < 1e-3, enrich>1)
```

```
##                                                        bio enrich observed
## 1  interspecies interaction between organisms (GO:0044419)   1.52      316
## 2                        response to stimulus (GO:0050896)   1.17     1201
## 3        response to external biotic stimulus (GO:0043207)   1.51      309
## 4                  response to other organism (GO:0051707)   1.51      309
## 5                 response to biotic stimulus (GO:0009607)   1.51      310
## 6          defense response to other organism (GO:0098542)   1.56      237
## 7               response to external stimulus (GO:0009605)   1.39      399
## 8                            defense response (GO:0006952)   1.49      250
## 9                          response to fungus (GO:0009620)   1.84      118
## 10     response to oxygen-containing compound (GO:1901700)   1.36      409
## 11                       response to chemical (GO:0042221)    1.3      674
## 12                 defense response to fungus (GO:0050832)    1.9       89
## 13               protein modification process (GO:0036211)   1.25      601
## 14      cellular protein modification process (GO:0006464)   1.25      601
## 15              response to organic substance (GO:0010033)   1.35      439
## 16                         response to stress (GO:0006950)   1.23      742
## 17     cellular response to chemical stimulus (GO:0070887)   1.41      298
##    expected     pval      fdr
## 1    207.88 8.66e-10 1.65e-06
## 2   1029.06 1.56e-07 6.84e-05
## 3    204.32 2.08e-09 1.98e-06
## 4    204.32 2.08e-09 1.70e-06
## 5    204.65 1.70e-09 1.94e-06
## 6    151.54 2.51e-08 1.79e-05
## 7    287.22 2.93e-08 1.86e-05
## 8    167.73 2.32e-07 9.46e-05
## 9     64.11 2.93e-07 1.11e-04
## 10   301.14 1.33e-07 6.31e-05
## 11   518.41 1.58e-09 2.26e-06
## 12    46.95 3.58e-06 9.29e-04
## 13   481.50 1.47e-06 4.19e-04
## 14   481.50 1.47e-06 3.99e-04
## 15   326.07 1.05e-07 5.43e-05
## 16   604.55 4.33e-07 1.55e-04
## 17   210.80 6.18e-07 2.07e-04
```

```r
##plastic down
plasticdowntable = data.frame(bio= dir1$GO.biological.process.complete, enrich = dir1$plasticdown.txt..fold.Enrichment., observed= dir1$plasticdown.txt..4129., expected = dir1$plasticdown.txt..expected., pval = dir1$plasticdown.txt..raw.P.value., fdr=dir1$plasticdown.txt..FDR.)
#dplyr::filter(plasticdowntable, fdr < 0.05, enrich>1)

dplyr::filter(plasticdowntable, fdr < 1e-3, enrich>1)
```

```
##                                           bio enrich observed expected     pval
## 1 photosynthesis, light reaction (GO:0019684)   2.55       78    30.55 2.07e-09
## 2                 photosynthesis (GO:0015979)   2.51      112    44.70 1.90e-12
## 3           plastid organization (GO:0009657)   2.06      153    74.28 1.37e-11
## 4       chloroplast organization (GO:0009658)   2.01      115    57.24 1.24e-08
##        fdr
## 1 3.93e-06
## 2 1.09e-08
## 3 3.92e-08
## 4 1.78e-05
```

```r
##plastic 8 down
plastic8downtable = data.frame(bio= dir1$GO.biological.process.complete, enrich = dir1$plastic8down.txt..fold.Enrichment., observed= dir1$plastic8down.txt..2285., expected = dir1$plastic8down.txt..expected., pval = dir1$plastic8down.txt..raw.P.value., fdr=dir1$plastic8down.txt..FDR.)
dplyr::filter(plastic8downtable, fdr < 0.05, enrich>1)
```

```
##                                                                  bio enrich
## 1            protoporphyrinogen IX biosynthetic process (GO:0006782)   5.06
## 2                                   plastid translation (GO:0032544)   4.76
## 3               protoporphyrinogen IX metabolic process (GO:0046501)   5.06
## 4     photosynthesis, light harvesting in photosystem I (GO:0009768)   4.57
## 5                         photosynthesis, dark reaction (GO:0019685)    4.5
## 6                     reductive pentose-phosphate cycle (GO:0019253)   4.42
## 7                           chloroplast rRNA processing (GO:1901259)   5.19
## 8                                 photosystem II repair (GO:0010206)   5.06
## 9                               photosystem II assembly (GO:0010207)   4.21
## 10                                      carbon fixation (GO:0015977)   4.12
## 11                     photosynthesis, light harvesting (GO:0009765)   4.46
## 12                               photosystem I assembly (GO:0048564)   4.76
## 13                     chlorophyll biosynthetic process (GO:0015995)   4.35
## 14                       photosynthesis, light reaction (GO:0019684)   3.96
## 15                      thylakoid membrane organization (GO:0010027)   3.53
## 16                                       photosynthesis (GO:0015979)   3.84
## 17                        plastid membrane organization (GO:0009668)   3.45
## 18         regulation of photosynthesis, light reaction (GO:0042548)   3.75
## 19                     protein targeting to chloroplast (GO:0045036)   3.05
## 20 establishment of protein localization to chloroplast (GO:0072596)   3.05
## 21                                       protein repair (GO:0030091)   3.97
## 22                      carotenoid biosynthetic process (GO:0016117)   3.48
## 23                  tetraterpenoid biosynthetic process (GO:0016109)   3.48
## 24                  protein localization to chloroplast (GO:0072598)   2.89
## 25                    tetrapyrrole biosynthetic process (GO:0033014)    3.7
## 26                         carotenoid metabolic process (GO:0016116)   3.12
## 27                     tetraterpenoid metabolic process (GO:0016108)   3.12
## 28   porphyrin-containing compound biosynthetic process (GO:0006779)   3.55
## 29              photosynthetic electron transport chain (GO:0009767)   3.14
## 30                             glucan catabolic process (GO:0009251)   3.65
## 31                        chlorophyll metabolic process (GO:0015994)   3.35
## 32                                 plastid organization (GO:0009657)   2.58
## 33                       tetrapyrrole metabolic process (GO:0033013)   3.06
## 34                             chloroplast organization (GO:0009658)   2.49
## 35      porphyrin-containing compound metabolic process (GO:0006778)   2.96
## 36                         pigment biosynthetic process (GO:0046148)   2.63
## 37                            glucose metabolic process (GO:0006006)   2.99
## 38                          ketone biosynthetic process (GO:0042181)   2.89
## 39                            pigment metabolic process (GO:0042440)   2.47
## 40                     monosaccharide metabolic process (GO:0005996)   2.23
## 41                               cell redox homeostasis (GO:0045454)   2.54
## 42       generation of precursor metabolites and energy (GO:0006091)   2.42
## 43                     response to high light intensity (GO:0009644)   2.81
## 44                       carbohydrate metabolic process (GO:0005975)   1.44
## 45                           lipid biosynthetic process (GO:0008610)   1.67
## 46                  small molecule biosynthetic process (GO:0044283)   1.83
## 47                          response to light intensity (GO:0009642)   2.05
## 48                       proton transmembrane transport (GO:1902600)   2.42
## 49                    organic acid biosynthetic process (GO:0016053)   1.74
## 50                 carboxylic acid biosynthetic process (GO:0046394)   1.74
## 51                           response to light stimulus (GO:0009416)   1.47
## 52                     heterocycle biosynthetic process (GO:0018130)   1.44
## 53                     small molecule metabolic process (GO:0044281)   1.65
## 54                          oxidation-reduction process (GO:0055114)   1.61
## 55                                response to radiation (GO:0009314)   1.41
## 56         organic cyclic compound biosynthetic process (GO:1901362)   1.39
## 57                            oxoacid metabolic process (GO:0043436)   1.61
## 58                       organic acid metabolic process (GO:0006082)    1.6
## 59                    carboxylic acid metabolic process (GO:0019752)   1.61
## 60                          ion transmembrane transport (GO:0034220)   1.85
## 61                cellular amino acid metabolic process (GO:0006520)   1.65
## 62                                        ion transport (GO:0006811)   1.53
## 63                                 biosynthetic process (GO:0009058)   1.35
## 64               organic substance biosynthetic process (GO:1901576)   1.36
## 65                              transmembrane transport (GO:0055085)   1.52
## 66                        cellular biosynthetic process (GO:0044249)   1.36
## 67      cellular nitrogen compound biosynthetic process (GO:0044271)   1.42
## 68         organonitrogen compound biosynthetic process (GO:1901566)   1.49
## 69                           amide biosynthetic process (GO:0043604)   1.53
## 70                     cellular amide metabolic process (GO:0043603)   1.45
##    observed expected     pval      fdr
## 1         9     1.78 8.81e-04 3.96e-02
## 2        11     2.31 3.37e-04 1.96e-02
## 3         9     1.78 8.81e-04 3.93e-02
## 4        13     2.85 1.29e-04 8.98e-03
## 5        12     2.67 2.59e-04 1.57e-02
## 6        11     2.49 5.20e-04 2.65e-02
## 7        12     2.31 1.04e-04 7.42e-03
## 8         9     1.78 8.81e-04 3.90e-02
## 9        18     4.27 1.44e-05 1.37e-03
## 10       11     2.67 7.79e-04 3.74e-02
## 11       23     5.16 4.83e-07 6.89e-05
## 12       11     2.31 3.37e-04 1.94e-02
## 13       24     5.52 3.72e-07 5.44e-05
## 14       67    16.90 4.79e-16 1.37e-12
## 15       27     7.65 1.40e-06 1.90e-04
## 16       95    24.73 1.65e-21 9.41e-18
## 17       27     7.83 1.93e-06 2.51e-04
## 18       12     3.20 8.29e-04 3.88e-02
## 19       19     6.23 3.52e-04 1.97e-02
## 20       19     6.23 3.52e-04 1.95e-02
## 21       12     3.03 5.76e-04 2.89e-02
## 22       13     3.74 8.43e-04 3.91e-02
## 23       13     3.74 8.43e-04 3.88e-02
## 24       19     6.58 4.61e-04 2.42e-02
## 25       27     7.30 7.11e-07 9.90e-05
## 26       15     4.80 7.98e-04 3.80e-02
## 27       15     4.80 7.98e-04 3.76e-02
## 28       24     6.76 4.94e-06 5.88e-04
## 29       19     6.05 1.58e-04 1.06e-02
## 30       13     3.56 6.03e-04 3.00e-02
## 31       28     8.36 1.93e-06 2.56e-04
## 32      106    41.11 1.97e-14 3.76e-11
## 33       31    10.14 2.88e-06 3.66e-04
## 34       79    31.67 1.69e-10 8.04e-08
## 35       29     9.79 8.81e-06 8.82e-04
## 36       36    13.70 7.48e-06 7.76e-04
## 37       17     5.69 8.61e-04 3.93e-02
## 38       17     5.87 9.76e-04 4.22e-02
## 39       40    16.19 6.63e-06 7.28e-04
## 40       31    13.88 3.76e-04 2.06e-02
## 41       24     9.43 3.29e-04 1.94e-02
## 42      105    43.42 6.75e-13 6.42e-10
## 43       23     8.19 1.34e-04 9.12e-03
## 44      127    88.44 3.12e-04 1.86e-02
## 45       94    56.41 2.50e-05 2.23e-03
## 46      141    76.87 1.36e-09 3.70e-07
## 47       38    18.51 3.39e-04 1.92e-02
## 48       22     9.08 9.16e-04 4.02e-02
## 49       93    53.56 6.32e-06 7.22e-04
## 50       93    53.56 6.32e-06 7.08e-04
## 51      127    86.13 1.34e-04 9.21e-03
## 52      123    85.24 3.80e-04 2.06e-02
## 53      282   170.47 9.43e-14 1.35e-10
## 54      122    75.80 6.19e-06 7.22e-04
## 55      127    90.04 6.47e-04 3.16e-02
## 56      142   101.96 4.90e-04 2.52e-02
## 57      179   111.39 3.98e-08 7.10e-06
## 58      179   111.93 4.31e-08 7.47e-06
## 59      171   106.06 6.28e-08 1.02e-05
## 60       67    36.12 2.35e-05 2.13e-03
## 61       75    45.38 2.37e-04 1.50e-02
## 62      103    67.44 2.15e-04 1.40e-02
## 63      425   314.61 4.12e-09 9.40e-07
## 64      413   303.75 3.69e-09 8.77e-07
## 65      133    87.37 2.25e-05 2.07e-03
## 66      396   291.12 9.52e-09 2.01e-06
## 67      196   138.09 1.28e-05 1.23e-03
## 68      226   151.25 8.14e-08 1.29e-05
## 69       98    64.06 2.56e-04 1.57e-02
## 70      113    78.12 5.42e-04 2.74e-02
```

```r
##plastic 8 up
plastic8uptable = data.frame(bio= dir2$GO.biological.process.complete, enrich = dir2$plastic8up.txt..fold.Enrichment., observed= dir2$plastic8up.txt..2465., expected = dir2$plastic8up.txt..expected., pval = dir2$plastic8up.txt..raw.P.value., fdr=dir2$plastic8up.txt..FDR.)
dplyr::filter(plastic8uptable, fdr < 0.05, enrich>1)
```

```
##                                                                            bio
## 1                         positive regulation of defense response (GO:0031349)
## 2              positive regulation of response to biotic stimulus (GO:0002833)
## 3                                              response to chitin (GO:0010200)
## 4                       secondary metabolite biosynthetic process (GO:0044550)
## 5                                     secondary metabolic process (GO:0019748)
## 6                                              response to fungus (GO:0009620)
## 7                                      defense response to fungus (GO:0050832)
## 8                                       response to oxygen levels (GO:0070482)
## 9                                             response to hypoxia (GO:0001666)
## 10                            response to decreased oxygen levels (GO:0036293)
## 11                    cellular response to abscisic acid stimulus (GO:0071215)
## 12                                   cellular response to alcohol (GO:0097306)
## 13                             cellular response to oxygen levels (GO:0071453)
## 14                   cellular response to decreased oxygen levels (GO:0036294)
## 15                                 regulation of defense response (GO:0031347)
## 16                                   cellular response to hypoxia (GO:0071456)
## 17                      regulation of response to biotic stimulus (GO:0002831)
## 18                             defense response to other organism (GO:0098542)
## 19                           regulation of innate immune response (GO:0045088)
## 20                    regulation of response to external stimulus (GO:0032101)
## 21                                           response to wounding (GO:0009611)
## 22                     defense response, incompatible interaction (GO:0009814)
## 23                                                immune response (GO:0006955)
## 24                            response to organonitrogen compound (GO:0010243)
## 25                     interspecies interaction between organisms (GO:0044419)
## 26                            regulation of immune system process (GO:0002682)
## 27                                         innate immune response (GO:0045087)
## 28                                    response to biotic stimulus (GO:0009607)
## 29                           response to external biotic stimulus (GO:0043207)
## 30                                     response to other organism (GO:0051707)
## 31                                          immune system process (GO:0002376)
## 32                                  defense response to bacterium (GO:0042742)
## 33                                          response to bacterium (GO:0009617)
## 34                                  regulation of immune response (GO:0050776)
## 35                                           response to ethylene (GO:0009723)
## 36                                               defense response (GO:0006952)
## 37                                    protein autophosphorylation (GO:0046777)
## 38                                                leaf senescence (GO:0010150)
## 39                    positive regulation of response to stimulus (GO:0048584)
## 40                                      response to abscisic acid (GO:0009737)
## 41                                                          aging (GO:0007568)
## 42                                            response to alcohol (GO:0097305)
## 43                                      response to acid chemical (GO:0001101)
## 44                                  response to external stimulus (GO:0009605)
## 45                                  response to nitrogen compound (GO:1901698)
## 46                               regulation of response to stress (GO:0080134)
## 47                         response to oxygen-containing compound (GO:1901700)
## 48                                        protein phosphorylation (GO:0006468)
## 49                                              response to lipid (GO:0033993)
## 50                         cellular response to chemical stimulus (GO:0070887)
## 51                                         protein ubiquitination (GO:0016567)
## 52                                   response to abiotic stimulus (GO:0009628)
## 53                cellular response to oxygen-containing compound (GO:1901701)
## 54              protein modification by small protein conjugation (GO:0032446)
## 55                                           response to chemical (GO:0042221)
## 56                                                phosphorylation (GO:0016310)
## 57                           ethylene-activated signaling pathway (GO:0009873)
## 58                                  response to organic substance (GO:0010033)
## 59                                response to endogenous stimulus (GO:0009719)
## 60                                     cellular response to lipid (GO:0071396)
## 61                                            response to hormone (GO:0009725)
## 62                                             cell communication (GO:0007154)
## 63                         cellular response to organic substance (GO:0071310)
## 64                                                      signaling (GO:0023052)
## 65                                            signal transduction (GO:0007165)
## 66                                   phosphorus metabolic process (GO:0006793)
## 67                phosphate-containing compound metabolic process (GO:0006796)
## 68                                           response to stimulus (GO:0050896)
## 69                       cellular response to endogenous stimulus (GO:0071495)
## 70                          cellular response to hormone stimulus (GO:0032870)
## 71                                   protein modification process (GO:0036211)
## 72                          cellular protein modification process (GO:0006464)
## 73                                             response to stress (GO:0006950)
## 74                             cellular protein metabolic process (GO:0044267)
## 75                             hormone-mediated signaling pathway (GO:0009755)
## 76                                      protein metabolic process (GO:0019538)
## 77                                  cellular response to stimulus (GO:0051716)
## 78                                     macromolecule modification (GO:0043412)
## 79                                 regulation of cellular process (GO:0050794)
## 80                               regulation of biological process (GO:0050789)
## 81                                          biological regulation (GO:0065007)
## 82                             regulation of response to stimulus (GO:0048583)
## 83                        regulation of primary metabolic process (GO:0080090)
## 84                       regulation of cellular metabolic process (GO:0031323)
## 85                                regulation of metabolic process (GO:0019222)
## 86              regulation of nitrogen compound metabolic process (GO:0051171)
## 87                             regulation of biosynthetic process (GO:0009889)
## 88                         regulation of RNA biosynthetic process (GO:2001141)
## 89                     regulation of transcription, DNA-templated (GO:0006355)
## 90             regulation of nucleic acid-templated transcription (GO:1903506)
## 91                    regulation of cellular biosynthetic process (GO:0031326)
## 92      regulation of cellular macromolecule biosynthetic process (GO:2000112)
## 93               regulation of macromolecule biosynthetic process (GO:0010556)
## 94                            regulation of RNA metabolic process (GO:0051252)
## 95                  regulation of macromolecule metabolic process (GO:0060255)
## 96 regulation of nucleobase-containing compound metabolic process (GO:0019219)
## 97                                  regulation of gene expression (GO:0010468)
##    enrich observed expected     pval      fdr
## 1    2.55       22     8.64 6.52e-04 3.26e-02
## 2    2.67       20     7.49 9.91e-04 4.49e-02
## 3    3.47       50    14.40 1.19e-10 7.52e-08
## 4    2.25       29    12.86 6.57e-04 3.26e-02
## 5    2.13       49    23.04 2.06e-05 1.68e-03
## 6    2.18       83    38.01 1.49e-08 3.16e-06
## 7    2.23       62    27.83 4.97e-07 6.60e-05
## 8     2.4       71    29.56 7.19e-09 2.16e-06
## 9    2.41       69    28.60 1.26e-08 2.88e-06
## 10   2.36       69    29.18 1.74e-08 3.56e-06
## 11   2.08       34    16.32 6.04e-04 3.08e-02
## 12   2.08       34    16.32 6.04e-04 3.05e-02
## 13   2.55       66    25.92 5.06e-09 1.81e-06
## 14   2.53       65    25.72 7.79e-09 2.22e-06
## 15    2.1       58    27.64 6.86e-06 6.75e-04
## 16   2.57       65    25.34 3.80e-09 1.55e-06
## 17   2.26       43    19.00 2.58e-05 2.02e-03
## 18   1.83      164    89.84 1.02e-10 7.29e-08
## 19   2.38       26    10.94 4.78e-04 2.53e-02
## 20   2.23       44    19.77 2.31e-05 1.83e-03
## 21   1.99       53    26.68 4.40e-05 3.11e-03
## 22   2.16       46    21.31 2.72e-05 2.10e-03
## 23   2.17       79    36.47 3.46e-08 6.38e-06
## 24   2.42       72    29.75 4.65e-09 1.77e-06
## 25   1.74      214   123.24 5.58e-12 6.37e-09
## 26   2.24       31    13.82 3.71e-04 2.02e-02
## 27   2.14       76    35.51 9.44e-08 1.63e-05
## 28   1.75      212   121.32 5.47e-12 7.81e-09
## 29   1.75      212   121.13 3.88e-12 1.11e-08
## 30   1.75      212   121.13 3.88e-12 7.38e-09
## 31   2.01       82    40.89 3.35e-07 4.78e-05
## 32   1.73       81    46.84 3.44e-05 2.58e-03
## 33   1.69       97    57.40 1.54e-05 1.29e-03
## 34   2.32       28    12.09 4.66e-04 2.51e-02
## 35    2.1       37    17.66 2.69e-04 1.52e-02
## 36   1.74      173    99.44 6.31e-10 3.27e-07
## 37   2.09       51    24.38 2.08e-05 1.68e-03
## 38   2.18       31    14.21 4.70e-04 2.51e-02
## 39   1.77       48    27.07 9.41e-04 4.41e-02
## 40    1.6       97    60.47 7.49e-05 4.92e-03
## 41      2       36    18.04 7.12e-04 3.48e-02
## 42    1.6       98    61.24 8.21e-05 5.26e-03
## 43   1.59       72    45.30 9.46e-04 4.39e-02
## 44    1.5      255   170.27 1.12e-08 2.67e-06
## 45   2.14       76    35.51 9.44e-08 1.59e-05
## 46   1.74       79    45.50 3.77e-05 2.76e-03
## 47   1.51      269   178.53 2.58e-09 1.13e-06
## 48   1.83      203   110.95 3.99e-13 2.28e-09
## 49   1.52      122    80.05 6.01e-05 3.99e-03
## 50    1.7      213   124.97 2.65e-11 2.16e-08
## 51   1.45      105    72.37 8.33e-04 3.93e-02
## 52   1.27      328   258.00 4.98e-05 3.42e-03
## 53   1.66      102    61.62 1.59e-05 1.32e-03
## 54   1.47      112    76.02 3.57e-04 1.96e-02
## 55    1.4      429   307.33 1.52e-10 8.69e-08
## 56   1.69      223   132.07 1.69e-11 1.61e-08
## 57   2.76       18     6.53 9.53e-04 4.39e-02
## 58   1.46      283   193.31 1.07e-08 2.79e-06
## 59   1.38      197   142.24 4.37e-05 3.12e-03
## 60   1.77       51    28.79 6.67e-04 3.28e-02
## 61   1.36      188   138.02 1.57e-04 9.34e-03
## 62   1.47      225   153.19 3.77e-07 5.12e-05
## 63   1.68      128    76.02 6.05e-07 7.85e-05
## 64   1.55      192   124.20 1.60e-07 2.48e-05
## 65   1.55      191   123.05 1.43e-07 2.27e-05
## 66   1.37      298   218.26 1.08e-06 1.34e-04
## 67   1.39      296   212.89 2.59e-07 3.80e-05
## 68   1.22      744   610.06 2.52e-08 4.80e-06
## 69    1.7      102    59.89 6.65e-06 6.66e-04
## 70   1.68       93    55.48 2.84e-05 2.17e-03
## 71   1.37      392   285.45 6.27e-09 2.11e-06
## 72   1.37      392   285.45 6.27e-09 1.99e-06
## 73    1.3      466   358.40 6.97e-08 1.24e-05
## 74    1.2      457   381.43 1.74e-04 1.03e-02
## 75   1.68       81    48.18 7.94e-05 5.15e-03
## 76   1.17      469   401.01 9.58e-04 4.38e-02
## 77   1.34      369   275.85 2.12e-07 3.18e-05
## 78    1.3      431   332.67 3.63e-07 5.06e-05
## 79   1.21      580   477.60 3.40e-06 3.66e-04
## 80   1.21      665   550.36 7.95e-07 1.01e-04
## 81   1.17      718   612.55 1.01e-05 9.34e-04
## 82   1.47      129    87.92 1.50e-04 9.13e-03
## 83   1.24      343   277.58 2.47e-04 1.41e-02
## 84   1.21      359   297.54 7.69e-04 3.72e-02
## 85   1.21      408   338.24 2.87e-04 1.60e-02
## 86   1.26      338   267.21 5.52e-05 3.71e-03
## 87   1.28      309   240.53 4.64e-05 3.23e-03
## 88   1.36      275   201.75 3.16e-06 3.47e-04
## 89   1.37      275   201.37 3.09e-06 3.53e-04
## 90   1.37      275   201.37 3.09e-06 3.46e-04
## 91   1.29      304   235.73 4.01e-05 2.90e-03
## 92   1.33      294   221.33 8.89e-06 8.32e-04
## 93   1.32      294   222.68 1.32e-05 1.17e-03
## 94   1.35      291   215.19 2.83e-06 3.30e-04
## 95   1.25      378   302.92 5.01e-05 3.41e-03
## 96   1.33      299   224.79 6.12e-06 6.24e-04
## 97   1.29      338   261.84 1.40e-05 1.23e-03
```

```r
##plastic 32 up
plastic32uptable = data.frame(bio= dir2$GO.biological.process.complete, enrich = dir2$plastic32up.txt..fold.Enrichment., observed= dir2$plastic32up.txt..4014., expected = dir2$plastic32up.txt..expected., pval = dir2$plastic32up.txt..raw.P.value., fdr=dir2$plastic32up.txt..FDR.)
dplyr::filter(plastic32uptable, fdr < 0.05, enrich>1)
```

```
##                                                                        bio
## 1                       monocarboxylic acid catabolic process (GO:0072329)
## 2                                 secondary metabolic process (GO:0019748)
## 3                                          response to fungus (GO:0009620)
## 4                                  defense response to fungus (GO:0050832)
## 5                                   response to oxygen levels (GO:0070482)
## 6                                         response to hypoxia (GO:0001666)
## 7                         response to decreased oxygen levels (GO:0036293)
## 8                              organic acid catabolic process (GO:0016054)
## 9                           carboxylic acid catabolic process (GO:0046395)
## 10                         defense response to other organism (GO:0098542)
## 11                        response to organonitrogen compound (GO:0010243)
## 12                 interspecies interaction between organisms (GO:0044419)
## 13                                        aerobic respiration (GO:0009060)
## 14                                response to biotic stimulus (GO:0009607)
## 15                       response to external biotic stimulus (GO:0043207)
## 16                                 response to other organism (GO:0051707)
## 17                           small molecule catabolic process (GO:0044282)
## 18                                           defense response (GO:0006952)
## 19                                  response to acid chemical (GO:0001101)
## 20                              response to external stimulus (GO:0009605)
## 21                              response to nitrogen compound (GO:1901698)
## 22                                      response to metal ion (GO:0010038)
## 23                                    response to cadmium ion (GO:0046686)
## 24                     response to oxygen-containing compound (GO:1901700)
## 25                     cellular response to chemical stimulus (GO:0070887)
## 26                            response to inorganic substance (GO:0010035)
## 27                         cellular protein catabolic process (GO:0044257)
## 28                                       response to chemical (GO:0042221)
## 29                                  protein catabolic process (GO:0030163)
## 30     modification-dependent macromolecule catabolic process (GO:0043632)
## 31                  organonitrogen compound catabolic process (GO:1901565)
## 32           modification-dependent protein catabolic process (GO:0019941)
## 33        energy derivation by oxidation of organic compounds (GO:0015980)
## 34                                 response to osmotic stress (GO:0006970)
## 35                                 cellular catabolic process (GO:0044248)
## 36                              response to organic substance (GO:0010033)
## 37 proteolysis involved in cellular protein catabolic process (GO:0051603)
## 38                                          catabolic process (GO:0009056)
## 39                        organic substance catabolic process (GO:1901575)
## 40              ubiquitin-dependent protein catabolic process (GO:0006511)
## 41                                         cell communication (GO:0007154)
## 42                     cellular response to organic substance (GO:0071310)
## 43                                                  signaling (GO:0023052)
## 44                                        signal transduction (GO:0007165)
## 45                               phosphorus metabolic process (GO:0006793)
## 46            phosphate-containing compound metabolic process (GO:0006796)
## 47                                       response to stimulus (GO:0050896)
## 48                  organonitrogen compound metabolic process (GO:1901564)
## 49                                       cellular respiration (GO:0045333)
## 50                                         response to stress (GO:0006950)
## 51                              cellular response to stimulus (GO:0051716)
##    enrich observed expected     pval      fdr
## 1    2.38       38    15.94 8.42e-05 7.51e-03
## 2    1.73       65    37.51 4.84e-04 3.07e-02
## 3    1.76      109    61.89 3.64e-06 5.78e-04
## 4    1.74       79    45.33 1.00e-04 8.55e-03
## 5    1.68       81    48.14 2.01e-04 1.47e-02
## 6     1.7       79    46.58 2.15e-04 1.56e-02
## 7    1.66       79    47.51 3.27e-04 2.33e-02
## 8    2.02       75    37.20 3.72e-06 5.74e-04
## 9    2.02       75    37.20 3.72e-06 5.59e-04
## 10   1.44      210   146.29 1.45e-05 1.88e-03
## 11   1.65       80    48.45 3.65e-04 2.51e-02
## 12   1.46      293   200.68 7.12e-08 3.70e-05
## 13   2.27       34    15.00 3.83e-04 2.60e-02
## 14   1.45      286   197.56 2.05e-07 8.35e-05
## 15   1.44      285   197.25 2.47e-07 8.29e-05
## 16   1.44      285   197.25 2.47e-07 7.83e-05
## 17   1.83       96    52.52 3.88e-06 5.68e-04
## 18   1.38      224   161.92 4.85e-05 4.86e-03
## 19   1.52      112    73.77 3.52e-04 2.45e-02
## 20   1.35      374   277.27 7.18e-07 1.71e-04
## 21   1.56       90    57.83 7.57e-04 4.41e-02
## 22    1.5      134    89.09 1.20e-04 9.69e-03
## 23    1.5      103    68.46 7.34e-04 4.41e-02
## 24   1.34      391   290.71 5.13e-07 1.39e-04
## 25    1.4      285   203.50 1.70e-06 3.35e-04
## 26   1.43      255   178.49 2.09e-06 3.73e-04
## 27   1.47      145    98.47 1.53e-04 1.18e-02
## 28   1.34      669   500.46 2.75e-11 7.85e-08
## 29   1.45      148   101.91 1.93e-04 1.45e-02
## 30   1.54      137    89.09 4.76e-05 4.85e-03
## 31   1.53      215   140.67 3.74e-07 1.12e-04
## 32   1.55      134    86.28 3.73e-05 4.02e-03
## 33   1.85       51    27.51 6.22e-04 3.82e-02
## 34   1.42      146   102.84 4.60e-04 2.95e-02
## 35   1.46      324   221.32 1.11e-08 7.94e-06
## 36   1.32      414   314.78 1.42e-06 3.24e-04
## 37    1.5      144    96.28 8.46e-05 7.32e-03
## 38   1.41      366   259.14 2.41e-08 1.53e-05
## 39   1.41      335   238.20 1.48e-07 6.52e-05
## 40   1.53      130    84.71 8.44e-05 7.41e-03
## 41   1.32      329   249.45 1.88e-05 2.39e-03
## 42   1.38      171   123.79 4.37e-04 2.83e-02
## 43   1.31      265   202.25 1.97e-04 1.46e-02
## 44   1.31      263   200.37 1.86e-04 1.41e-02
## 45   1.21      429   355.42 6.07e-04 3.77e-02
## 46   1.21      418   346.67 7.45e-04 4.43e-02
## 47   1.17     1159   993.42 2.30e-07 8.77e-05
## 48   1.17      997   851.19 1.60e-06 3.52e-04
## 49   1.96       46    23.44 5.32e-04 3.34e-02
## 50   1.24      721   583.61 2.33e-07 8.32e-05
## 51   1.21      544   449.20 6.62e-05 6.52e-03
```

```r
##plastic 32 down
plastic32downtable = data.frame(bio= dir2$GO.biological.process.complete, enrich = dir2$plastic32down.txt..fold.Enrichment., observed= dir2$plastic32down.txt..3826., expected = dir2$plastic32down.txt..expected., pval = dir2$plastic32down.txt..raw.P.value., fdr=dir2$plastic32down.txt..FDR.)
dplyr::filter(plastic32downtable, fdr < 0.05, enrich>1)
```

```
##                                           bio enrich observed expected     pval
## 1                 photosynthesis (GO:0015979)   2.58      107    41.42 6.84e-13
## 2 photosynthesis, light reaction (GO:0019684)   2.61       74    28.31 1.77e-09
## 3           plastid organization (GO:0009657)   1.99      137    68.83 5.22e-10
## 4       chloroplast organization (GO:0009658)   1.94      103    53.04 1.95e-07
##        fdr
## 1 3.91e-09
## 2 3.38e-06
## 3 1.49e-06
## 4 2.78e-04
```

```r
eud = read.csv('data/GOlists/evolupdownGOresults.csv')
##evolved up
evolveduptable = data.frame(bio= eud$GO.biological.process.complete, enrich = eud$evolvedup.txt..fold.Enrichment., observed= eud$evolvedup.txt..84., expected = eud$evolvedup.txt..expected., pval = eud$evolvedup.txt..raw.P.value., fdr=eud$evolvedup.txt..FDR.)
dplyr::filter(evolveduptable, pval<0.01, enrich>1)
```

```
##                                                       bio enrich observed
## 1                    response to desiccation (GO:0009269)  23.52        2
## 2  cellular response to gibberellin stimulus (GO:0071370)  13.49        3
## 3                  response to abscisic acid (GO:0009737)    3.4        7
## 4                        response to alcohol (GO:0097305)   3.35        7
## 5                          response to lipid (GO:0033993)   2.93        8
## 6 interspecies interaction between organisms (GO:0044419)   2.62       11
## 7       response to external biotic stimulus (GO:0043207)   2.42       10
## 8                 response to other organism (GO:0051707)   2.42       10
## 9                response to biotic stimulus (GO:0009607)   2.42       10
##   expected    pval fdr
## 1     0.09 0.00415   1
## 2     0.22 0.00175   1
## 3     2.06 0.00491   1
## 4     2.09 0.00524   1
## 5     2.73 0.00630   1
## 6     4.20 0.00324   1
## 7     4.13 0.00846   1
## 8     4.13 0.00846   1
## 9     4.13 0.00855   1
```

```r
##evolved down
evolveddowntable = data.frame(bio= eud$GO.biological.process.complete, enrich = eud$evolveddown.txt..fold.Enrichment., observed= eud$evolveddown.txt..85., expected = eud$evolveddown.txt..expected., pval = eud$evolveddown.txt..raw.P.value., fdr=eud$evolveddown.txt..FDR.)
dplyr::filter(evolveddowntable, pval<0.01, enrich>1)
```

```
##                                                      bio enrich observed
## 1                   protein phosphorylation (GO:0006468)   2.61       10
## 2                           phosphorylation (GO:0016310)   2.63       12
## 3   positive regulation of nuclear division (GO:0051785)  30.21        2
## 4    regulation of meiotic nuclear division (GO:0040020)  30.21        2
## 5 positive regulation of meiotic cell cycle (GO:0051446)  60.43        2
## 6          regulation of meiotic cell cycle (GO:0051445)  23.24        2
## 7 positive regulation of cell cycle process (GO:0090068)   15.9        2
## 8               regulation of cell division (GO:0051302)  10.79        3
##   expected     pval fdr
## 1     3.83 0.005100   1
## 2     4.55 0.002010   1
## 3     0.07 0.002700   1
## 4     0.07 0.002700   1
## 5     0.03 0.000878   1
## 6     0.09 0.004240   1
## 7     0.13 0.008280   1
## 8     0.28 0.003190   1
```

Look at most common GO terms (not enrichment)

```r
dplyr::filter(evolveduptable, observed > 10)
```

```
##                                                             bio enrich observed
## 1       interspecies interaction between organisms (GO:0044419)   2.62       11
## 2                    response to external stimulus (GO:0009605)    1.9       11
## 3  phosphate-containing compound metabolic process (GO:0006796)   1.52       11
## 4                     phosphorus metabolic process (GO:0006793)   1.48       11
## 5                             response to chemical (GO:0042221)   1.43       15
## 6                     response to abiotic stimulus (GO:0009628)   1.36       12
## 7                             response to stimulus (GO:0050896)   1.25       26
## 8                     protein modification process (GO:0036211)   1.23       12
## 9            cellular protein modification process (GO:0006464)   1.23       12
## 10                   cellular response to stimulus (GO:0051716)   1.17       11
## 11                              response to stress (GO:0006950)   1.15       14
## 12                                  Unclassified (UNCLASSIFIED)   1.14       16
## 13                       protein metabolic process (GO:0019538)    1.1       15
## 14              cellular protein metabolic process (GO:0044267)   1.08       14
## 15                      macromolecule modification (GO:0043412)   1.06       12
## 16                multicellular organismal process (GO:0032501)   0.99       11
## 17                              biological_process (GO:0008150)   0.97       68
## 18                               metabolic process (GO:0008152)   0.96       33
## 19                                cellular process (GO:0009987)   0.94       44
## 20                      cellular metabolic process (GO:0044237)   0.93       29
## 21                  regulation of cellular process (GO:0050794)   0.92       15
## 22                           biological regulation (GO:0065007)   0.91       19
## 23             organic substance metabolic process (GO:0071704)   0.89       28
## 24                regulation of biological process (GO:0050789)   0.85       16
## 25        cellular macromolecule metabolic process (GO:0044260)   0.85       15
## 26                 macromolecule metabolic process (GO:0043170)   0.85       19
## 27       organonitrogen compound metabolic process (GO:1901564)   0.84       15
## 28                       primary metabolic process (GO:0044238)   0.83       24
## 29             nitrogen compound metabolic process (GO:0006807)   0.76       19
##    expected    pval fdr
## 1      4.20 0.00324   1
## 2      5.80 0.04730   1
## 3      7.25 0.16800   1
## 4      7.44 0.17600   1
## 5     10.47 0.13700   1
## 6      8.79 0.28000   1
## 7     20.79 0.20500   1
## 8      9.73 0.39600   1
## 9      9.73 0.39600   1
## 10     9.40 0.60100   1
## 11    12.21 0.53700   1
## 12    14.00 0.55700   1
## 13    13.67 0.65700   1
## 14    13.00 0.76200   1
## 15    11.34 0.75100   1
## 16    11.13 1.00000   1
## 17    70.00 0.55700   1
## 18    34.38 0.82400   1
## 19    46.76 0.58200   1
## 20    31.16 0.65200   1
## 21    16.28 0.89000   1
## 22    20.87 0.70500   1
## 23    31.49 0.49800   1
## 24    18.75 0.51400   1
## 25    17.65 0.59100   1
## 26    22.42 0.45900   1
## 27    17.81 0.50500   1
## 28    29.02 0.30000   1
## 29    25.11 0.15300   1
```

```r
dplyr::filter(evolveddowntable, observed > 10)
```

```
##                                                              bio enrich
## 1                                   phosphorylation (GO:0016310)   2.63
## 2   phosphate-containing compound metabolic process (GO:0006796)   1.63
## 3                      phosphorus metabolic process (GO:0006793)   1.86
## 4                              response to stimulus (GO:0050896)   0.95
## 5                      protein modification process (GO:0036211)   1.52
## 6             cellular protein modification process (GO:0006464)   1.52
## 7                                response to stress (GO:0006950)   1.05
## 8                                    Unclassified (UNCLASSIFIED)   0.85
## 9                         protein metabolic process (GO:0019538)    1.3
## 10               cellular protein metabolic process (GO:0044267)   1.37
## 11                       macromolecule modification (GO:0043412)   1.31
## 12                             reproductive process (GO:0022414)   1.68
## 13                 multicellular organismal process (GO:0032501)   1.42
## 14                                     reproduction (GO:0000003)   1.66
## 15                               biological_process (GO:0008150)   1.03
## 16                                metabolic process (GO:0008152)   1.15
## 17                                 cellular process (GO:0009987)   1.04
## 18                       cellular metabolic process (GO:0044237)   1.21
## 19                   regulation of cellular process (GO:0050794)   0.97
## 20                            biological regulation (GO:0065007)   0.99
## 21                               system development (GO:0048731)   1.41
## 22               multicellular organism development (GO:0007275)   1.27
## 23              organic substance metabolic process (GO:0071704)   1.16
## 24                 anatomical structure development (GO:0048856)   1.37
## 25                 regulation of biological process (GO:0050789)      1
## 26         cellular macromolecule metabolic process (GO:0044260)    1.4
## 27                  macromolecule metabolic process (GO:0043170)   1.23
## 28        organonitrogen compound metabolic process (GO:1901564)   1.22
## 29                            developmental process (GO:0032502)   1.32
## 30                        primary metabolic process (GO:0044238)   1.09
## 31              nitrogen compound metabolic process (GO:0006807)   1.22
## 32     cellular aromatic compound metabolic process (GO:0006725)   1.21
## 33        organic cyclic compound metabolic process (GO:1901360)   1.09
## 34   developmental process involved in reproduction (GO:0003006)   1.81
## 35                    cellular biosynthetic process (GO:0044249)   1.48
## 36           organic substance biosynthetic process (GO:1901576)   1.33
## 37 nucleobase-containing compound metabolic process (GO:0006139)   1.32
## 38                             biosynthetic process (GO:0009058)   1.45
## 39      cellular macromolecule biosynthetic process (GO:0034645)   2.02
## 40                    heterocycle metabolic process (GO:0046483)   1.26
## 41               macromolecule biosynthetic process (GO:0009059)   1.95
## 42     cellular nitrogen compound metabolic process (GO:0034641)   1.24
##    observed expected    pval fdr
## 1        12     4.55 0.00201   1
## 2        12     7.34 0.08060   1
## 3        14     7.53 0.02070   1
## 4        20    21.04 0.90000   1
## 5        15     9.84 0.08870   1
## 6        15     9.84 0.08870   1
## 7        13    12.36 0.87700   1
## 8        12    14.17 0.66100   1
## 9        18    13.83 0.23700   1
## 10       18    13.15 0.17400   1
## 11       15    11.47 0.26500   1
## 12       12     7.16 0.07480   1
## 13       16    11.26 0.14700   1
## 14       12     7.21 0.07620   1
## 15       73    70.83 0.66100   1
## 16       40    34.79 0.26900   1
## 17       49    47.32 0.74400   1
## 18       38    31.53 0.17600   1
## 19       16    16.47 1.00000   1
## 20       21    21.12 1.00000   1
## 21       11     7.78 0.25400   1
## 22       13    10.23 0.32000   1
## 23       37    31.87 0.26200   1
## 24       16    11.66 0.20300   1
## 25       19    18.98 1.00000   1
## 26       25    17.86 0.06220   1
## 27       28    22.68 0.21800   1
## 28       22    18.02 0.28800   1
## 29       16    12.15 0.21700   1
## 30       32    29.36 0.56800   1
## 31       31    25.41 0.19200   1
## 32       14    11.52 0.42600   1
## 33       13    11.96 0.75400   1
## 34       11     6.08 0.05380   1
## 35       16    10.83 0.10200   1
## 36       15    11.30 0.26000   1
## 37       13     9.87 0.30600   1
## 38       17    11.70 0.11300   1
## 39       11     5.43 0.02350   1
## 40       14    11.15 0.33600   1
## 41       11     5.65 0.02840   1
## 42       16    12.94 0.36300   1
```

Write out tables for paper

```r
write.csv( dplyr::filter(plasticuptable, fdr < 0.1, enrich>1) %>% arrange(desc(enrich)), file = "TableS1.csv")


write.csv( dplyr::filter(plasticdowntable, fdr < 0.1, enrich>1) %>% arrange(desc(enrich)), file = "TableS2.csv")
```
