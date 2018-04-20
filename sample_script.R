### This script may be used to reproduce most of the results in
### Hofmeyr, Pavlidis and Eckley (2017) "Minimum Spectral Connectivity Projection Pursuit", https://arxiv.org/pdf/1509.01546.pdf
### minor adjustments to the implementation since the latest ArXiv update mean a few results differ slightly

### Install and load the SCPP package:
if(!("devtools"%in%installed.packages())) install.packages("devtools")
library(devtools)
if(!("SCPP"%in%installed.packages())) install_github("DavidHofmeyr/SCPP")
library(SCPP)

### Some of the datasets considered in the paper were too large to fit in the R package. These can be obtained from
# smartphone, isolet: https://archive.ics.uci.edu/ml/datasets.html
# Yale faces B: https://cervisia.org/machine_learning_data.php/
# phoneme: https://web.stanford.edu/~hastie/ElemStatLearn/

### For the remaining datasets:

## Optical recognition of handwritten digits dataset

data("optidigits")

sol <- SCPP_cluster(optidigits$x, 10)

cluster_performance(sol$cluster, optidigits$c)

## Pen based recognition of handwritten digits dataset

data("pendigits")

sol <- SCPP_cluster(pendigits$x, 10)

cluster_performance(sol$cluster, pendigits$c)

## Multiple feature digits dataset

data("mfdigits")

sol <- SCPP_cluster(mfdigits$x, 10)

cluster_performance(sol$cluster, mfdigits$c)

## Satellite dataset

data("satellite")

sol <- SCPP_cluster(satellite$x, 6)

cluster_performance(sol$cluster, satellite$c)

## Image segmentation dataset

data("imageseg")

sol <- SCPP_cluster(imageseg$x, 7)

cluster_performance(sol$cluster, imageseg$c)

## Breast cancer dataset

data("breastcancer")

sol <- SCPP_cluster(breastcancer$x, 2)

cluster_performance(sol$cluster, breastcancer$c)

## Synthetic control chart dataset

data("chart")

sol <- SCPP_cluster(chart$x, 6)

cluster_performance(sol$cluster, chart$c)

## Dermatology dataset

data("dermatology")

sol <- SCPP_cluster(dermatology$x, 6)

cluster_performance(sol$cluster, dermatology$c)

## Yeast dataset

data("yeast")

sol <- SCPP_cluster(yeast$x, 5)

cluster_performance(sol$cluster, yeast$c)
