#Analysis of untargeted LC-MS metabolome data collected across the succession of marine biofilm 

##Pre-processing of LC-MS/MS data
Raw HRMS data was converted to mzML using MSConvert (ProteoWizard) and preprocessed using MZmine 3 [74]. Molecular networking (MN) was all completed within the GNPS platform [75], which includes: Feature Based Molecular Networking [76] and Ion Identity Molecular Networking [77]. Ion Identity MN workflow can be found here: (To be filled out). Formula predictions, compound class annotations and NPclass annotations were run through SIRIUS/v5.5.7 and include: SIRIUS, ZODIAC, CSI:FingerID and CANOPUS. All data run in this study can be found in GNPS MassIVE (To be filled out). 

Formula predictions (SIRIUS) and the feature LC-MS table was imported and analysed further in R. 

##Load dependencies

```
library("dplyr")
library("phyloseq")
library(tidyverse)
library("devtools")
library("Biostrings")
library(ggplot2)
library(vegan)
library(ggvenn)
library(zCompositions)
library(compositions)
library(ggbeeswarm)
library(FSA)
```
