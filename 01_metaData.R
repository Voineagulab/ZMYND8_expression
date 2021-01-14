##### Adding additional information in the metadata file before perfoming any analysis 


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(tools)
  library(magrittr)
  library(tibble)
  library(plyr)
  library(pheatmap)
  library(readxl)
  library(gdata)
  library(biomaRt)
  library(data.table)
  library(pander)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)})

rm(list=ls())

##### FILES #####
path.data <- path.expand("BrainSpan/Kang")


### Load metadata file 
columns_metadata <- read.csv(file.path(path.data, "columns_metadata.csv"), header = TRUE)


### Add stages based on the technical white paper 
columns_metadata$age <-  gsub(" ","-", columns_metadata$age )

columns_metadata$stage[columns_metadata$age == "8-pcw" | columns_metadata$age == "9-pcw"] <- "s02a" 
columns_metadata$stage[columns_metadata$age == "12-pcw"] <- "s02b"
columns_metadata$stage[columns_metadata$age == "13-pcw"] <- "s03a"
columns_metadata$stage[columns_metadata$age == "16-pcw" | columns_metadata$age == "17-pcw"] <-"s03b"
columns_metadata$stage[columns_metadata$age == "19-pcw" | columns_metadata$age == "21-pcw" | columns_metadata$age == "24-pcw"] <- "s04"
columns_metadata$stage[columns_metadata$age == "25-pcw" | columns_metadata$age == "26-pcw" | columns_metadata$age == "35-pcw" | columns_metadata$age == "37-pcw"] <- "s05"
columns_metadata$stage[columns_metadata$age == "4-mos"] <- "s06"
columns_metadata$stage[columns_metadata$age == "10-mos" | columns_metadata$age == "1-yrs"] <- "s07"
columns_metadata$stage[columns_metadata$age == "2-yrs" | columns_metadata$age == "3-yrs" | columns_metadata$age == "4-yrs"] <- "s08"
columns_metadata$stage[columns_metadata$age == "8-yrs"| columns_metadata$age == "11-yrs"] <- "s09"
columns_metadata$stage[columns_metadata$age == "13-yrs" | columns_metadata$age == "15-yrs" | columns_metadata$age == "18-yrs" | columns_metadata$age == "19-yrs"] <-  "s10"
columns_metadata$stage[columns_metadata$age == "21-yrs" |
                         columns_metadata$age == "23-yrs"] <- "s11"
columns_metadata$stage[columns_metadata$age == "30-yrs" |
                         columns_metadata$age == "36-yrs" | columns_metadata$age == "37-yrs"] <- "s12"
columns_metadata$stage[columns_metadata$age == "40-yrs"] <- "s13"

levels(columns_metadata$stage) <- c("s02a", "s02b", "s03a", "s03b", "s04", "s05", "s06", 
                                    "s07", "s08", "s09", "s10", "s11", "s12", "s13")


levels(columns_metadata$age) <- c("8-pcw","12-pcw", "13-pcw", "16-pcw", 
                                  "17-pcw", "19-pcw", "21-pcw",
                                  "24-pcw", "25-pcw", "26-pcw", "35-pcw", 
                                  "37-pcw", "4-mos", "10-mos", "1-yrs", 
                                  "2-yrs","3-yrs", "4-yrs", "8-yrs", "11-yrs", 
                                  "13-yrs", "15-yrs", "18-yrs", "19-yrs", "21-yrs", 
                                  "23-yrs", "30-yrs", "36-yrs", "37-yrs", "40-yrs")


### Period: Either prenatal or postnatal 

columns_metadata$period[columns_metadata$stage == "s02a"|
                          columns_metadata$stage == "s02b"|
                          columns_metadata$stage == "s03a"|
                          columns_metadata$stage == "s03b"|
                          columns_metadata$stage == "s04"|
                          columns_metadata$stage == "s05"] <- "prenatal"

columns_metadata$period[columns_metadata$stage == "s06"|
                          columns_metadata$stage == "s07"|
                          columns_metadata$stage == "s08"|
                          columns_metadata$stage == "s09"|
                          columns_metadata$stage == "s10"|
                          columns_metadata$stage == "s11" |
                          columns_metadata$stage == "s12" |
                          columns_metadata$stage == "s13"] <- "postnatal"

#### Regions: Either neocortex (NCX) or all other regions (other regions)


#########  regions 

columns_metadata$regions[columns_metadata$structure_acronym == "OFC" | 
                           columns_metadata$structure_acronym == "DFC" |
                           columns_metadata$structure_acronym == "VFC" |
                           columns_metadata$structure_acronym ==  "MFC" | 
                           columns_metadata$structure_acronym == "M1C" |  ## frontal cortex
                           columns_metadata$structure_acronym == "M1C-S1C" | 
                           columns_metadata$structure_acronym == "IPC" | 
                           columns_metadata$structure_acronym == "S1C" |  ### parietal cortex 
                           columns_metadata$structure_acronym == "A1C" | 
                           columns_metadata$structure_acronym == "STC" |  
                           columns_metadata$structure_acronym == "ITC" | ##temporal cortex
                           columns_metadata$structure_acronym == "V1C" |
                           columns_metadata$structure_acronym == "Ocx" |
                           columns_metadata$structure_acronym == "PCx" | 
                           columns_metadata$structure_acronym == "TCx" ] <- "NCX"


columns_metadata$regions[columns_metadata$structure_acronym == "MGE" | 
                           columns_metadata$structure_acronym == "LGE" | 
                           columns_metadata$structure_acronym == "CGE" | 
                           columns_metadata$structure_acronym == "STR" |
                           columns_metadata$structure_acronym ==  "HIP" | 
                           columns_metadata$structure_acronym == "AMY" |
                           columns_metadata$structure_acronym == "DTH" | 
                           columns_metadata$structure_acronym == "MD" |
                           columns_metadata$structure_acronym == "URL" | 
                           columns_metadata$structure_acronym == "CBC" | 
                           columns_metadata$structure_acronym == "CB"] <- "Other_regions"




#### NUMERIC AGES 
###### If age contains -pcw , then remove -pcw and divide age by 52 and add negative 
###### If age contains mos, then remove -mos and divide age by 12
##### else remove -yrs



columns_metadata$ages[columns_metadata$age == "8-pcw"] <- -0.15

columns_metadata$ages[columns_metadata$age == "9-pcw"] <- -0.17

columns_metadata$ages[columns_metadata$age == "12-pcw"] <- -0.23

columns_metadata$ages[columns_metadata$age == "13-pcw"] <- -0.25

columns_metadata$ages[columns_metadata$age == "16-pcw"] <- -0.31

columns_metadata$ages[columns_metadata$age == "17-pcw"] <- -0.32

columns_metadata$ages[columns_metadata$age == "19-pcw"] <- -0.36

columns_metadata$ages[columns_metadata$age == "21-pcw"] <- -0.40

columns_metadata$ages[columns_metadata$age == "24-pcw"] <- -0.46

columns_metadata$ages[columns_metadata$age == "25-pcw"] <- -0.48

columns_metadata$ages[columns_metadata$age == "26-pcw"] <- -0.498

columns_metadata$ages[columns_metadata$age == "35-pcw"] <- -0.67

columns_metadata$ages[columns_metadata$age == "37-pcw"] <- -0.71

columns_metadata$ages[columns_metadata$age == "4-mos"] <- 0.33
columns_metadata$ages[columns_metadata$age == "10-mos"] <- 0.83
columns_metadata$ages[columns_metadata$age == "1-yrs"] <- 1
columns_metadata$ages[columns_metadata$age == "2-yrs"] <- 2
columns_metadata$ages[columns_metadata$age == "3-yrs"] <- 3
columns_metadata$ages[columns_metadata$age == "4-yrs"] <- 4
columns_metadata$ages[columns_metadata$age == "8-yrs"] <- 8
columns_metadata$ages[columns_metadata$age == "11-yrs"] <- 11
columns_metadata$ages[columns_metadata$age == "13-yrs"] <- 13
columns_metadata$ages[columns_metadata$age == "15-yrs"] <- 15
columns_metadata$ages[columns_metadata$age == "18-yrs"] <- 18
columns_metadata$ages[columns_metadata$age == "19-yrs"] <-  19
columns_metadata$ages[columns_metadata$age == "21-yrs"] <- 21
columns_metadata$ages[columns_metadata$age == "22-yrs"] <- 22
columns_metadata$ages[columns_metadata$age == "23-yrs"] <- 23
columns_metadata$ages[columns_metadata$age == "30-yrs"] <- 30
columns_metadata$ages[columns_metadata$age == "36-yrs"] <- 36
columns_metadata$ages[columns_metadata$age == "37-yrs"] <- 37
columns_metadata$ages[columns_metadata$age == "40-yrs"] <- 40




##### Create name for column data for the matrix 


columns_metadata <- columns_metadata %>% 
  dplyr::mutate(name_for_matrix = paste(donor_id, age, structure_acronym, stage, sep = "_"))



write.csv(columns_metadata, "updated_metadata.csv")
