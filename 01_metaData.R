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
  library(RColorBrewer)
  library(stringr)})

rm(list=ls())

##### FILES #####
path.data <- path.expand("BrainSpan/Kang")


### Load metadata file 
columns_metadata <- read.csv(file.path(path.data, "columns_metadata.csv"), header = TRUE)
counts_matrix <- read.csv(file.path(path.data, "expression_matrix.csv"), header= FALSE, row.names= 1)
rows_metadata <- read.csv(file.path(path.data, "rows_metadata.csv"))

### Add stages based on the technical white paper 
columns_metadata$age <-  gsub(" ","-", columns_metadata$age )



levels(columns_metadata$age) <- c("8-pcw","12-pcw", "13-pcw", "16-pcw", 
                                  "17-pcw", "19-pcw", "21-pcw",
                                  "24-pcw", "25-pcw", "26-pcw", "35-pcw", 
                                  "37-pcw", "4-mos", "10-mos", "1-yrs", 
                                  "2-yrs","3-yrs", "4-yrs", "8-yrs", "11-yrs", 
                                  "13-yrs", "15-yrs", "18-yrs", "19-yrs", "21-yrs", 
                                  "23-yrs", "30-yrs", "36-yrs", "37-yrs", "40-yrs")


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


columns_metadata$numeric_age[grepl("pcw", columns_metadata$age, ignore.case = TRUE)]<-
   columns_metadata$age[grepl("pcw", columns_metadata$age)] %>%  str_remove("-pcw")%>% 
  as.numeric() %>% divide_by(-52)
                              


columns_metadata$numeric_age[grepl("mos", columns_metadata$age, ignore.case = TRUE)] <- 
  columns_metadata$age[grepl("-mos", columns_metadata$age)] %>%  str_remove("-mos") %>% 
  as.numeric() %>% divide_by(12)

columns_metadata$numeric_age[grepl("yrs", columns_metadata$age, ignore.case = TRUE)] <- 
  columns_metadata$age[grepl("-yrs", columns_metadata$age)] %>%  str_remove("-yrs") %>% 
  as.numeric 



##### Create name for column data for the matrix 


columns_metadata <- columns_metadata %>% 
  dplyr::mutate(name_for_matrix = paste(donor_id, age, structure_acronym, stage, sep = "_"))




### Add ZMYND8 expression 

colnames(counts_matrix) <- columns_metadata$name_for_matrix
rownames(counts_matrix) <- rows_metadata$ensembl_gene_id



ZMYND8_metaData <- counts_matrix["ENSG00000101040",] %>%
  as.data.frame() %>% 
  rownames_to_column("geneID") %>%
  melt() %>% 
  cbind(columns_metadata) %>% 
  dplyr::rename(ZMYND8_expression = value) %>%
  dplyr::select(-c(variable, column_num))



write.csv(ZMYND8_metaData, "ZMYND8_metadata.csv")
