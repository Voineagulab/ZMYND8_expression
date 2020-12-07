#### Load required libraries 

library(tibble)
library(dplyr)
library(ggplot2)
library(viridis)
library(reshape2)
library(magrittr)
library(gridExtra)


# LOAD BRAINSPAN EXPRESSSION DATA
## expand the path where the BrainSpan data has been stored 
path.data <- path.expand("BrainSpan/Kang")

### expression matrix 
counts_matrix <- read.csv(file.path(path.data, "expression_matrix.csv"), header= FALSE, row.names= 1)

### columns metadata 
columns_metadata <- read.csv(file.path(path.data, "columns_metadata.csv"), header = TRUE)

### rows metadata
rows_metadata <- read.csv(file.path(path.data, "rows_metadata.csv"))


# ORGANISE COLUMN NAMES - stages and anatomy 

columns_metadata$age <-  gsub(" ","-", columns_metadata$age )

## STAGES 
### organise prenatal stages 
columns_metadata$stage[columns_metadata$age == "8-pcw" | columns_metadata$age == "9-pcw"] <- "s02a" 
columns_metadata$stage[columns_metadata$age == "12-pcw"] <- "s02b"
columns_metadata$stage[columns_metadata$age == "13-pcw"] <- "s03a"
columns_metadata$stage[columns_metadata$age == "16-pcw" | columns_metadata$age == "17-pcw"] <-"s03b"
columns_metadata$stage[columns_metadata$age == "19-pcw" | columns_metadata$age == "21-pcw" | columns_metadata$age == "24-pcw"] <- "s04"
columns_metadata$stage[columns_metadata$age == "25-pcw" | columns_metadata$age == "26-pcw" | columns_metadata$age == "35-pcw" | columns_metadata$age == "37-pcw"] <- "s05"


### organise postnatal stages 
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



### levels 
levels(columns_metadata$stage) <- c("s02a", "s02b", "s03a", "s03b", "s04", "s05", "s06", 
                                    "s07", "s08", "s09", "s10", "s11", "s12", "s13")




## REGIONS 

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


### Mutate nnames to create unique names based on structure acronym, stage and donor ID for 
### column names for the expression matrix 
columns_metadata <- columns_metadata %>% 
  dplyr::mutate(name_for_matrix = paste(donor_id, age, structure_acronym, stage, sep = "_"))


### Add gene names as rownames and 
### mutated names as colnames 
colnames(counts_matrix) <- columns_metadata$name_for_matrix
rownames(counts_matrix) <- rows_metadata$ensembl_gene_id


### PLOT ZMYND8 EXPRESSION 


#### neocortex
NCX <- counts_matrix["ENSG00000101040",] %>%
  as.data.frame() %>% 
  rownames_to_column("geneID") %>%
  melt() %>% 
  cbind(columns_metadata) %>%
  dplyr::filter(regions == "NCX") %>%
  ggplot(aes(stage, value))  + 
  geom_boxplot(aes(fill=stage), alpha=0.65)  + 
  geom_jitter(size=2.5) +
  scale_color_viridis(discrete = TRUE, option = "magma")+ theme_bw() +
  scale_fill_viridis(discrete = TRUE, direction = 1)  +
  xlab("Stage") + ylab("RPKM") +
  theme(legend.position = "none", 
        axis.text=element_text(size=11, face="bold", color="black"), 
        axis.title = element_text(siz=12, face="bold", color="black"), 
        axis.text.x = element_text(angle=90)) +
  ggtitle("Neocortex") + scale_x_discrete(labels=c("s02a" = "8 - 9pcw", 
                                                   "s02b" = "10 - 12pcw", 
                                                   "s03a" = "13 - 15pcw", 
                                                   "s03b" = "16 - 18pcw",
                                                   "s04" =  "19 - 24pcw", 
                                                   "s05" = "25 - 38pcw", 
                                                   "s06" = "Birth - 5 mos", 
                                                   "s07" = "6 - 18 mos", 
                                                   "s08" = "19 mos - 5yrs", 
                                                   "s09" = "6 - 11 yrs",
                                                   "s10" = "12 - 19 yrs",
                                                   "s11" = "20 - 29 yrs", 
                                                   "s12" = "30 - 39 yrs", 
                                                   "s13" = "40 - 49 yrs"))


### Other regions 
OR <- counts_matrix["ENSG00000101040",] %>%
  as.data.frame() %>% 
  rownames_to_column("geneID") %>%
  melt() %>% 
  cbind(columns_metadata) %>%
  dplyr::filter(regions == "Other_regions") %>%
  ggplot(aes(stage, value))  + 
  geom_boxplot(aes(fill=stage), alpha=0.65)  + 
  geom_jitter(size=2.5) +
  scale_color_viridis(discrete = TRUE, option = "magma")+ theme_bw() +
  scale_fill_viridis(discrete = TRUE, direction = 1)  +
  xlab("Stage") + ylab("RPKM") +
  theme(legend.position = "none", 
        axis.text=element_text(size=11, face="bold", color="black"), 
        axis.title = element_text(siz=12, face="bold", color="black"), 
        axis.text.x = element_text(angle=90)) +
  ggtitle("Other regions") + scale_x_discrete(labels=c("s02a" = "8 - 9pcw", 
                                                   "s02b" = "10 - 12pcw", 
                                                   "s03a" = "13 - 15pcw", 
                                                   "s03b" = "16 - 18pcw",
                                                   "s04" =  "19 - 24pcw", 
                                                   "s05" = "25 - 38pcw", 
                                                   "s06" = "Birth - 5 mos", 
                                                   "s07" = "6 - 18 mos", 
                                                   "s08" = "19 mos - 5yrs", 
                                                   "s09" = "6 - 11 yrs",
                                                   "s10" = "12 - 19 yrs",
                                                   "s11" = "20 - 29 yrs", 
                                                   "s12" = "30 - 39 yrs", 
                                                   "s13" = "40 - 49 yrs"))

### SAVE AS PDF 
pdf("ZMYND8_stage_expression_brainSpan.pdf",width =14, height =12)
grid.arrange(NCX, OR, ncol=2, nrow=2)
dev.off()
