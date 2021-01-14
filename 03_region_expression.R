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
  library(gridExtra)})




path.data <- path.expand("BrainSpan/Kang")

#### upload all data
counts_matrix <- read.csv(file.path(path.data, "expression_matrix.csv"), header= FALSE, row.names= 1)
columns_metadata <- read.csv("updated_metadata.csv", header = TRUE, row.names = 2)
rows_metadata <- read.csv(file.path(path.data, "rows_metadata.csv"))



### matrix 

colnames(counts_matrix) <- columns_metadata$name_for_matrix
rownames(counts_matrix) <- rows_metadata$ensembl_gene_id

####

### subset matrices based on period 

prenatal <- columns_metadata %>%
  dplyr::filter(period == "prenatal")

postnatal <- columns_metadata %>%
  dplyr::filter(period == "postnatal")



### subset matrices 
prenatal_matrix <- counts_matrix %>% 
  dplyr::select(prenatal$name_for_matrix)


postnatal_matrix <- counts_matrix %>% 
  dplyr::select(postnatal$name_for_matrix)

#### check how many regions are in each period 

pren_R<- as.data.frame(table(columns_metadata$structure_acronym[grep("pcw", columns_metadata$age)]))
postn_R <- as.data.frame(table(columns_metadata$structure_acronym[-grep("pcw", columns_metadata$age)]))
which(pren_R[,2]>=10)
which(postn_R[,2]>=10)

pren_R_structres <- pren_R %>% 
  dplyr::filter(Freq >=10)


postn_R_structures <- postn_R %>% 
  dplyr::filter(Freq >= 10)


prenatal_structures <- prenatal[prenatal$structure_acronym %in% pren_R_structres$Var1,]

table(prenatal_structures$structure_acronym)

prenatal_S_matrix <- prenatal_matrix %>% 
  dplyr::select(prenatal_structures$name_for_matrix)
dim(prenatal_S_matrix)


postnatal_structures <- postnatal[postnatal$structure_acronym %in% postn_R_structures$Var1,]

postnatal_S_matrix <- postnatal_matrix %>% 
  dplyr::select(postnatal_structures$name_for_matrix)

#### plot 

PreN_fig <- prenatal_S_matrix["ENSG00000101040",] %>%
  as.data.frame() %>% 
  rownames_to_column("geneID") %>%
  melt() %>% 
  cbind(prenatal_structures) %>%
  ggplot(aes(structure_acronym, value))  + 
  geom_boxplot(aes(fill=structure_acronym), alpha=0.65)  + 
  geom_jitter(size=2.5) +
  scale_color_viridis(discrete = TRUE, option = "magma")+ theme_bw() +
  scale_fill_viridis(discrete = TRUE, direction = 1)  +
  xlab("Stage") + ylab("RPKM") +
  theme(legend.position = "none", 
        axis.text=element_text(size=11, face="bold", color="black"), 
        axis.title = element_text(siz=12, face="bold", color="black"), 
        axis.text.x = element_text(angle=90)) + ggtitle("Prenatal Expression")


PN_fig <- postnatal_S_matrix["ENSG00000101040",] %>%
  as.data.frame() %>% 
  rownames_to_column("geneID") %>%
  melt() %>% 
  cbind(postnatal_structures) %>%
  ggplot(aes(structure_acronym, value))  + 
  geom_boxplot(aes(fill=structure_acronym), alpha=0.65)  + 
  geom_jitter(size=2.5) +
  scale_color_viridis(discrete = TRUE, option = "magma")+ theme_bw() +
  scale_fill_viridis(discrete = TRUE, direction = 1)  +
  xlab("Stage") + ylab("RPKM") +
  theme(legend.position = "none", 
        axis.text=element_text(size=11, face="bold", color="black"), 
        axis.title = element_text(siz=12, face="bold", color="black"), 
        axis.text.x = element_text(angle=90)) + ggtitle("Postnatal Expression")




pdf("ZMYND8_region_expression_brainSpan.pdf",width =14, height =12)
grid.arrange(PreN_fig, PN_fig, ncol=2, nrow=2)
dev.off()
