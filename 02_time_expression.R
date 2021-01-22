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
ZMYND8_data <- read.csv("ZMYND8_metadata.csv", header = TRUE, row.names = 1)
counts_matrix <- read.csv(file.path(path.data, "expression_matrix.csv"), header= FALSE, row.names= 1)
rows_metadata <- read.csv(file.path(path.data, "rows_metadata.csv"))


#### Neocortex expression 

NCX <- ZMYND8_data %>%
  dplyr::filter(regions == "NCX") %>%
  ggplot(aes(stage, ZMYND8_expression))  + 
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
                                                   "s13" = ">=40 yrs"))

### Other regions 
OR <- ZMYND8_data %>%
  dplyr::filter(regions == "Other_regions") %>%
  ggplot(aes(stage, ZMYND8_expression))  + 
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
                                                   "s13" = ">=40 yrs"))




###### Makee all plots on the same page 

pdf("ZMYND8_stage_expression_brainSpan.pdf",width =14, height =12)
grid.arrange(NCX, OR, ncol=2, nrow=2)
dev.off()




#### Correlation with developmental age 


##### neocortex 
NCX_samples <- columns_metadata %>% 
  dplyr::filter(regions == "NCX") 

NCX_matrix <- counts_matrix %>% 
  dplyr::select(NCX_samples$name_for_matrix)

counts_ZMYND8_NCX <- NCX_matrix[rownames(NCX_matrix) %in% c("ENSG00000101040"),]

NCX_pval <- apply(t(counts_ZMYND8_NCX),2, function(x) cor.test(x, NCX_samples$ages, method="s"))
developmental_stages_NCX <- do.call("rbind", NCX_pval)


##### Other regions 
OtherRegions_samples <-  columns_metadata %>% 
  dplyr::filter(regions == "Other_regions")
OtherRegions_matrix <- counts_matrix %>%
  dplyr::select(OtherRegions_samples$df_name)

counts_ZMYND8_OR <- OtherRegions_matrix[rownames(OtherRegions_matrix) %in% c("ENSG00000101040"),]
OR_pval <- apply(t(counts_ZMYND8_OR),2, function(x) cor.test(x, OtherRegions_samples$ages, method="s"))
developmental_stages_OR <- do.call("rbind", OR_pval)
