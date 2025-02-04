library("ggplot2")
library("here")
library("RColorBrewer")
require("gplots")
require(knitr)
require(pheatmap)
require(readr)
require(cowplot)
require(data.table)
require(dplyr)
require(tidyverse)
require(ggrepel)


# Collection of senescence sets
sene_db = readRDS(here("Data","senescence_signatures.rds"))
length(sene_db)
#sapply(sene_db,length) %>% sort()
names(sene_db)


# Fig 3A
lapply(sene_db, length) %>% plyr::ldply(.id = "signature") %>% 
  ggplot(aes(reorder(signature, V1), V1)) +
  geom_col(fill="lightblue") + 
  geom_text(aes(label=round(V1, 2)) )+
  labs(x="", y="No. of genes") +
  coord_flip() +
  theme_minimal_grid() -> p
print(p)

# Fig 3B
shared_genes = expand.grid(names(sene_db), names(sene_db))%>% as.data.frame()
colnames(shared_genes) <- c("V1", "V2")

shared_genes$shared = 0; shared_genes$jaccard = 0
for(i in 1:nrow(shared_genes)){
  n_op = 0; jac = 0
  n1 = sene_db[[shared_genes[i,1]]]
  n2 = sene_db[[shared_genes[i,2]]]
  n_op = intersect(n1, n2)
  if(length(n_op) > 0) {
    jac = length(n_op)/(length(n1) + length(n2) - length(n_op))
  }
  shared_genes[i,"shared"] = length(n_op)
  shared_genes[i,"jaccard"] = jac
}

shared_genes %>% 
  ggplot(aes(reorder(V1, jaccard), reorder(V2, jaccard))) +
  geom_tile(aes(fill=jaccard)) + 
  scale_fill_gradient2(low = "blue", mid = "#FFFFFF", high = "#FF0000", midpoint = 0.1, space = "Lab", na.value = "grey50", guide = "colourbar")+
  geom_text(aes(label=round(jaccard,2)), size=3, color="black")+
  labs(x="", y="") +  coord_equal() +
  theme_minimal_grid() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) -> p

print(p)

# Fig 3C-J
auc_all <- read.csv(here("Data","Combine_4dataset_AUC.csv"),row.names = 1, header = T, check.names = F)

auc_all %>% rownames_to_column() %>% 
  reshape2::melt() -> dat_auc 

plot_list <- lapply(unique(dat_auc$variable), function(itype){
  p <- dat_auc %>% 
    dplyr::filter(variable == itype) %>% 
    ggplot(aes(reorder(rowname, value), value, fill = value)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(value, 2))) +
    coord_flip() +
    labs(x = "", y = "AUC", title = itype) +
    ggsci::scale_fill_gsea() +
    theme_minimal_grid()
  
  return(p)  
})


#Example: p[[1]]


# Fig 3K
# Comparisons of 14 signatures

apply(-auc_all, 2, rank) %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  reshape2::melt() %>% 
  ggplot(aes(reorder(rowname, -value, FUN=median ), value, color= variable)) +
  geom_boxplot(outliers = F, color="darkblue") +
  ggbeeswarm::geom_beeswarm(size=2) +
  ylim(0, 20) +
  labs(y="Rank of AUCs", x="") +
  coord_flip()+ 
  scale_color_brewer(palette = "Dark2") +
  ggprism::theme_prism() -> p
print(p)

