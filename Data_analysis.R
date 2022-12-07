library(ape)
library(ggtree)
library(phylolm)
library(ggplot2)
library(ggpubr)

## LOADING TREE ----------------------------------------------------------------

setwd("~/euklen/")

tree <- read.tree("busco_tree_pretrim.newick")
tree$tip.label <- gsub(".cds.filtered","",tree$tip.label)
tree$tip.label <- tree$tip.label %>% substr(1,15)

## LOADING DATA ----------------------------------------------------------------

setwd("~/euklen/")
growth_data <- read.csv("euk_Weissman.csv") %>% subset(Gene.Length>1e6)
growth_data_deduplicated <- subset(growth_data,!duplicated(Accession))
growth_data_sub <- subset(growth_data_deduplicated,Gene.Length>1e6)
growth_data_sub$Accession <- substr(growth_data_sub$Accession,1,15)
tree_sub <- drop.tip(tree, which(!(tree$tip.label %in% growth_data_sub$Accession)))
rownames(growth_data_sub) <- growth_data_sub$Accession

## LOADING TREE MMETSP ---------------------------------------------------------

setwd("~/euklen/")
tree.mmestp <- read.tree("SpeciesTree_rooted.txt")
tree.mmestp$tip.label <- gsub(".cds.filtered","",tree.mmestp$tip.label)

## LOADING DATA MMESTP ---------------------------------------------------------

setwd("~/euklen/")
growth_data.mmestp <- read.csv("growth.csv")
growth_data_deduplicated.mmestp <- subset(growth_data.mmestp,!duplicated(Accession))
growth_data_sub.mmestp <- subset(growth_data_deduplicated.mmestp,Number.Genes.Filtered>2000)
tree_sub.mmestp <- drop.tip(tree.mmestp, which(!(tree.mmestp$tip.label %in% growth_data_sub.mmestp$Accession)))
rownames(growth_data_sub.mmestp) <- growth_data_sub.mmestp$Accession

## Genome Length Data ----------------------------------------------------------

setwd("~/euklen/")
nvl_data <- read.csv("euk_genomevsnogene.csv")

## Plot ------------------------------------------------------------------------

NvD <- ggplot(growth_data.mmestp, 
              aes(x = Number.Genes.Filtered, y = Doubling.Time)) +
  geom_point() + 
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method="lm",color="gray") +
  theme_pubclean() +
  xlab("Number of Expressed Genes") +
  ylab("Doubling Time (Hours)") +
  ggtitle("MMETSP")

LvD <- ggplot(growth_data, 
              aes(x = Gene.Length, y = Doubling.Time)) +
  geom_point() + 
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method="lm",color="gray") +
  theme_pubclean() +
  xlab("Genome Length") +
  ylab("Doubling Time (Hours)") +
  ggtitle("GenBank")

NvL <- ggplot(nvl_data, 
       aes(x = Genome.Length*1e6, y = No.gene)) +
  geom_point() + 
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method="lm",color="gray") +
  theme_pubclean() +
  xlab("Genome Length") +
  ylab("Number of Genes")

setwd("~/euklen/")
png("Fig1.png",width=700,height=500)
ggarrange(LvD,
          ggarrange(NvL,NvD,nrow=2,labels=c("(b)","(c)")),
          labels=c("(a)",""),
          hjust=-1,
          ncol=2)
dev.off()

setwd("~/euklen/")
pdf("Fig1.pdf",width=8,height=5)
ggarrange(LvD,
          ggarrange(NvL,NvD,nrow=2,labels=c("(b)","(c)")),
          labels=c("(a)",""),
          hjust=-1,
          ncol=2)
dev.off()


## Statistics ------------------------------------------------------------------

### GenBank

nophylo_model <- lm(log10(Doubling.Time) ~ log10(Gene.Length) + OGT,
                    data=growth_data_sub)
summary(nophylo_model)

phylo_model<- phylolm(log10(Doubling.Time) ~ log10(Gene.Length) + OGT,
                              data=growth_data_sub,
                              phy =  tree_sub,
                              model="BM")
summary(phylo_model)

### MMETSP

nophylo_model.mmestp <- lm(log10(Doubling.Time) ~ log10(Number.Genes.Filtered) + OGT,
                    data=growth_data_sub.mmestp)
summary(nophylo_model.mmestp)

phylo_model.mmestp <- phylolm(log10(Doubling.Time) ~ log10(Number.Genes.Filtered) + OGT,
                              data=growth_data_sub.mmestp,
                              phy=tree_sub.mmestp,
                              model="BM")
summary(phylo_model.mmestp)

### NvL

nvl_model <- lm(log10(Genome.Length) ~ log10(No.gene), data =nvl_data)
summary(nvl_model)

