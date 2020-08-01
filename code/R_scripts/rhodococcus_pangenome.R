######################################
### Rhodococcus pangenome analysis ###
######################################

# Ryan Blaustein
# ryan.blaustein@northwestern.edu
# ryan.blaustein1@gmail.com
# last edited: 07-31-2020

######################################

library(ggplot2)
library(reshape2)
library(heatmap.plus)
library(vegan)
library(ape)
library(phangorn)
library(seqinr)
library(gridExtra)

#####
### load data
#####


### pangenome files

## pangenome matrix
pangenome_mat = read.table("pangenome/roary_out_all/gene_presence_absence.Rtab", h=T, row.names=1)

## pangenome csv (for function binning)
pangenome_csv = read.csv("pangenome/roary_out_all/gene_presence_absence.csv", h=T, row.names=1)
# pull annotations from csv file
pangenome_anno = pangenome_csv[,1:2]

## total gene clusters
totalgenes = read.table("pangenome/roary_out_all/number_of_genes_in_pan_genome.Rtab")


### TG9 clade core gene alignment from pangenome subset

## core gene alignment
core_tree <- read.tree("pangenome/roary_out_3B/ad_hoc/core_gene_aln_tree.newick")


### TG9 files

## subset pangenome matrix around TG9 clade
pangenome_mat_3B = read.table("pangenome/roary_out_3B/gene_presence_absence.Rtab", h=T, row.names=1)

## fasta file -- amino acid seqs
TG9_aa <- read.fasta(file = "TG9_aa_seqs/TG9_prok.faa") # from prokka output for GCA_003288095 (TG9)

## eggnog cog category output
TG9_cog =  read.table("TG9_aa_seqs/TG9_cog_cat.txt", h=T, row.names=1, sep="\t") # from running eggnog


#####
### pangenome gene clusters summary
#####

### quick stats (mean, sd, se)
mean(apply(pangenome_mat, 2, sum))
sd(apply(pangenome_mat, 2, sum))
sd(apply(pangenome_mat, 2, sum))/sqrt(dim(pangenome_mat)[2])
# TG9
sum(pangenome_mat[, grep("GCA_003288095", colnames(pangenome_mat))])

# hypothetical protein
length(grep("hypothetical", pangenome_anno$Annotation)) # total
100*length(grep("hypothetical", pangenome_anno$Annotation))/dim(pangenome_anno)[1] # percent genes in pangenome
mean(apply(pangenome_mat[grep("hypothetical", pangenome_anno$Annotation),], 2, sum)) # average hp per genome


### total genes (FIGURE S1)
colnames(totalgenes) = c(1:dim(totalgenes)[2])
totalgenes_m = melt(totalgenes)
# plot
ggplot(data = totalgenes_m, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA, col = "darkgray", width = 0.5) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.7, fatten = 1.5) +
  xlab("No. of genomes") +
  ylab("Total genes") +
  scale_x_discrete(breaks=c(0, 25, 50, 75, 100, 125, 150)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=15, color = "black"), 
        axis.title = element_text(size=16, color = "black"))


### histogram disribution (FIGURE 2A)
genes_by_genome <- apply(pangenome_mat, 1, sum)
gene_freq = data.frame(genome = c(1:dim(pangenome_mat)[2]),
                       count = c(1:dim(pangenome_mat)[2]),
                       group = c(rep("Cloud", floor(dim(pangenome_mat)[2]*0.15)),
                                 rep("Shell", floor(dim(pangenome_mat)[2]*0.8)),
                                 rep("Core", ceiling(dim(pangenome_mat)[2]*0.05))))
gene_freq$group = factor(gene_freq$group, levels = c("Cloud", "Shell", "Core"))
for (i in 1:dim(pangenome_mat)[2]) {
  gene_freq[i,2] = length(which(genes_by_genome == i))
}
# plot
ggplot(gene_freq, aes(x = genome, y = count, fill = group
)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("No. of genomes") +
  ylab("Total genes") +
  xlim(0,150) +
  theme_classic() +
  scale_fill_manual(values=c("red", "orange","blue")) +
  theme(axis.text = element_text(size=15, color = "black"), 
        axis.title = element_text(size=16, color = "black")) +
  theme(legend.position="none")


### pie distribution (FIGURE 2A)
pie.df <- data.frame(group = c("Cloud", "Shell", "Core"),
                     count = tapply(gene_freq$count, gene_freq$group, sum),
                     percent = 100*tapply(gene_freq$count, gene_freq$group, sum)/sum(gene_freq$count))
pie.df$percent = round(pie.df$percent, 1)
pie.df$group = factor(pie.df$group, levels = c("Cloud", "Shell", "Core"))
# plot
ggplot(pie.df, aes(x="", y=percent, fill=group)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("red", "orange","blue")) +
  #geom_text(aes(y = c(68.2, 28.5, 3.3), 
  #              label = count), size=3) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.border = element_blank(), 
        panel.grid=element_blank(), axis.ticks = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=22, color = "black"))
# total count/percent per fraction
pie.df


### pan fractions for average genome and TG9 (FIGURE 2B)
# Assign genes to core/cloud/shell
gene_type = data.frame(genes_by_genome,
                       category = rep("Core", length(genes_by_genome)))
gene_type$category = as.character(gene_type$category)
gene_type$category[which(gene_type$genes_by_genome < 137)] = 
  rep("Shell", length(gene_type$category[which(gene_type$genes_by_genome < 137)]))
gene_type$category[which(gene_type$genes_by_genome < 22)] = 
  rep("Cloud", length(gene_type$category[which(gene_type$genes_by_genome < 22)]))
# Determine fraction count per genome                              
genes_per_genome = matrix(nrow = 144, ncol = 3)
for (i in 1:144) {
  genes_per_genome[i,] = rbind(tapply(pangenome_mat[,i], gene_type$category, sum))
}
colnames(genes_per_genome) = c("Cloud", "Core", "Shell")
rownames(genes_per_genome) = colnames(pangenome_mat)
head(genes_per_genome)
# stats for avegares and TG9
GPG_stat = data.frame(strain = c("Rhodococcus", "Rhodococcus", "Rhodococcus", "TG9", "TG9", "TG9"),
                      count = c(apply(genes_per_genome, 2, mean),
                                genes_per_genome[grep("GCA_003288095", rownames(genes_per_genome)),]),
                      category = c("Cloud", "Core", "Shell", "Cloud", "Core", "Shell"))
GPG_stat$category = factor(GPG_stat$category,
                           levels = c("Cloud", "Shell", "Core"))
# barplot
ggplot(GPG_stat, aes(x = strain, y = count, fill = category)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() + 
  scale_fill_manual(values = c("red", "orange","blue")) +
  ylab("Genes per genome") +
  theme(axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(size=14, color = "black"), 
        axis.title = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size=14, color = "black"), 
        legend.title = element_blank())
tapply(GPG_stat$count, GPG_stat$strain, sum)


#####
### binned gene products summary
#####

## bin annotations for presence/absence
fxn_df = as.data.frame(matrix(nrow = length(levels(pangenome_anno$Annotation)), ncol = ncol(pangenome_mat)))
for (i in 1:ncol(pangenome_mat)) {
  fxn_df[,i] = tapply(pangenome_mat[,i], pangenome_anno$Annotation,sum)
}
rownames(fxn_df) = names(tapply(pangenome_mat[,1], pangenome_anno$Annotation,sum))
colnames(fxn_df) = colnames(pangenome_mat)

## convert fxn_df to binary presence/absence data
fxn_df_bin = ifelse(fxn_df>0, 1, 0)
fxn_df_bin = as.data.frame(fxn_df_bin)

## quick stats (mean, sd, se)
dim(fxn_df_bin)
mean(apply(fxn_df_bin, 2, sum))
sd(apply(fxn_df_bin, 2, sum))
sd(apply(fxn_df_bin, 2, sum))/sqrt(dim(fxn_df_bin)[2])
# TG9
sum(fxn_df_bin[, grep("GCA_003288095", colnames(fxn_df_bin))])


### histogram disribution (FIGURE 2C)
fxns_by_genome <- apply(fxn_df_bin, 1, sum)
fxn_freq = data.frame(genome = c(1:dim(fxn_df_bin)[2]),
                      count = c(1:dim(fxn_df_bin)[2]),
                      group = c(rep("Cloud", floor(dim(fxn_df_bin)[2]*0.15)),
                                rep("Shell", floor(dim(fxn_df_bin)[2]*0.8)),
                                rep("Core", ceiling(dim(fxn_df_bin)[2]*0.05))))
fxn_freq$group = factor(fxn_freq$group, levels = c("Cloud", "Shell", "Core"))
for (i in 1:dim(fxn_df_bin)[2]) {
  fxn_freq[i,2] = length(which(fxns_by_genome == i))
}
# plot
ggplot(fxn_freq, aes(x = genome, y = count, fill = group
)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("No. of genomes") +
  ylab("Total gene products") +
  xlim(0,150) +
  theme_classic() +
  scale_fill_manual(values=c("red", "orange","blue")) +
  theme(axis.text = element_text(size=15, color = "black"), 
        axis.title = element_text(size=16, color = "black")) +
  theme(legend.position="none")


### pie distribution (FIGURE 2C)
pie.df_fxn <- data.frame(group = c("Cloud", "Shell", "Core"),
                         count = tapply(fxn_freq$count, fxn_freq$group, sum),
                         percent = 100*tapply(fxn_freq$count, fxn_freq$group, sum)/sum(fxn_freq$count))
pie.df_fxn$percent = round(pie.df_fxn$percent, 1)
pie.df_fxn$group = factor(pie.df_fxn$group, levels = c("Cloud", "Shell", "Core"))
# plot
ggplot(pie.df_fxn, aes(x="", y=percent, fill=group)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("red", "orange","blue")) +
  #geom_text(aes(y = c(68.2, 28.5, 3.3), 
  #              label = count), size=3) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.border = element_blank(), 
        panel.grid=element_blank(), axis.ticks = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=22, color = "black"))
# total count/percent per fraction
pie.df_fxn


### pan fractions for average genome and TG9 (FIGURE 2D)
# Assign genes to core/cloud/shell
fxn_type = data.frame(fxns_by_genome,
                      category = rep("Core", length(fxns_by_genome)))
fxn_type$category = as.character(fxn_type$category)
fxn_type$category[which(fxn_type$fxns_by_genome < 137)] = 
  rep("Shell", length(fxn_type$category[which(fxn_type$fxns_by_genome < 137)]))
fxn_type$category[which(fxn_type$fxns_by_genome < 22)] = 
  rep("Cloud", length(fxn_type$category[which(fxn_type$fxns_by_genome < 22)]))
# Determine fraction count per genome
tapply(fxn_df_bin[,1], fxn_type$category, sum)
fxns_per_genome = matrix(nrow = 144, ncol = 3)
for (i in 1:144) {
  fxns_per_genome[i,] = rbind(tapply(fxn_df_bin[,i], fxn_type$category, sum))
}
colnames(fxns_per_genome) = c("Cloud", "Core", "Shell")
rownames(fxns_per_genome) = colnames(fxn_df_bin)
head(fxns_per_genome)
# stats for avegares and TG9
FPG_stat = data.frame(strain = c("Rhodococcus", "Rhodococcus", "Rhodococcus", "TG9", "TG9", "TG9"),
                      count = c(apply(fxns_per_genome, 2, mean),
                                fxns_per_genome[grep("GCA_003288095", rownames(fxns_per_genome)),]),
                      category = c("Cloud", "Core", "Shell", "Cloud", "Core", "Shell"))
FPG_stat$category = factor(FPG_stat$category,
                           levels = c("Cloud", "Shell", "Core"))
# barplot
ggplot(FPG_stat, aes(x = strain, y = count, fill = category)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() + 
  scale_fill_manual(values = c("red", "orange","blue")) +
  ylab("Gene products per genome") +
  theme(axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(size=14, color = "black"), 
        axis.title = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size=14, color = "black"), 
        legend.title = element_blank())


### heatmap (FIGURE 3A)
heat_fxn = as.matrix(fxn_df_bin)
colnames(heat_fxn) = substr(colnames(heat_fxn), start = 1, stop = nchar(colnames(heat_fxn)) - 5)
dim(heat_fxn)
rowcol = rep("white", 144)
rowcol[grep("GCA_003288095", colnames(heat_fxn))] = c("red")
# plot
heatmap.plus(heat_fxn,
             distfun = function(x) vegdist(x, method = 'bray', binary = TRUE), 
             labRow = F, 
             cexCol = 0.45, 
             #labCol = F,
             col = c("yellow", "black"), 
             #ColSideColors = as.matrix(data.frame(rowcol, rowcol)),
             #margins = c(1,10),
             scale = "none")


#####
### TG9 clade subset
#####

### Core gene alignment (FIGURE 3C)

## set labels
core_tree$tip.label = substr(core_tree$tip.label, start = 1, stop = 13)
core_tree$tip.label = paste0("   ", core_tree$tip.label)

## set colors
col.tree = rep("gray", 15)
# gordoniae
col.tree[6] = c("blue")
# pyridinivorans
col.tree[c(1, 4, 13:15)] = rep("yellow", 5)
# rhodochrous
col.tree[c(7, 8, 10)] = rep("orange", 3)
# biphenylivorans
col.tree[11] = c("red")

## plot phylogram
plot(midpoint(core_tree),
     cex = 0.9,
     font = 1,
     no.margin = F,
     direction = "downwards",
     edge.width = 1.2)
tiplabels(bg = col.tree, pch = 22, cex = 2.2)
add.scale.bar(ask = TRUE, cex = 1.2)

## plot legend (print and add to tree figure)
plot(1:50, col = "white")
legend(1, 49, legend=c("Rhodococcus sp.", "R. biphenylivorans", "R. gordoniae",
                       "R. pyridinivorans", "R. rhodochrous"),
       c("gray", "red", "blue", "yellow", "orange"), 
       #pch=21, 
       text.font = 3,
       cex=1.5)

### histogram disribution of clade 3B (FIGURE 3B)
genes_by_genome <- apply(pangenome_mat_3B, 1, sum)
gene_freq = data.frame(genome = c(1:dim(pangenome_mat_3B)[2]),
                       count = c(1:dim(pangenome_mat_3B)[2]))
for (i in 1:dim(pangenome_mat_3B)[2]) {
  gene_freq[i,2] = length(which(genes_by_genome == i))
}
gene_freq
# plot
ggplot(gene_freq, aes(x = genome, y = count)) +
  geom_bar(stat = "identity", width = 0.95, fill = "maroon") +
  xlab("No. of genomes") +
  ylab("Total genes") +
  theme_classic() +
  theme(axis.text = element_text(size=15, color = "black"), 
        axis.title = element_text(size=16, color = "black")) +
  theme(legend.position="none")


#####
### TG9 gene product screen
#####

## set df for fxn list totals, presence in TG9, and total in cluster 3B
fxns_TG9 = data.frame(fxn_type, 
                      TG9 = fxn_df_bin[, grep("GCA_003288095",colnames(fxn_df_bin))],
                      Cluster_3B = apply(fxn_df_bin[,grep("GCA_900455725|GCA_004153645|GCA_001465325|GCA_009831075|GCA_003004765|GCA_002215235|GCA_900177695|GCA_001886355|GCA_001665495|GCA_003288095|GCA_000763325|GCA_000236965|GCA_003633655|GCA_001894905|GCA_900105195",
                                                          colnames(fxn_df_bin))], 1, sum))
fxns_TG9 = fxns_TG9[order(fxns_TG9$fxns_by_genome),]

## add totals by group
fxns_TG9$genomes_other = fxns_TG9$fxns_by_genome - fxns_TG9$Cluster_3B
fxns_TG9$Cluster_3B_other = fxns_TG9$Cluster_3B - fxns_TG9$TG9

##add percentages by group
#fxns_TG9$genome_percent = 100*fxns_TG9$fxns_by_genome/144
#fxns_TG9$Cluster_3B_percent = 100*fxns_TG9$Cluster_3B/15
fxns_TG9$genomes_other_percent = 100*fxns_TG9$genomes_other/129
fxns_TG9$Cluster_3B_other_percent = 100*fxns_TG9$Cluster_3B_other/14
fxns_TG9$genomes_total_percent = 100*fxns_TG9$fxns_by_genome/144

## write table
#write.table(fxns_TG9, 'pangenome/roary_out_all/ad_hoc/fxns_TG9.txt',sep="\t", col.names=NA, quote = FALSE)

## subset only fxns in TG9
fxns_TG9_sub = fxns_TG9[which(fxns_TG9$TG9 > 0),]
fxns_TG9_sub = fxns_TG9_sub[order(fxns_TG9_sub$fxns_by_genome),]

## screen for potential antibiotic resistance
fxns_TG9_sub[grep("esistance|lactam|acrolide|etracyclin|mycin|olymyxin|ifampin|uinolone|enicillin|inoglycoside|reptogramin|ricin", rownames(fxns_TG9_sub)),]
fxns_TG9[grep("esistance|lactam|acrolide|etracyclin|mycin|olymyxin|ifampin|uinolone|enicillin|inoglycoside|reptogramin|ricin", rownames(fxns_TG9)),]

# screen totals
dim(fxns_TG9_sub[grep("esistance|lactam|acrolide|etracyclin|mycin|olymyxin|ifampin|uinolone|enicillin|inoglycoside|reptogramin|ricin", rownames(fxns_TG9_sub)),])
dim(fxns_TG9[grep("esistance|lactam|acrolide|etracyclin|mycin|olymyxin|ifampin|uinolone|enicillin|inoglycoside|reptogramin|ricin", rownames(fxns_TG9)),])

#write.table(fxns_TG9[grep("esistance|lactam|acrolide|etracyclin|mycin|olymyxin|ifampin|uinolone|enicillin|inoglycoside|reptogramin|ricin", rownames(fxns_TG9)),], 
#            'pangenome/roary_out_all/ad_hoc/ARG_list_putative.txt',sep="\t", col.names=NA, quote = FALSE)
# NOTE: needs manual edit based on screening genes on UniProtKB (http://www.uniprot.org) to confirm antibiotic resistance as the biological process (i.e., ensure that it was not antibiotic biosynthesis)

## screen for other functions of interest
fxns_TG9[grep("rpf|Rpf", rownames(fxns_TG9)),]
fxns_TG9[grep("Biphenyl|biphenyl", rownames(fxns_TG9)),]
fxns_TG9[grep("Toxin|toxin", rownames(fxns_TG9)),]
fxns_TG9_sub[grep("rpf|Rpf", rownames(fxns_TG9_sub)),]
fxns_TG9_sub[grep("Biphenyl|biphenyl", rownames(fxns_TG9_sub)),]
fxns_TG9_sub[grep("Toxin|toxin", rownames(fxns_TG9_sub)),]


#####
### Resistome
#####

## total ARG gene clusters
resistome = read.table("pangenome/roary_out_all/ad_hoc/ARG_list_corrected.txt", sep='\t',h=T,row.names=1,check=F,comment='',quote="\"")

## barplot (FIGURE 5)
resistome$gene = factor(resistome$gene, levels = resistome$gene)
resistome$TG9_category = factor(resistome$TG9_category, 
                                levels = c("Core", "Accessory (TG9)", "Accessory (other strains)"))
resistome$putative_resistance_class = factor(resistome$putative_resistance_class,
                                             levels = c("aminoglycosides", "B-lactams", "macrolides", "polymyxins", 
                                                        "quinolones", "rifamycins", "streptogramins", "tetracyclines",  
                                                        "multidrug", "antimicrobial_other", "antitumor*", "heavy metal*"))

ggplot(resistome, aes(x=gene, y=percent_genomes, fill=putative_resistance_class)) +
  geom_bar(position=position_dodge(width=0.9), stat = "identity") +
  facet_grid(.~TG9_category, scale="free_x", space="free") +
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = c("#000033", "blue", "skyblue", "green", "yellow", "red", "purple", 
                               "magenta", "pink", "gray", "orange", "brown")) +
  ylab(expression(paste("Percent of ", italic("Rhodococcus"), " genomes"))) +
  xlab("Antimicrobial Resistance Genes") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank())


######
### Subsetting TG9 amino-acid sequences by pangenome fraction
#####

## amino acid sequences corresponding to pangenome components (core, shell, cloud)
head(TG9_aa)
head(fxn_df_bin)
head(fxn_type)
TG9_fxn_list = data.frame(presence = fxn_df_bin[,grep("GCA_003288095", colnames(fxn_df_bin))],
                          cat = fxn_type$category,
                          fxn_name = rownames(fxn_type),
                          fxn_count = fxn_type$fxns_by_genome)

## inspect functions
TG9_fxn_list = TG9_fxn_list[order(TG9_fxn_list$fxn_count, decreasing = T),]
TG9_fxn_list = TG9_fxn_list[order(TG9_fxn_list$presence, decreasing = T),]
head(TG9_fxn_list)
dim(TG9_fxn_list) # total unique annotated gene products in Rhodococcus pangenome

## remove hypothetical protein from TG9 aa file (cannot screen function)
TG9_aa_trim = TG9_aa[-c(grep("hypothetical", getAnnot(TG9_aa)))]
head(TG9_aa_trim)
length(TG9_aa_trim) # TG9 gene products with annotated function

## TG9 gene sequences with gene products not specifically characterized in the pangenome (followed with manual search)
TG9_aa_unchar = TG9_aa_trim[c(121, 147, 160, 184, 188, 203, 212, 304, 307, 368, 
                              378, 503, 515, 563, 605, 686, 688, 791, 794, 803,
                              853, 867, 883, 938, 984, 999, 1009, 1048, 1051, 1073,
                              1077, 1111, 1209, 1232, 1300, 1317, 1340, 1354, 1389, 1403,
                              1404, 1439, 1463, 1514, 1548, 1550, 1576, 1603, 1613, 1616,
                              1666, 1686, 1704, 1888, 1904, 1910, 1944, 1957, 1998, 2001,
                              2023, 2064, 2090, 2097, 2098, 2158, 2226, 2344, 2345, 2354,
                              2373, 2405, 2448, 2481, 2495, 2583)]

## trim out uncharacterized seqs
TG9_aa_trim = TG9_aa_trim[-c(121, 147, 160, 184, 188, 203, 212, 304, 307, 368, 
                             378, 503, 515, 563, 605, 686, 688, 791, 794, 803,
                             853, 867, 883, 938, 984, 999, 1009, 1048, 1051, 1073,
                             1077, 1111, 1209, 1232, 1300, 1317, 1340, 1354, 1389, 1403,
                             1404, 1439, 1463, 1514, 1548, 1550, 1576, 1603, 1613, 1616,
                             1666, 1686, 1704, 1888, 1904, 1910, 1944, 1957, 1998, 2001,
                             2023, 2064, 2090, 2097, 2098, 2158, 2226, 2344, 2345, 2354,
                             2373, 2405, 2448, 2481, 2495, 2583)]

## assign pangenome fraction
fxn_aa_cat = c(1:length(getAnnot(TG9_aa_trim)))
for (i in 1:length(getAnnot(TG9_aa_trim))) {
  fxn_aa_cat[i] = as.character(TG9_fxn_list[grep(substr(unlist(getAnnot(TG9_aa_trim)[i]), start = 18, 
                                                        stop = nchar(unlist(getAnnot(TG9_aa_trim)[i]))),
                                                 TG9_fxn_list$fxn_name, fixed = T),]$cat)
}
fxn_aa_cat
length(TG9_aa_trim)
head(TG9_aa_trim)

## write out fasta files
# cloud
#write.fasta(getSequence(TG9_aa_trim)[which(fxn_aa_cat == "Cloud")], getAnnot(TG9_aa_trim)[which(fxn_aa_cat == "Cloud")], 'TG9_aa_seqs/TG9_aa_cloud.fa')
# shell
#write.fasta(getSequence(TG9_aa_trim)[which(fxn_aa_cat == "Shell")], getAnnot(TG9_aa_trim)[which(fxn_aa_cat == "Shell")], 'TG9_aa_seqs/TG9_aa_shell.fa')
# core
#write.fasta(getSequence(TG9_aa_trim)[which(fxn_aa_cat == "Core")], getAnnot(TG9_aa_trim)[which(fxn_aa_cat == "Core")], 'TG9_aa_seqs/TG9_aa_core.fa')
# uncharacterized
#write.fasta(getSequence(TG9_aa_unchar), getAnnot(TG9_aa_unchar), 'ad_hoc/aa_seqs/TG9_aa_unchar.fa')
# all, with "hypothetical protein" removed
#write.fasta(getSequence(TG9_aa_trim), getAnnot(TG9_aa_trim), 'TG9_aa_seqs/TG9_aa_no_hp.fa')
## NOTE: manually modify the above aa seq files with additional uncharacterized genes that receive an assignment after syntax correction and use for "blast koala" inputs

# totals
unlist(getAnnot(TG9_aa_trim)[which(fxn_aa_cat == "Cloud")]) # 108 + 2 (from manual search) = 110 "cloud" gene sequences
unlist(getAnnot(TG9_aa_trim)[which(fxn_aa_cat == "Shell")]) # 591 + 14 (from manual search) = 604 "shell" gene sequences
unlist(getAnnot(TG9_aa_trim)[which(fxn_aa_cat == "Core")]) # 1838 + 15 (from manual search) = 1,853 "core" gene sequences
unlist(getAnnot(TG9_aa_unchar)) # 76 "uncharacterized"; note that the syntax had to be exact match (31 are assigned as cloud/shell/core after manual search)


#####
### link cog categories to TG9 genome pangenome fractions
#####

## load fasta files
# fasta file -- cloud
TG9_aa_cloud <- read.fasta(file = "TG9_aa_seqs/pan_fractions/blast_koala_inputs/TG9_aa_cloud_edit.fa")
length(TG9_aa_cloud)
# fasta file -- shell
TG9_aa_shell <- read.fasta(file = "TG9_aa_seqs/pan_fractions/blast_koala_inputs/TG9_aa_shell_edit.fa")
length(TG9_aa_shell)
# fasta file -- accessory
TG9_aa_accessory <- read.fasta(file = "TG9_aa_seqs/pan_fractions/blast_koala_inputs/TG9_aa_accessory.fa") # concatenated cloud and shell fasta files
length(TG9_aa_accessory)
# fasta file -- core
TG9_aa_core <- read.fasta(file = "TG9_aa_seqs/pan_fractions/blast_koala_inputs/TG9_aa_core_edit.fa")
length(TG9_aa_core)
# fasta file -- no hp
TG9_aa_nohp <- read.fasta(file = "TG9_aa_seqs/pan_fractions/blast_koala_inputs/TG9_aa_no_hp.fa")
length(TG9_aa_nohp)

##  subset data from cog table
cog_cloud = TG9_cog[grep(paste(substr(unlist(getAnnot(TG9_aa_cloud)), start = 3, stop = 16), collapse = "|"), rownames(TG9_cog)),]
cog_shell = TG9_cog[grep(paste(substr(unlist(getAnnot(TG9_aa_shell)), start = 3, stop = 16), collapse = "|"), rownames(TG9_cog)),]
cog_accessory = TG9_cog[grep(paste(substr(unlist(getAnnot(TG9_aa_accessory)), start = 3, stop = 16), collapse = "|"), rownames(TG9_cog)),]
cog_core = TG9_cog[grep(paste(substr(unlist(getAnnot(TG9_aa_core)), start = 3, stop = 16), collapse = "|"), rownames(TG9_cog)),]

## counts table
cog_TG9_pan = data.frame(count_core = tapply(cog_core$COG_1, cog_core$COG_1, length),
                         count_accessory = tapply(cog_accessory$COG_1, cog_accessory$COG_1, length),
                         count_shell = tapply(cog_shell$COG_1, cog_shell$COG_1, length),
                         count_cloud = tapply(cog_cloud$COG_1, cog_cloud$COG_1, length))
cog_TG9_pan[is.na(cog_TG9_pan)] = 0

## manual adjust by multiple categories: visualize corrections needed
cog_core[which(nchar(as.character(cog_core$COG_2)) > 0),]
cog_accessory[which(nchar(as.character(cog_accessory$COG_2)) > 0),]
cog_shell[which(nchar(as.character(cog_shell$COG_2)) > 0),]
cog_cloud[which(nchar(as.character(cog_cloud$COG_2)) > 0),]

# modify accessory/shell set: account for multiple assignments by removal
tapply(cog_accessory[which(nchar(as.character(cog_accessory$COG_2)) > 0),]$COG_1,
       cog_accessory[which(nchar(as.character(cog_accessory$COG_2)) > 0),]$COG_1, length)
cog_TG9_pan$count_accessory[grep("C", rownames(cog_TG9_pan))] = c(80)
cog_TG9_pan$count_accessory[grep("F", rownames(cog_TG9_pan))] = c(15.5)
cog_TG9_pan$count_accessory[grep("I", rownames(cog_TG9_pan))] = c(43)
cog_TG9_pan$count_shell[grep("C", rownames(cog_TG9_pan))] = c(75)
cog_TG9_pan$count_shell[grep("F", rownames(cog_TG9_pan))] = c(14.5)
cog_TG9_pan$count_shell[grep("I", rownames(cog_TG9_pan))] = c(42)

# modify accessory/shell set: account for multiple assignments by addition
tapply(cog_accessory[which(nchar(as.character(cog_accessory$COG_2)) > 0),]$COG_2,
       cog_accessory[which(nchar(as.character(cog_accessory$COG_2)) > 0),]$COG_2, length)
cog_TG9_pan$count_accessory[grep("H", rownames(cog_TG9_pan))] = c(16.5)
cog_TG9_pan$count_accessory[grep("Q", rownames(cog_TG9_pan))] = c(57)
cog_TG9_pan$count_shell[grep("H", rownames(cog_TG9_pan))] = c(15.5)
cog_TG9_pan$count_shell[grep("Q", rownames(cog_TG9_pan))] = c(55)

# modify core set: account for multiple assignments by removal
c_set = tapply(cog_core[which(nchar(as.character(cog_core$COG_2)) > 0),]$COG_1,
               cog_core[which(nchar(as.character(cog_core$COG_2)) > 0),]$COG_1, length)
c_set[is.na(c_set)] = 0
cog_TG9_pan$count_core = cog_TG9_pan$count_core - (c_set/2)

# modify core set: account for multiple assignments by addition
tapply(cog_core[which(nchar(as.character(cog_core$COG_2)) > 0),]$COG_2,
       cog_core[which(nchar(as.character(cog_core$COG_2)) > 0),]$COG_2, length)
cog_TG9_pan$count_core[grep("G", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("G", rownames(cog_TG9_pan))]+3)
cog_TG9_pan$count_core[grep("H", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("H", rownames(cog_TG9_pan))]+0.5)
cog_TG9_pan$count_core[grep("J", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("J", rownames(cog_TG9_pan))]+1)
cog_TG9_pan$count_core[grep("M", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("M", rownames(cog_TG9_pan))]+0.5)
cog_TG9_pan$count_core[grep("O", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("O", rownames(cog_TG9_pan))]+0.5)
cog_TG9_pan$count_core[grep("P", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("P", rownames(cog_TG9_pan))]+0.5) 
cog_TG9_pan$count_core[grep("Q", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("Q", rownames(cog_TG9_pan))]+3.5)
cog_TG9_pan$count_core[grep("T", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("T", rownames(cog_TG9_pan))]+1.5)
cog_TG9_pan$count_core[grep("W", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("W", rownames(cog_TG9_pan))]+0.5)
# add 0.5 more for "P" since cog_3
cog_TG9_pan$count_core[grep("P", rownames(cog_TG9_pan))] = c(cog_TG9_pan$count_core[grep("P", rownames(cog_TG9_pan))]+0.5) 

## set fractions
cog_TG9_pan_frac = data.frame(cog_TG9_pan,
                              frac_core = 100*cog_TG9_pan$count_core/1851,
                              frac_accessory = 100*cog_TG9_pan$count_accessory/714,
                              frac_shell = 100*cog_TG9_pan$count_shell/604,
                              frac_cloud = 100*cog_TG9_pan$count_cloud/110)
rownames(cog_TG9_pan_frac)[1] = c("R")
cog_TG9_pan_frac$letter = rownames(cog_TG9_pan_frac)
cog_TG9_pan_frac$name = c("General function prediction only", "RNA processing and modification", "Chromatin structure and dynamics",
                          "Energy production and conversion", "Cell cycle control, cell division, chromosome partitioning",
                          "Amino acid transport and metabolism", "Nucleotide transport and metabolism", "Carbohydrate transport and metabolism",
                          "Coenzyme transport and metabolism", "Lipid transport and metabolism","Translation, ribosomal structure and biogenesis",
                          "Transcription", "Replication, recombination and repair", "Cell wall/membrane/envelope biogenesis",
                          "Cell motility", "Posttranslational modification, protein turnover, chaperones", "Inorganic ion transport and metabolism",
                          "Secondary metabolites biosynthesis, transport and catabolism", "Function unknown",
                          "Signal transduction mechanisms", "Intracellular trafficking, secretion, and vesicular transport", "Defense mechanisms",
                          "Cytoskeleton")
cog_TG9_pan_frac$cat = c("P", 
                         "I", "I", "M", "C", "M", "M", "M", "M", "M",
                         "I", "I", "I", "C", "C", "C", "M",
                         "M", "P", "C", "C", "C", "C")

## set category proportions
cat_prop = data.frame(core = tapply(cog_TG9_pan_frac$frac_core, cog_TG9_pan_frac$cat, sum),
                      acc = tapply(cog_TG9_pan_frac$frac_accessory, cog_TG9_pan_frac$cat, sum),
                      shell = tapply(cog_TG9_pan_frac$frac_shell, cog_TG9_pan_frac$cat, sum),
                      cloud = tapply(cog_TG9_pan_frac$frac_cloud, cog_TG9_pan_frac$cat, sum),
                      cat = c("Cellular Processes and Signaling", "Information Storage and Processing",
                              "Metabolism", "Poorly Characterized"))


### plot pie fractions (FIGURE 4B)
# core
t_cog_core = ggplot(cat_prop, aes(x="", y=core, fill = cat)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0, direction=-1) +
  scale_fill_manual(values = c("yellow", "skyblue", "violet", "gray")) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.border = element_blank(), 
        panel.grid=element_blank(), axis.ticks = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=22, color = "black")) +
  theme(legend.position = "none")
# shell
t_cog_shell = ggplot(cat_prop, aes(x="", y=shell, fill = cat)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0, direction=-1) +
  scale_fill_manual(values = c("yellow", "skyblue", "violet", "gray")) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.border = element_blank(), 
        panel.grid=element_blank(), axis.ticks = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=22, color = "black")) +
  theme(legend.position = "none")
# cloud
t_cog_cloud = ggplot(cat_prop, aes(x="", y=cloud, fill = cat)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0, direction=-1) +
  scale_fill_manual(values = c("yellow", "skyblue", "violet", "gray")) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.border = element_blank(), 
        panel.grid=element_blank(), axis.ticks = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12, color = "black")) +
  #theme(legend.position = "bottom") +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "none")

# plot 
grid.arrange(t_cog_cloud, t_cog_shell, t_cog_core, nrow = 1)


### plot bars by cog category (FIGURE 4C)
# set data
cog_pan_bars = data.frame(count = c(cog_TG9_pan_frac$count_core, cog_TG9_pan_frac$count_shell, cog_TG9_pan_frac$count_cloud),
                          percent = c(cog_TG9_pan_frac$frac_core,  cog_TG9_pan_frac$frac_shell, cog_TG9_pan_frac$frac_cloud),
                          p_fraction = c(rep("core", 23), rep("shell", 23), rep("cloud", 23)),
                          name = rep(cog_TG9_pan_frac$name, 3),
                          letter = rep(cog_TG9_pan_frac$letter, 3),
                          cat = rep(cog_TG9_pan_frac$cat, 3))
cog_pan_bars$p_fraction = factor(cog_pan_bars$p_fraction,
                                 levels = c("cloud", "shell", "core")) 
# plot
ggplot(cog_pan_bars, aes(x = letter, y = percent/100, fill = p_fraction)) +
  geom_bar(position=position_dodge(width=0.9), stat = "identity") +
  facet_grid(.~cat, scale="free_x", space="free") +
  ylab("TG9 genome component proportion") +
  xlab("COG Category") +
  theme_bw() +
  scale_fill_manual(values = c("red", "orange", "blue")) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")


### TG9 fractions from pangenome (FIGURE 4A)
tg9_pie = data.frame(count = apply(cog_TG9_pan, 2, sum),
                     cat = c("Core", "Accessory", "Shell", "Cloud"))
tg9_pie$cat = factor(tg9_pie$cat,
                     levels = c("Cloud", "Shell", "Core")) 
# plot
ggplot(tg9_pie[-c(2),], aes(x="", y=count, fill = cat)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() + 
  scale_fill_manual(values = c("red", "orange", "blue")) +
  ylab("COG-annotated genes") +
  theme(axis.text.x = element_text(size=14, color = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(size=14, color = "black"), 
        axis.title = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size=14, color = "black"), 
        legend.title = element_blank())
tapply(GPG_stat$count, GPG_stat$strain, sum)

