setwd("/Users/gaoyang/Documents/COMMENbT/Rscripts/experiment/PCA图/")
library(openxlsx)
install.packages("svglite")
install.packages(c("umap", "Rtsne", "cluster"))

library(vegan)
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(ggsignif)
library(svglite)
library(umap)
library(cluster)
library(Rtsne)

otu_raw <- read.xlsx("demoData1_rpkg.xlsx", rowNames = TRUE)
otu <- t(otu_raw)
head(otu)

#Calculate weighted Bray-Curtis distance
otu.distance <- vegdist(otu)
#pcoa
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
#Data format conversion and integration
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
#Simple plotting
p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p

#Add group information
group <- read.xlsx("group_mp_hs_NPSW.xlsx", rowNames = FALSE)
colnames(group) <- c("samples","group")
df <- merge(pc12,group,by="samples")
color=c("#1597A5","#FFC24B","#FEB3AE")


#####################################################################################

pdf("pcoa_DSR_NPSW_HS2.pdf", height = 11, width = 11)
ggplot(data=df,aes(x=V1,y=V2,
                       color=group,shape=group))+
  # Draw confidence ellipses:
  stat_ellipse(aes(fill = group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+
  theme_bw()+
  geom_point(size=3)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
  #guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PCoA1 ",pc[1],"%"),
       y=paste0("PCoA2 ",pc[2],"%"))+
  #scale_color_manual(values = color) +
  #scale_fill_manual(values = c("#1597A5","#FFC24B","#FEB3AE"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 14),axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),legend.title = element_text(size = 17),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),'cm')) 
#  theme(axis.title.x=element_text(size=12),
#        axis.title.y=element_text(size=12,angle=90),
#        axis.text.y=element_text(size=10),
#        axis.text.x=element_text(size=10),
#        panel.grid=element_blank())

dev.off()
#ggsave("pcoa_DSR_NPSW_HS1.pdf",height = 5,width = 7)

####################################################################################
############## Use PCoA and NMDS (beta diversity) to compare plasmid and phage mediated ARGs matrix differences ##############
# NMDS
# NMDS analysis
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_water.xlsx", rowNames = TRUE) 
otu <- t(otu_raw)
nmds1 <- metaMDS(otu, distance = 'bray', k = 2)
summary(nmds1)
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
group <- group[match(rownames(otu), rownames(group)),]
# Extract data
nmds1.stress <- nmds1$stress
nmds1.point <- data.frame(nmds1$point)
nmds1.species <- data.frame(nmds1$species)
sample_site <- nmds1.point[1:2]
sample_site$names <- rownames(sample_site)
colnames(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
# Combine group data
sample_site <- cbind(sample_site, group)

# Distance matrix
distance_matrix <- vegdist(otu, method = "bray")

# PERMANOVA analysis
adonis_result <- adonis(distance_matrix ~ Abu_method, data = group)
adonis_p_value <- adonis_result$aov.tab$`Pr(>F)`[1]

# SIMPER analysis
simper_result <- simper(distance_matrix, group$Abu_method)
simper_summary <- summary(simper_result)

# Mantel Test
original_distance_matrix <- vegdist(otu, method = "bray")
nmds_distance_matrix <- vegdist(nmds1$points)
mantel_result <- mantel(original_distance_matrix, nmds_distance_matrix)
mantel_r <- mantel_result$statistic
mantel_p_value <- mantel_result$signif

# Set plot boundaries
x_limits <- range(sample_site$NMDS1) + c(-0.1, 0.1)
y_limits <- range(sample_site$NMDS2) + c(-0.1, 0.1)
svglite("NMDS_plot_plasmid_phage_ARG_class_water.svg", width = 10, height = 8)
# Plot and add statistical information
ggplot(data=sample_site, aes(x=NMDS1, y=NMDS2, colour=Abu_method)) +
  geom_line(aes(group=Group), colour='Gray75', size=0.75) +
  geom_point(size=4, position="identity", alpha = 0.8) +
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title="NMDS_wood_MGE_ARG_class") + xlab("NMDS1") + ylab("NMDS2") +
  xlim(x_limits) + ylim(y_limits) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # Adjust plot margins
  annotate("text", x = min(x_limits), y = max(y_limits), 
           label = paste("Stress =", round(nmds1.stress, 3)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.05, 
           label = paste("PERMANOVA p =", formatC(adonis_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.1, 
           label = paste("Mantel r =", round(mantel_r, 3), "p =", formatC(mantel_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1)
dev.off()

# Plot and add statistical information and point names
ggplot(data=sample_site, aes(x=NMDS1, y=NMDS2, colour=Abu_method)) +
  geom_line(aes(group=Group), colour='Gray75', size=0.75) +
  geom_point(size=4, position="identity", alpha = 0.8) +
  geom_text(aes(label=names), size=3, vjust=1.5) +  # Add point names
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title="NMDS_wood_MGE_ARG_class") + xlab("NMDS1") + ylab("NMDS2") +
  xlim(x_limits) + ylim(y_limits) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # Adjust plot margins
  annotate("text", x = min(x_limits), y = max(y_limits), 
           label = paste("Stress =", round(nmds1.stress, 3)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.05, 
           label = paste("PERMANOVA p =", formatC(adonis_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.1, 
           label = paste("Mantel r =", round(mantel_r, 3), "p =", formatC(mantel_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1)


# PCoA
# PCoA analysis
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_nonbio.xlsx", rowNames = TRUE)
otu <- t(otu_raw)
# Calculate Bray-Curtis distance
otu.distance <- vegdist(otu)
pcoa <- cmdscale(otu.distance, eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100, digits=2)
# Group information
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
group <- group[match(rownames(otu), rownames(group)),]
# Data transformation and integration
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
sample_site <- pc12[1:2]
sample_site$names <- rownames(sample_site)
colnames(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
# Combine group data
sample_site <- cbind(sample_site, group)

# PERMANOVA analysis
adonis_result <- adonis(otu.distance ~ Abu_method, data = group)
adonis_p_value <- adonis_result$aov.tab$`Pr(>F)`[1]

# Set plot boundaries
x_limits <- range(sample_site$PCoA1) + c(-0.2, 0.2)
y_limits <- range(sample_site$PCoA2) + c(-0.2, 0.2)

# Save as SVG file
svglite("PCoA_plot_plasmid_phage_ARG_class_All_nonbio.svg", width = 10, height = 8)
ggplot(data = sample_site, aes(x = PCoA1, y = PCoA2, colour = Abu_method)) +
  geom_line(aes(group = Group), colour = 'Gray75', size = 0.75) +
  geom_point(size = 4, position = "identity", alpha = 0.8) +
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title = "PCOA_wood_MGE_ARG_class") +
  xlab(paste0("PCoA1 (", pc[1], "%)")) + 
  ylab(paste0("PCoA2 (", pc[2], "%)")) +
  xlim(x_limits) + ylim(y_limits) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # Adjust plot margins
  annotate("text", x = min(x_limits), y = max(y_limits), 
           label = paste("PERMANOVA p =", formatC(adonis_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1)
dev.off()

#UMAP
# read data
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_water.xlsx", rowNames = TRUE)
otu <- t(otu_raw)

# Group information
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
group <- group[match(rownames(otu), rownames(group)),]

# UMAP analysis
umap_result <- umap(otu)
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c('UMAP1', 'UMAP2')
umap_coords$names <- rownames(umap_coords)
sample_site_umap <- cbind(umap_coords, group)

# PERMANOVA analysis
adonis_result_umap <- adonis(vegdist(otu) ~ Abu_method, data = group)
adonis_p_value_umap <- adonis_result_umap$aov.tab$`Pr(>F)`[1]

# Silhouette analysis
silhouette_umap <- silhouette(as.numeric(factor(group$Abu_method)), dist(umap_coords))
avg_silhouette_width_umap <- mean(silhouette_umap[, 'sil_width'])

# Mantel Test
umap_distance_matrix <- dist(umap_coords)
mantel_result_umap <- mantel(vegdist(otu), umap_distance_matrix)
mantel_r_umap <- mantel_result_umap$statistic
mantel_p_value_umap <- mantel_result_umap$signif

# Set plot boundaries
x_limits_umap <- range(sample_site_umap$UMAP1) + c(-0.2, 0.2)
y_limits_umap <- range(sample_site_umap$UMAP2) + c(-0.2, 0.2)

# save SVG 
svglite("UMAP_plot_plasmid_phage_ARG_class_water.svg", width = 10, height = 8)
ggplot(data = sample_site_umap, aes(x = UMAP1, y = UMAP2, colour = Abu_method)) +
  geom_line(aes(group = Group), colour = 'Gray75', size = 0.75) +
  geom_point(size = 4, position = "identity", alpha = 0.8) +
  #geom_text(aes(label = names), size = 3, vjust = 1.5) +  # 添加点的名字
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title = "UMAP_wood_MGE_ARG_class") + xlab("UMAP1") + ylab("UMAP2") +
  xlim(x_limits_umap) + ylim(y_limits_umap) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # 调整图形边距
  annotate("text", x = max(x_limits_umap), y = min(y_limits_umap), 
           label = paste("PERMANOVA p =", formatC(adonis_p_value_umap, format = "e", digits = 2)), hjust = 1, vjust = 0) +
  annotate("text", x = max(x_limits_umap), y = min(y_limits_umap) + 0.05, 
           label = paste("Mantel r =", round(mantel_r_umap, 3), "p =", formatC(mantel_p_value_umap, format = "e", digits = 2)), hjust = 1, vjust = 0) +
  annotate("text", x = max(x_limits_umap), y = min(y_limits_umap) + 0.1, 
           label = paste("Avg Silhouette Width =", round(avg_silhouette_width_umap, 3)), hjust = 1, vjust = 0)
dev.off()


#tsne analysis
# Load data
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_bio.xlsx", rowNames = TRUE)
otu <- t(otu_raw)

# Group information
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
group <- group[match(rownames(otu), rownames(group)),]

# t-SNE analysis
set.seed(42)  # Seed for reproducibility
sample_count <- nrow(otu)
perplexity_value <- min(30, floor(sample_count / 3))
while (perplexity_value >= sample_count / 2) {
  perplexity_value <- perplexity_value - 1
}

if (perplexity_value < 5) {
  perplexity_value <- 5
  warning("Perplexity value is too small. Setting to minimum value of 5.")
}

if (perplexity_value >= sample_count / 4) {
  stop("Perplexity value is too large for the number of samples. Please check your sample size.")
}

tsne_result <- Rtsne(otu, dims = 2, perplexity = perplexity_value)

# Extract t-SNE coordinates
tsne_coords <- as.data.frame(tsne_result$Y)
colnames(tsne_coords) <- c('tSNE1', 'tSNE2')
tsne_coords$names <- rownames(tsne_coords)
sample_site_tsne <- cbind(tsne_coords, group)

# PERMANOVA analysis
adonis_result_tsne <- adonis(vegdist(otu) ~ Abu_method, data = group)
adonis_p_value_tsne <- adonis_result_tsne$aov.tab$`Pr(>F)`[1]

# Silhouette analysis
silhouette_tsne <- silhouette(as.numeric(factor(group$Abu_method)), dist(tsne_coords))
avg_silhouette_width_tsne <- mean(silhouette_tsne[, 'sil_width'])

# Mantel Test
tsne_distance_matrix <- dist(tsne_coords)
mantel_result_tsne <- mantel(vegdist(otu), tsne_distance_matrix)
mantel_r_tsne <- mantel_result_tsne$statistic
mantel_p_value_tsne <- mantel_result_tsne$signif

# Set plot boundaries
x_limits_tsne <- range(sample_site_tsne$tSNE1) + c(-0.2, 0.2)
y_limits_tsne <- range(sample_site_tsne$tSNE2) + c(-0.2, 0.2)

# Save as SVG file
svglite("tSNE_plot_plasmid_phage_ARG_class_MPs.svg", width = 10, height = 8)
ggplot(data = sample_site_tsne, aes(x = tSNE1, y = tSNE2, colour = Abu_method)) +
  geom_line(aes(group = Group), colour = 'Gray75', size = 0.75) +
  geom_point(size = 4, position = "identity", alpha = 0.8) +
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title = "tSNE_wood_MGE_ARG_class") + xlab("tSNE1") + ylab("tSNE2") +
  xlim(x_limits_tsne) + ylim(y_limits_tsne) +
  theme(plot.margin = margin(10, 50, 10, 10)) +
  annotate("text", x = min(x_limits_tsne), y = max(y_limits_tsne), 
           label = paste("PERMANOVA p =", formatC(adonis_p_value_tsne, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits_tsne), y = max(y_limits_tsne) - 0.05, 
           label = paste("Mantel r =", round(mantel_r_tsne, 3), "p =", formatC(mantel_p_value_tsne, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits_tsne), y = max(y_limits_tsne) - 0.1, 
           label = paste("Avg Silhouette Width =", round(avg_silhouette_width_tsne, 3)), hjust = 0, vjust = 1)
dev.off()

####################################################################################
##############To compare plasmid and phage-mediated ARGs using PCoA and NMDS###############

# Load necessary libraries
library(RColorBrewer)
library(vegan)
library(openxlsx)

# Load data
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
otu <- t(otu_raw)

# Calculate Shannon index with base e
shannon_e <- diversity(otu, index = 'shannon', base = exp(1))

# Calculate Shannon index with base 2
shannon_2 <- diversity(otu, index = 'shannon', base = 2)

# Output Shannon index results
shannon_df <- data.frame(Shannon_e = shannon_e, Shannon_2 = shannon_2)
print(shannon_df)

# Calculate species count
species_count <- specnumber(otu)

# Calculate Pielou’s Evenness Index
pielou_evenness <- shannon_e / log(species_count)
pielou_evenness_df <- data.frame(Pielou_Evenness = pielou_evenness)
print(pielou_evenness_df)

# Calculate Gini-Simpson index
gini_simpson <- diversity(otu, index = 'simpson')
gini_simpson_df <- data.frame(Gini_Simpson = gini_simpson)
print(gini_simpson_df)

# Classic Simpson index
classic_simpson <- 1 - gini_simpson
classic_simpson_df <- data.frame(Classic_Simpson = classic_simpson)
print(classic_simpson_df)
simpson_index <- 1 - Gini_simpson
simpson_index <- as.data.frame(simpson_index)
simpson_index

#Read plot data
Plot.Adiversity <- read.table(file = "Figure4_metadata_and_organized_alpha_diversity.txt",header = T,sep="\t", row.names=1)
Plot.Deviation <- read.table("Deviation_alpha.txt",header = T, row.names = 1,sep="\t")
Plot.Adiversity <-read.xlsx("demoData2_plasmid_phage_ARG_class_alpha_shannon_总表.xlsx", rowNames = TRUE)
#Plot
svglite("/Users/gaoyang/Documents/微塑料抗性基因——毕业论文/2毕业论文/图片/alpha_beta_diversity/Alpha_diversity_Shannon_Simpson_plasmid_phage_ARG_class_compare2.svg", width = 10, height = 10)
ggplot(data=Plot.Adiversity, aes(x=MGE_type, y=value, fill=MGE_type))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_line(aes(group=Group), colour='Gray75', size=0.75)+
  geom_boxplot() + theme_bw()+
  facet_wrap(Method~Env, scales = "free", ncol = 4)+
  scale_fill_manual(values = brewer.pal(8, "Set1"))+
  geom_signif(comparisons = list(c("plasmid", "phage")),
              step_increase = 0.1, map_signif_level = T,
              test = wilcox.test)
dev.off()


