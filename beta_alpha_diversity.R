setwd("/Users/gaoyang/Documents/COMMENbT/Rscripts/experiment/PCA图/")
library(openxlsx)
install.packages("svglite")
install.packages(c("umap", "Rtsne", "cluster"))
#加载包
library(vegan)#计算距离时需要的包
library(ggplot2)#绘图包
library(ggprism)
library(RColorBrewer)
library(ggsignif)
library(svglite)
library(umap)
library(cluster)
library(Rtsne)
#Dr.Jiang
#otu_raw <- read.table(file="DrXiaotao/oral_argsoap.normalize_cellnumber.type.tab.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
#otu_raw <- read.table(file="otu.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
otu_raw <- read.xlsx("demoData1_rpkg.xlsx", rowNames = TRUE)
otu <- t(otu_raw)
head(otu)

#计算加权bray_curtis距离
otu.distance <- vegdist(otu)
#pcoa分析
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
#数据格式转换及数据整合
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
#简单绘图
p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p

#添加分组信息
group <- read.xlsx("group_mp_hs_NPSW.xlsx", rowNames = FALSE)
colnames(group) <- c("samples","group")
df <- merge(pc12,group,by="samples")
color=c("#1597A5","#FFC24B","#FEB3AE")

#####################################################################################
#参考
ggplot(data = df,aes(x=V1,y=V2,
                     color=group,shape=group))+
  # 绘制置信椭圆：
  stat_ellipse(aes(fill = group),#dsr_arg_input$type
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  # 绘制散点：
  geom_point(size = 3.5)+
  labs(x = df$V1,y = df$V2,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  #scale_fill_manual(values = c("purple","orange","pink"))+
  #scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm')) +coord_fixed(ratio = 3)
#####################################################################################

pdf("pcoa_DSR_NPSW_HS2.pdf", height = 11, width = 11)
ggplot(data=df,aes(x=V1,y=V2,
                       color=group,shape=group))+
  # 绘制置信椭圆：
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

#anosim
#数据处理及PCoA分析
#otu_raw <- read.table(file="otu.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
otu_raw <- read.xlsx("demoData1.xlsx", rowNames = TRUE)
otu <- t(otu_raw)
otu.distance <- vegdist(otu)
PCoA <- cmdscale (otu.distance,eig=TRUE)
pc12 <- PCoA$points[,1:2]
pc <- round(PCoA$eig/sum(PCoA$eig)*100,digits=2)#解释度
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
#group <- read.table("group.txt", sep='\t', header=T)
group <- read.xlsx("demoData2.xlsx", rowNames = FALSE)
colnames(group) <- c("samples","group")
df <- merge(pc12,group,by="samples")
head(df)
#绘图
ggplot(df,aes(x=V1, y=V2,color=group,shape=group))+#指定数据、X轴、Y轴
  geom_point(size=3)+
  theme_bw()


#使用vegan包中的anosim函数进行anosim分析
group <- read.xlsx("demoData2_noWood.xlsx", rowNames = TRUE)
df_anosim <- anosim(otu,df$group,permutations = 999)#数据也可以是原始otu数据
df_anosim 
#df_anosim <- anosim(otu,df$group,permutations = 999)
#整理出作图数据
df1<-data.frame(
  x=df_anosim$class.vec,
  y=df_anosim$dis.rank
)
df1
#绘图 
par(family='Arial') #'STKaiti'
ggplot(df1,aes(x=x,y=y))+
  stat_boxplot(geom = "errorbar", width=0.15,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=x), 
               outlier.colour="white",size=0.5)+
  theme_bw()+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  #scale_fill_manual(values = brewer.pal(8, "Set1"))+
  scale_fill_manual(values=c("#1597A5","red","#FEB3AE","#FFC24B","gray")) + #指定颜色
  ggtitle("Bray-Curtis Anosim")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 14,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+
  theme(legend.position = 'none')+
  theme(text=element_text(family='Arial'))+
  geom_signif(comparisons = list(c("biodegradable_plastic", "non-biodegradable_plastic")),
              step_increase = 0.1, map_signif_level = T,
              test = wilcox.test)+
  labs(x = paste("R=",df_anosim$statistic,", ","p=", df_anosim$signif),
       y = "Rank of Distance (Bray_Curtis)")

####################################### Anosim的解释与缺陷 #####################################
#Anderson and Walsh (2013) considered this issue. I'll summarise their findings with respect to your question.

#ANOSIM is very sensitive to heterogeneity and the results of Anderson & Walsh would suggest that don't trust the ANOSIM results; 
#they'll basically just tell you that there is some difference 
#(be it in terms of location (differences in mean), 
#dispersion (variances) or correlation structure), 
#not that there is a location difference were a significant ANOSIM result be obtained.
#PERMANOVA (which is basically adonis()) was found to be largely unaffected by heterogeneity in Anderson & Walsh's simulations but only for balanced designs.

#For unbalanced designs PERMANOVA and ANOSIM were

#too liberal if the smaller group had greater dispersion, and
#too conservative if the larger group had greater dispersion.
#This result was especially so for ANOSIM.

#Basically, how much you can trust the results of your PERMANOVA depends on the balance in the design.

#Anderson MJ, Walsh DCI. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological monographs [Internet] 2013; 83: 557. Available from: http://doi.org/10.1890/12-2010.1
########################################################################################



#############Nat Method Code#############
#ggplot(df1,aes(x=x,y=y))+
#stat_boxplot(geom = "errorbar",width=0.15)+
  #geom_line(aes(group=SampleID), colour='Gray75', size=0.75)+
  #geom_boxplot() + theme_bw()+
  #facet_wrap(Method~Env, scales = "free", ncol = 5)+
  #scale_fill_manual(values = brewer.pal(8, "Set1"))+
  #geom_signif(comparisons = list(c("Taxonomic", "Sequence")),
              #step_increase = 0.1, map_signif_level = T,
              #test = wilcox.test)
################################################################


#MRPP分析
MRPP <- mrpp(otu.distance,df$group,permutations = 999)
MRPP

#Adonis
#分组信息
group <- read.xlsx("demoData2.xlsx", rowNames = TRUE)

# 基于bray-curtis距离进行计算
#下面的代码检验了microplastic type对ARG组成差异影响的显著程度，
#获得P-value=0.001 < 0.05，表示塑料种类方式对ARG组成有显著影响。
#data(dune)
#data(dune.env)
set.seed(1)
otu.div <- adonis2(otu ~ Group, data = group, permutations = 999, method="bray")
otu.div


#我们还需要利用betadisper评估下每组样本ARG组成的多元一致性 (Multivariate homogeneity of groups dispersions (variances))。
#如下代码计算出P=0.168表示不同分组样品检测指标的离散度(方差)没有显著差异。
#那么，adonis检测出的差异就是因为每组数据在空间的中心点不同造成的，进一步说明Management对物种组成有显著影响。
# 计算加权bray-curtis距离
otu.distance <- vegdist(otu, method="bray")
otu.betadisper <- adonis2(otu.distance~Group,data=group,distance = "bray",permutations = 999)
otu.betadisper

dispersion <- betadisper(otu.distance, group=group$Group)
permutest(dispersion)
#画个图看看
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
ggplot(dispersion,aes(x=V1, y=V2,color=group,shape=group))+#指定数据、X轴、Y轴
  geom_point(size=3)+
  theme_bw()

################以下是对adonis显示出p有显著性，但betadiverse分析中p值没有显著性的解释#################
#Q: When running adonis (vegan package) I got an r2 = 0.45, and p = 0.001. 
#When I ran the betadisper and ran a subsequent permutation test I got an F = 1 and p = 0.3.
#A: A non-significant result in betadisper is not necessarily related to a significant/non-significant result in adonis. 
#The two tests are testing different hypothesis. The former tests homogeneity of dispersion among groups (regions in your case), 
#which is a condition (assumption) for adonis. The latter tests no difference in ‘location’, 
#that is, tests whether composition among groups is similar or not. 
#You may have the centroids of two groups in NMS at a very similar position in the ordination space, 
#but if their dispersions are quite different, adonis will give you a significant p-value, 
#thus, the result is heavily influenced not by the difference in composition between groups
#but by differences in composition within groups (heterogeneous dispersion, and thus a measure of beta diversity). 
#In short, your results are fine, you are meeting the ‘one assumption’ for adonis (homogeneous dispersion) and 
#thus you are certain that results from adonis are ‘real’ and not an artifact of heterogeneous dispersions. 
#For more information you can read Anderson (2006) Biometrics 62(1):245-253 and Anderson (2006) Ecology Letters 9(6):683-693. Hope this helps!
#https://stats.stackexchange.com/questions/212137/betadisper-and-adonis-in-r-am-i-interpreting-my-output-correctly



####################################################################################
##############利用PCoA和NMDS(beta diversity)比较plasmid和phage介导的ARGs矩阵差异性##############
#NMDS
#nmds分析
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_water.xlsx", rowNames = TRUE) #"demoData1_16s-rpkm_bio.xlsx"
otu <- t(otu_raw)
nmds1 <- metaMDS(otu, distance = 'bray', k = 2)
summary(nmds1)
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE) #"demoData2_noWood_noWater.xlsx"
group <- group[match(rownames(otu),rownames(group)),]
#提取数据
nmds1.stress <- nmds1$stress
nmds1.point <- data.frame(nmds1$point)
nmds1.species <- data.frame(nmds1$species)
sample_site <- nmds1.point[1:2]
sample_site$names <- rownames(sample_site)
colnames(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
#合并分组数据
sample_site <- cbind(sample_site,group)

#距离矩阵
distance_matrix <- vegdist(otu, method = "bray")

# PERMANOVA 分析
adonis_result <- adonis(distance_matrix ~ Abu_method, data = group)
adonis_p_value <- adonis_result$aov.tab$`Pr(>F)`[1]

# SIMPER 分析
simper_result <- simper(distance_matrix, group$Abu_method)
simper_summary <- summary(simper_result)

# Mantel Test
original_distance_matrix <- vegdist(otu, method = "bray")
nmds_distance_matrix <- vegdist(nmds1$points)
mantel_result <- mantel(original_distance_matrix, nmds_distance_matrix)
mantel_r <- mantel_result$statistic
mantel_p_value <- mantel_result$signif

# 设置图形边界
x_limits <- range(sample_site$NMDS1) + c(-0.1, 0.1)
y_limits <- range(sample_site$NMDS2) + c(-0.1, 0.1)
#pdf("NMDS_plot_plasmid_phage_ARG_class.svg", width = 10, height = 8)   #8,9
svglite("NMDS_plot_plasmid_phage_ARG_class_water.svg", width = 10, height = 8)
# 绘图并添加统计学信息
ggplot(data=sample_site, aes(x=NMDS1, y=NMDS2, colour=Abu_method)) +
  geom_line(aes(group=Group), colour='Gray75', size=0.75) +
  geom_point(size=4, position="identity", alpha = 0.8) +
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title="NMDS_wood_MGE_ARG_class") + xlab("NMDS1") + ylab("NMDS2") +
  xlim(x_limits) + ylim(y_limits) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # 调整图形边距
  annotate("text", x = min(x_limits), y = max(y_limits), 
           label = paste("Stress =", round(nmds1.stress, 3)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.05, 
           label = paste("PERMANOVA p =", formatC(adonis_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.1, 
           label = paste("Mantel r =", round(mantel_r, 3), "p =", formatC(mantel_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1)
dev.off()

# 绘图并添加统计学信息和点的名字
ggplot(data=sample_site, aes(x=NMDS1, y=NMDS2, colour=Abu_method)) +
  geom_line(aes(group=Group), colour='Gray75', size=0.75) +
  geom_point(size=4, position="identity", alpha = 0.8) +
  geom_text(aes(label=names), size=3, vjust=1.5) +  # 添加点的名字
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title="NMDS_wood_MGE_ARG_class") + xlab("NMDS1") + ylab("NMDS2") +
  xlim(x_limits) + ylim(y_limits) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # 调整图形边距
  annotate("text", x = min(x_limits), y = max(y_limits), 
           label = paste("Stress =", round(nmds1.stress, 3)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.05, 
           label = paste("PERMANOVA p =", formatC(adonis_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits), y = max(y_limits) - 0.1, 
           label = paste("Mantel r =", round(mantel_r, 3), "p =", formatC(mantel_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1)

#绘图
#P1<-ggplot(data=sample_site,aes(x=NMDS1,y=NMDS2,colour=sample_site$group))+theme_bw()+theme(panel.grid= element_line(color =NA),
#                                                               panel.grid.minor = element_line(color = NA),
#                                                              panel.border = element_rect(fill = NA, colour ="black"),
#                                                               axis.text.x  = element_text(size=15,family="A", colour="black",hjust = 0.7),
#                                                               axis.title.x = element_text(vjust=0.2, size = 15,family="A"),
#                                                               axis.text.y  = element_text(size=15,family="A", colour="black",hjust = 0.7),
#                                                               axis.title.y = element_text(vjust=1, size = 15,family="A", colour="black"))+
#  geom_point(aes(color=Biodiversity),size=2,alpha=0.9)+theme(legend.position= "top")+theme(panel.grid=element_blank())+
#  #下面的#00FF00等分别是颜色的16进制代码，可自己百度，需要注意的是，点的颜色用scale_color_manual（）函数，置信椭圆用scale_fill_manual（）
#  scale_color_manual(values=c(R1 = "#00FF00", R2 = "#F74ED6", R4 = "#AD07E3"))+
#  stat_ellipse(data=sample_site,geom = "polygon",aes(fill=Biodiversity),alpha=0.35,level = 0.95)+
#  scale_fill_manual(values=c(R1 = "#00FF00", R2 = "#F74ED6", R4 = "#AD07E3"))
#P1

#原来的图
ggplot(data=sample_site, aes(x=NMDS1, y=NMDS2, colour=Abu_method))+
  geom_line(aes(group=Group), colour='Gray75', size=0.75)+
  geom_point(size=4, position="identity", alpha = 0.8)+
  theme_bw()+
  scale_colour_manual(values = brewer.pal(8, "Set1"))+
  labs(title="NMDS_wood_MGE_ARG_class") + xlab("NMDS1") + ylab("NMDS2")

#PCoA
#pcoa分析
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_nonbio.xlsx", rowNames = TRUE) #"demoData1_16s-rpkm_bio.xlsx"
otu <- t(otu_raw)
#计算bray_curtis距离
otu.distance <- vegdist(otu)
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
#分组信息
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
group <- group[match(rownames(otu),rownames(group)),]
#数据格式转换及数据整合
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
sample_site <- pc12[1:2]
sample_site$names <- rownames(sample_site)
colnames(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
#合并分组数据
sample_site <- cbind(sample_site,group)

# PERMANOVA 分析
adonis_result <- adonis(otu.distance ~ Abu_method, data = group)
adonis_p_value <- adonis_result$aov.tab$`Pr(>F)`[1]

# 设置图形边界
x_limits <- range(sample_site$PCoA1) + c(-0.2, 0.2)
y_limits <- range(sample_site$PCoA2) + c(-0.2, 0.2)

# 保存为 SVG 文件
svglite("PCoA_plot_plasmid_phage_ARG_class_All_nonbio.svg", width = 10, height = 8)
ggplot(data = sample_site, aes(x = PCoA1, y = PCoA2, colour = Abu_method)) +
  geom_line(aes(group = Group), colour = 'Gray75', size = 0.75) +
  geom_point(size = 4, position = "identity", alpha = 0.8) +
  #geom_text(aes(label = names), size = 3, vjust = 1.5) +  # 添加点的名字
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title = "PCOA_wood_MGE_ARG_class") + xlab(paste0("PCoA1 (", pc[1], "%)")) + ylab(paste0("PCoA2 (", pc[2], "%)")) +
  xlim(x_limits) + ylim(y_limits) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # 调整图形边距
  annotate("text", x = min(x_limits), y = max(y_limits), 
           label = paste("PERMANOVA p =", formatC(adonis_p_value, format = "e", digits = 2)), hjust = 0, vjust = 1)
dev.off()

#原来的图
ggplot(data=sample_site , aes(x=PCoA1, y=PCoA2, colour=Abu_method))+
  geom_line(aes(group=Group), colour='Gray75', size=0.75)+
  geom_point(size=4, position="identity", alpha = 0.8)+
  theme_bw()+
  scale_colour_manual(values = brewer.pal(8, "Set1"))+
  labs(title="PCOA_wood_MGE_ARG_class") + xlab("PCOA1") + ylab("PCOA2")


#UMAP
# 读取数据
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_water.xlsx", rowNames = TRUE)
otu <- t(otu_raw)

# 分组信息
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
group <- group[match(rownames(otu), rownames(group)),]

# UMAP 分析
umap_result <- umap(otu)
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c('UMAP1', 'UMAP2')
umap_coords$names <- rownames(umap_coords)
sample_site_umap <- cbind(umap_coords, group)

# PERMANOVA 分析
adonis_result_umap <- adonis(vegdist(otu) ~ Abu_method, data = group)
adonis_p_value_umap <- adonis_result_umap$aov.tab$`Pr(>F)`[1]

# Silhouette 分析
silhouette_umap <- silhouette(as.numeric(factor(group$Abu_method)), dist(umap_coords))
avg_silhouette_width_umap <- mean(silhouette_umap[, 'sil_width'])

# Mantel Test
umap_distance_matrix <- dist(umap_coords)
mantel_result_umap <- mantel(vegdist(otu), umap_distance_matrix)
mantel_r_umap <- mantel_result_umap$statistic
mantel_p_value_umap <- mantel_result_umap$signif

# 设置图形边界
x_limits_umap <- range(sample_site_umap$UMAP1) + c(-0.2, 0.2)
y_limits_umap <- range(sample_site_umap$UMAP2) + c(-0.2, 0.2)

# 保存为 SVG 文件
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


#tsne分析
# 读取数据
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_bio.xlsx", rowNames = TRUE)
otu <- t(otu_raw)
# 分组信息
group <- read.xlsx("demoData2_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE)
group <- group[match(rownames(otu), rownames(group)),]
# t-SNE 分析
set.seed(42)  # 设置随机种子以确保结果可重复
sample_count <- nrow(otu)
# 确保 perplexity 小于样本数量的一半，并且是合理的值
perplexity_value <- min(30, floor(sample_count / 3))
while (perplexity_value >= sample_count / 2) {
  perplexity_value <- perplexity_value - 1
}

# 检查 perplexity 是否合理
if (perplexity_value < 5) {
  perplexity_value <- 5
  warning("Perplexity value is too small. Setting to minimum value of 5.")
}

# 检查 perplexity 是否小于样本数量的一半
if (perplexity_value >= sample_count / 4) {
  stop("Perplexity value is too large for the number of samples. Please check your sample size.")
}

tsne_result <- Rtsne(otu, dims = 2, perplexity = perplexity_value)
# 提取 t-SNE 坐标
tsne_coords <- as.data.frame(tsne_result$Y)
colnames(tsne_coords) <- c('tSNE1', 'tSNE2')
tsne_coords$names <- rownames(tsne_coords)
sample_site_tsne <- cbind(tsne_coords, group)

# PERMANOVA 分析
adonis_result_tsne <- adonis(vegdist(otu) ~ Abu_method, data = group)
adonis_p_value_tsne <- adonis_result_tsne$aov.tab$`Pr(>F)`[1]

# Silhouette 分析
silhouette_tsne <- silhouette(as.numeric(factor(group$Abu_method)), dist(tsne_coords))
avg_silhouette_width_tsne <- mean(silhouette_tsne[, 'sil_width'])

# Mantel Test
tsne_distance_matrix <- dist(tsne_coords)
mantel_result_tsne <- mantel(vegdist(otu), tsne_distance_matrix)
mantel_r_tsne <- mantel_result_tsne$statistic
mantel_p_value_tsne <- mantel_result_tsne$signif

# 设置图形边界
x_limits_tsne <- range(sample_site_tsne$tSNE1) + c(-0.2, 0.2)
y_limits_tsne <- range(sample_site_tsne$tSNE2) + c(-0.2, 0.2)

# 保存为 SVG 文件
svglite("tSNE_plot_plasmid_phage_ARG_class_MPs.svg", width = 10, height = 8)
ggplot(data = sample_site_tsne, aes(x = tSNE1, y = tSNE2, colour = Abu_method)) +
  geom_line(aes(group = Group), colour = 'Gray75', size = 0.75) +
  geom_point(size = 4, position = "identity", alpha = 0.8) +
  #geom_text(aes(label = names), size = 3, vjust = 1.5) +  # 添加点的名字
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8, "Set1")) +
  labs(title = "tSNE_wood_MGE_ARG_class") + xlab("tSNE1") + ylab("tSNE2") +
  xlim(x_limits_tsne) + ylim(y_limits_tsne) +
  theme(plot.margin = margin(10, 50, 10, 10)) +  # 调整图形边距
  annotate("text", x = min(x_limits_tsne), y = max(y_limits_tsne), 
           label = paste("PERMANOVA p =", formatC(adonis_p_value_tsne, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits_tsne), y = max(y_limits_tsne) - 0.05, 
           label = paste("Mantel r =", round(mantel_r_tsne, 3), "p =", formatC(mantel_p_value_tsne, format = "e", digits = 2)), hjust = 0, vjust = 1) +
  annotate("text", x = min(x_limits_tsne), y = max(y_limits_tsne) - 0.1, 
           label = paste("Avg Silhouette Width =", round(avg_silhouette_width_tsne, 3)), hjust = 0, vjust = 1)
dev.off()

####################################################################################
##############利用PCoA和NMDS(Alpha diversity)比较plasmid和phage介导的ARGs矩阵差异性###############

#Library R packages
library(RColorBrewer)

#loading data
otu_raw <- read.xlsx("demoData1_plasmid_phage_ARG_class_总表.xlsx", rowNames = TRUE) #"demoData1_16s-rpkm_bio.xlsx"
otu <- t(otu_raw)


#Shannon 指数,通常使用2、e作为底数
#以e作为底数表示方法
Shannon <- diversity(otu, index = 'shannon', base = exp(1))
#以2作为底数表示方法
Shannon <- diversity(otu, index = 'shannon', base = 2)
#输出Shannon_index结果
Shannon <- as.data.frame(Shannon)
Shannon
# 计算物种数量
species_count <- specnumber(otu)
# 计算 Pielou’s Evenness Index
pielou_evenness <- Shannon / log(species_count)
pielou_evenness <- as.data.frame(pielou_evenness)
#Simpson指数分为经典 Simpson 指数和Gini-Simpson 指数，不过平时常用的 Simpson 指数即为 Gini-Simpson 指数
#Gini-Simpson 指数代码
Gini_simpson  <- diversity(otu, index = 'simpson')
#经常使用
Gini_simpson
#经典 Simpson 指数
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
#ggsave("/Users/gaoyang/Documents/微塑料抗性基因——毕业论文/2毕业论文/图片/alpha_beta_diversity/Alpha_diversity_Shannon_Simpson_plasmid_phage_ARG_class_compare.pdf",height = 8,width = 7)

