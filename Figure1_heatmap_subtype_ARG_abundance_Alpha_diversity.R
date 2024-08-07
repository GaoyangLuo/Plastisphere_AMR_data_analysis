if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("circlize")

install.packages("reshape2")
install.packages("tibble")
install.packages("tidyr")
install.packages("data.table")
if (!requireNamespace("ggsignif", quietly = TRUE)) install.packages("ggsignif")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")

library(pheatmap)
library(readr)
library(dplyr)
library(readxl)
library(reshape2)
library(tibble)
library(tidyr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)


# read ARG abundance metadata
arg_abundance <- fread("DSR_merged_samples_with_class_fillZero_multidrgu_adjust.tsv", header = TRUE, sep = "\t")

# read group metadata
group_info <- read_excel("group_顺序.xlsx")

# use melt func change data to long type
arg_abundance_long <- melt(arg_abundance, id.vars = c("ARG_name", "Class", "Database", "MGE_type"), 
                           variable.name = "Sample", value.name = "Abundance")

# cat group
arg_abundance_long <- arg_abundance_long %>%
  left_join(group_info, by = c("Sample" = "Sample"))

# calculate ARG total abundance
arg_total_abundance <- arg_abundance_long %>%
  group_by(ARG_name) %>%
  summarise(Total_Abundance = sum(Abundance))

# choose top 50 abundant ARG
top_n <- 50
top_args <- arg_total_abundance %>%
  top_n(n = top_n, wt = Total_Abundance) %>%
  pull(ARG_name)

# 过滤原始数据以仅包含高丰度的ARG filter high abundant ARG
arg_abundance_filtered <- arg_abundance_long %>%
  filter(ARG_name %in% top_args)

# 计算每个ARG的相对丰度 calculate relative ARG abundance
arg_relative_abundance <- arg_abundance_filtered %>%
  group_by(ARG_name) %>%
  summarise(Relative_Abundance = sum(Abundance) / sum(arg_abundance_filtered$Abundance))

 
# 提取分组信息并按 Class 排序 sort according to Class order
row_annotation <- arg_abundance_filtered %>%
  filter(ARG_name %in% top_args) %>%
  select(ARG_name, Class, MGE_type) %>%
  distinct() %>%
  left_join(arg_relative_abundance, by = "ARG_name") %>%
  arrange(Class) %>%
  column_to_rownames(var = "ARG_name")

# 按Class排序后的ARG名称顺序 sorted col of ARGs names
sorted_args <- rownames(row_annotation)

# 转换为矩阵格式 change to matrice
arg_matrix <- arg_abundance_long %>%
  select(ARG_name, Sample, Abundance) %>%
  spread(key = Sample, value = Abundance) %>%
  column_to_rownames(var = "ARG_name")

# 按丰度顺序选择前50个ARG choose top 50 abundant ARGs
arg_matrix <- arg_matrix[sorted_args, , drop = FALSE]

# 确保列顺序与 group_info 中的顺序一致 ensure the same sort as that in group_info
sample_order <- group_info$Sample
arg_matrix <- arg_matrix[, sample_order, drop = FALSE]

# 计算每个样本的总抗性基因丰度 calculate total ARG abundance of each sample
sample_total_abundance <- colSums(arg_matrix)

# 将所有值转换为数值类型并处理无效值 process NA and change all value to as.matrix
arg_matrix <- as.matrix(arg_matrix)
arg_matrix[is.na(arg_matrix)] <- 0       # 将NA值替换为0 change NA to 0
arg_matrix[is.nan(arg_matrix)] <- 0      # 将NaN值替换为0 change NaN to 0
arg_matrix[is.infinite(arg_matrix)] <- 0 # 将Inf值替换为0 change Inf to 0

# 检查是否还有无效值 check there is invalid values
invalid_values <- which(is.na(arg_matrix) | is.nan(arg_matrix) | is.infinite(arg_matrix), arr.ind = TRUE)

if (nrow(invalid_values) > 0) {
  stop("Matrix contains invalid values at positions: ", paste(apply(invalid_values, 1, paste, collapse = ", "), collapse = "; "))
}

# 提取分组信息 acquire catagory information
#row_annotation <- arg_abundance_long %>%
#  select(ARG_name, Class, MGE_type) %>%
#  distinct() %>%
#  column_to_rownames(var = "ARG_name")

# 提取列注释信息 acquire col information
col_annotation <- group_info %>%
  column_to_rownames(var = "Sample")

# 使用RColorBrewer生成Set3调色板 use RColorBrewer and Set3
set3_colors <- brewer.pal(12, "Set3")

# 函数：生成颜色映射，允许颜色重复 map color and alow repeated colour
generate_colors <- function(categories) {
  n <- length(unique(categories))
  colors <- rep(set3_colors, length.out = n)
  names(colors) <- unique(categories)
  return(colors)
}

# 创建颜色映射 create mapped colour
group_colors <- generate_colors(group_info$Group)
class_colors <- generate_colors(row_annotation$Class)
mge_type_colors <- generate_colors(row_annotation$MGE_type)

# 创建 ComplexHeatmap 所需的注释对象 create subject for annotation with ComplexHeatmap
# 创建包含条形图和分组信息的顶部分注释 create top annotation 
top_annotation <- HeatmapAnnotation(
  Total_Abundance = anno_barplot(sample_total_abundance,
                                 which = "column",
                                 height = unit(4, "cm")),
  Group = group_info$Group,
  col = list(Group = group_colors),
  annotation_name_side = "left",
  annotation_legend_param = list(Group = list(title_gp = gpar(fontsize = 20,fontfamily = "sans"), 
                                              labels_gp = gpar(fontsize = 20,fontfamily = "sans")))
)

# 创建包含条形图的左侧行注释 create left barplot annotation
left_barplot_annotation <- rowAnnotation(
  Relative_Abundance = anno_barplot(
    row_annotation$Relative_Abundance, 
    which = "row", 
    gp = gpar(fill = "grey"),  # 设置柱子颜色为灰色),  set color of coloum
    bar_width = 0.8, 
    border = FALSE, 
    axis_param = list(direction = "reverse"),
    width = unit(4, "cm"),
  )
)

# 创建包含分组信息的左侧行注释 create left group annotation
left_group_annotation <- rowAnnotation(
  MGE_type = row_annotation$MGE_type,
  Class = row_annotation$Class,
  col = list(Group = group_colors),
  annotation_legend_param = list(
    MGE_type = list(title_gp = gpar(fontsize = 20), 
                    labels_gp = gpar(fontsize = 20)),
    Class = list(title_gp = gpar(fontsize = 20), 
                 labels_gp = gpar(fontsize = 20))
  )
)

# 合并左侧行注释 group left annotation and barplot
left_annotation <- c(left_barplot_annotation, left_group_annotation)


# 检查矩阵的数据类型和内容 check data type and content
str(arg_matrix)
head(arg_matrix)

# 设置颜色渐变 set colour
pdf("heatmap_ARG_subtype_withAbundance_multidrug_adjust.pdf", height = 40, width = 40)   #8,9
pdf("heatmap_ARG_subtype_withoutAbundance_multidrug_adjust_1.pdf", height = 40, width = 40)   #8,9
color_scheme <- colorRampPalette(c("#FDEBEA", "#D5281F"))(100)
# 绘制热图 draw heatmap
Heatmap(
  arg_matrix,
  name = "Abundance",
  col = color_scheme,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 18,fontfamily = "sans"),
  column_names_gp = gpar(fontsize = 18,fontfamily = "sans"),
  rect_gp = gpar(col = "black", lwd = 1),  # 设置格子线颜色和宽度
  #cell_fun = function(j, i, x, y, width, height, fill) {
   # grid.text(sprintf("%.1e", arg_matrix[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
  #},
  top_annotation = top_annotation,  # 顶部分注释
  left_annotation = left_annotation,  # 左侧行注释
  heatmap_legend_param = list(
    title = "Abundance",
    title_gp = gpar(fontsize = 15),fontfamily = "sans",  # 设置热图图例标题的字体大小
    legend_gp = gpar(fontsize = 20,fontfamily = "sans")  # 设置图例字体大小
  ),
  width = unit(10 * ncol(arg_matrix), "mm"),
  height = unit(8 * nrow(arg_matrix), "mm"),
  column_title = "Heatmap",
  column_title_gp = gpar(fontsize = 20, fontface = "bold",fontfamily = "sans")
)


dev.off()

#默认设置，不画上左柱状图
pdf("heatmap_ARG_subtype_test3.pdf", height = 40, width = 20)   #8,9
pheatmap(arg_matrix,
         # 是否聚类：
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         # 加color bar：列注释信息；
         annotation_row = row_annotation,
         annotation_col = col_annotation,
         top_annotation = top_annotation,
         left_annotation = left_annotation,
         # color bar 颜色设定：
         annotation_colors = ann_colors,
         # 设置单元格颜色渐变；(100)表示分100段渐变；
         color = colorRampPalette(c("#FDEBEA","#D5281F"))(100), 
         # 行、列标签的字体大小 8,10
         fontsize_col = 10,
         fontsize_row = 10,
         # 设置每个单元格的宽度和高度 15,20
         cellwidth = 15, 
         cellheight = 15,
         # 行、列聚类树的高度：
         treeheight_row = 50, 
         treeheight_col = 30,
         #display_numbers = TRUE参数设定在每个热图格子中显示相应的数值，
         # number_color参数设置数值字体的颜色
         #display_numbers = arg_matrix,
         number_format = "%.1e",
         fontsize_number = 8,
         number_color = "black",
         # 设置标题：
         main = "Heatmap",
         heatmap_legend_param = list(title = "Abundance"),
         name = "Abundance")

dev.off()

# 绘制热图
pheatmap(arg_matrix,
         annotation_row = row_annotation,
         annotation_col = col_annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6)
###############################计算Alpha diversity###########################

# 检查arg_matrix的结构
print(dim(arg_matrix))
print(colnames(arg_matrix))

# 转置矩阵以确保样本为行，ARG为列
arg_matrix_for_alpha_diversity <- t(arg_matrix)

# 确保数据是整数
arg_matrix_for_alpha_diversity_int <- apply(arg_matrix_for_alpha_diversity, 2, as.integer)

# 计算每个样本的Shannon指数
shannon_index <- diversity(arg_matrix_for_alpha_diversity, index = "shannon")

# 计算每个样本的Simpson指数
simpson_index <- diversity(arg_matrix_for_alpha_diversity, index = "simpson")

# 计算每个样本的Chao1指数
chao1_result <- estimateR(arg_matrix_for_alpha_diversity_int)
chao1_index <- chao1_result["S.chao1", ]

# 计算Pielou's均匀度指数
pielou_evenness <- shannon_index / log(specnumber(arg_matrix_for_alpha_diversity))

# 创建一个数据框来存储结果
alpha_diversity <- data.frame(
  Sample = rownames(arg_matrix_for_alpha_diversity),
  Shannon = shannon_index,
  Simpson = simpson_index,
  Chao1 = chao1_index,
  Pielou = pielou_evenness
)

# 合并分组信息
alpha_diversity <- alpha_diversity %>%
  left_join(group_info, by = c("Sample" = "Sample"))

# 仅保留感兴趣的组别
alpha_diversity <- alpha_diversity %>%
  filter(Group %in% c("BD-MPs", "nBD-MPs"))

# 输出结果
print(alpha_diversity)

# 创建函数来生成带有显著性标记的箱线图
create_boxplot_with_significance <- function(index, y_label) {
  groups <- unique(alpha_diversity$Group)
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  # 计算每个比较的显著性
  p_value_list <- sapply(comparisons, function(comp) {
    group_data <- alpha_diversity %>% filter(Group %in% comp)
    t.test(group_data[[index]] ~ group_data$Group)$p.value
  })
  
  # 根据 p 值添加显著性标记
  #annotations <- ifelse(p_value_list < 0.05, "*", "ns")
  
  # 为显著的比较添加星号标记
  significance_levels <- ifelse(p_value_list < 0.001, "***",
                                ifelse(p_value_list < 0.01, "**",
                                       ifelse(p_value_list < 0.05, "*", "ns")))
  # 只保留显著的比较和其对应的星号标记
  significant_comparisons <- comparisons[p_value_list < 0.05]
  significant_significance_levels <- significance_levels[p_value_list < 0.05]
  
  # 获取每组数据的最大值，并稍微向上调整以确保空间足够放置星号
  max_y <- max(alpha_diversity[[index]], na.rm = TRUE)
  increment <- (max_y - min(alpha_diversity[[index]])) * 0.05
  y_positions <- seq(max_y + increment, by = increment, length.out = length(significant_comparisons))
  
  # 创建箱线图并添加显著性标记
  p <- ggplot(alpha_diversity, aes(x = Group, y = get(index), fill = Group)) +
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = significant_comparisons, 
                          annotations = significant_significance_levels,
                          map_signif_level = TRUE,
                          y_position = y_positions) +
    theme_minimal() +
    labs(title = paste(index, "Index by Group"), y = y_label, x = "Group") +
    theme(legend.position = "none")
  
  return(p)
}

# 创建各指数的箱线图
p_shannon <- create_boxplot_with_significance("Shannon", "Shannon Index")
p_simpson <- create_boxplot_with_significance("Simpson", "Simpson Index")
p_chao1 <- create_boxplot_with_significance("Chao1", "Chao1 Index")
p_pielou <- create_boxplot_with_significance("Pielou", "Pielou's Evenness")

# 将多个图表合并成一个图
combined_plot <- ggpubr::ggarrange(p_shannon, p_simpson, p_chao1, p_pielou, 
                                   ncol = 2, nrow = 2, 
                                   common.legend = TRUE, legend = "bottom")

# 打印合并后的图表
print(combined_plot)
ggsave("alpha_diversity_plots_shannon_simpson_chao1_pielou.pdf", plot = combined_plot, device = cairo_pdf, width = 10, height = 8)


# 统计分析：比较两组的 Shannon 指数
shannon_t_test <- t.test(Shannon ~ Group, data = alpha_diversity)
print(shannon_t_test)

# 统计分析：比较两组的 Simpson 指数
simpson_t_test <- t.test(Simpson ~ Group, data = alpha_diversity)
print(simpson_t_test)

# 绘制 Shannon 指数的箱线图
p1 <- ggplot(alpha_diversity, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Index by Group", y = "Shannon Index", x = "Group")
print(p1)
# 保存 Shannon 指数图为 PDF
# 保存合并后的图表为可编辑的 PDF


  


# 绘制Shannon和Simpson指数的箱线图
p1 <- ggplot(alpha_diversity, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Index by Group", y = "Shannon Index", x = "Group")

print(p1)
ggsave("shannon_index.pdf", plot = p1, width = 8, height = 6)

p2 <- ggplot(alpha_diversity, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Simpson Index by Group", y = "Simpson Index", x = "Group")

print(p2)

# 绘制Chao1指数的箱线图
p3 <- ggplot(alpha_diversity, aes(x = Group, y = Chao1, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Chao1 Index by Group", y = "Chao1 Index", x = "Group")

print(p3)

# 绘制Pielou’s均匀度指数的箱线图
p4 <- ggplot(alpha_diversity, aes(x = Group, y = Pielou, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Pielou's Evenness by Group", y = "Pielou's Evenness", x = "Group")

print(p4)

print(alpha_diversity)