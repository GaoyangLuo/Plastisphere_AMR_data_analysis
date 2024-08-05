setwd("/Users/gaoyang/Documents/COMMENbT/Rscripts/experiment/散点图线性回归")
install.packages("compositions")
install.packages(c("coda", "MASS", "Matrix", "mvtnorm", "pracma"))
install.packages("remotes")
remotes::install_github("jimhester/bayesm")
install.packages("bayesm")
install.packages("philentropy") 

# 加载必要的库
library(readxl)
library(ggplot2)
library(reshape2)
library(compositions)
library(philentropy)
# 读取数据并转换为数值型矩阵
read_and_convert_to_numeric <- function(file_path, sheet) {
  data <- read_excel(file_path, sheet = sheet)
  # 去掉第一列（类别列）
  data <- data[, -1]
  # 查看数据类型
  print(sapply(data, class))
  # 将所有数据转换为数值型
  data <- data %>% mutate(across(everything(), as.numeric))
  # 再次查看数据类型
  print(sapply(data, class))
  return(as.matrix(data))
}

# 读取数据
matrix1 <- read_and_convert_to_numeric("DSR_ARG_abundance_cellcopy.xlsx", sheet = 1)
matrix2 <- read_and_convert_to_numeric("DSR_ARG_abundance_cellkk2.xlsx", sheet = 1)
matrix3 <- read_and_convert_to_numeric("DSR_ARG_abundance_rpkg.xlsx", sheet = 1)
matrix4 <- read_and_convert_to_numeric("DSR_ARG_abundance_rpkm.xlsx", sheet = 1)

# 定义计算相似度的函数
# 定义计算相似度的函数（皮尔逊相关系数）
calculate_similarity <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  return(cor(vec1, vec2, method = "pearson"))
}

calculate_cosine_similarity <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  return(sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2))))
}

calculate_euclidean_distance <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  return(sqrt(sum((vec1 - vec2)^2)))
}

calculate_manhattan_distance <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  return(sum(abs(vec1 - vec2)))
}

calculate_jaccard_similarity <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  intersection <- sum(vec1 & vec2)
  union <- sum(vec1 | vec2)
  return(intersection / union)
}

calculate_bray_curtis_dissimilarity <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  return(vegdist(rbind(vec1, vec2), method = "bray")[1])
}

# 定义计算 Jensen-Shannon divergence (JSD) 的函数
calculate_jsd <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  # 转换为概率分布，确保每个向量的元素之和为 1
  vec1 <- vec1 / sum(vec1)
  vec2 <- vec2 / sum(vec2)
  prob_matrix <- rbind(vec1, vec2)
  jsd_value <- distance(prob_matrix, method = "jensen-shannon")
  return(sqrt(jsd_value))
}

# 计算所有组合之间的相似度
similarities <- data.frame(
  pair = c("matrix1 vs matrix2", "matrix1 vs matrix3", "matrix1 vs matrix4",
           "matrix2 vs matrix3", "matrix2 vs matrix4", "matrix3 vs matrix4"),
  similarity = c(
    calculate_similarity(matrix1, matrix2),
    calculate_similarity(matrix1, matrix3),
    calculate_similarity(matrix1, matrix4),
    calculate_similarity(matrix2, matrix3),
    calculate_similarity(matrix2, matrix4),
    calculate_similarity(matrix3, matrix4)
  ),
  cosine_similarity = c(
    calculate_cosine_similarity(matrix1, matrix2),
    calculate_cosine_similarity(matrix1, matrix3),
    calculate_cosine_similarity(matrix1, matrix4),
    calculate_cosine_similarity(matrix2, matrix3),
    calculate_cosine_similarity(matrix2, matrix4),
    calculate_cosine_similarity(matrix3, matrix4)
  ),
  euclidean_distance = c(
    calculate_euclidean_distance(matrix1, matrix2),
    calculate_euclidean_distance(matrix1, matrix3),
    calculate_euclidean_distance(matrix1, matrix4),
    calculate_euclidean_distance(matrix2, matrix3),
    calculate_euclidean_distance(matrix2, matrix4),
    calculate_euclidean_distance(matrix3, matrix4)
  ),
  manhattan_distance = c(
    calculate_manhattan_distance(matrix1, matrix2),
    calculate_manhattan_distance(matrix1, matrix3),
    calculate_manhattan_distance(matrix1, matrix4),
    calculate_manhattan_distance(matrix2, matrix3),
    calculate_manhattan_distance(matrix2, matrix4),
    calculate_manhattan_distance(matrix3, matrix4)
  ),
  jaccard_similarity = c(
    calculate_jaccard_similarity(matrix1, matrix2),
    calculate_jaccard_similarity(matrix1, matrix3),
    calculate_jaccard_similarity(matrix1, matrix4),
    calculate_jaccard_similarity(matrix2, matrix3),
    calculate_jaccard_similarity(matrix2, matrix4),
    calculate_jaccard_similarity(matrix3, matrix4)
  ),
  bray_curtis_dissimilarity = c(
    calculate_bray_curtis_dissimilarity(matrix1, matrix2),
    calculate_bray_curtis_dissimilarity(matrix1, matrix3),
    calculate_bray_curtis_dissimilarity(matrix1, matrix4),
    calculate_bray_curtis_dissimilarity(matrix2, matrix3),
    calculate_bray_curtis_dissimilarity(matrix2, matrix4),
    calculate_bray_curtis_dissimilarity(matrix3, matrix4)
  ),
  jsd_sqrt = c(
    calculate_jsd(matrix1, matrix2),
    calculate_jsd(matrix1, matrix3),
    calculate_jsd(matrix1, matrix4),
    calculate_jsd(matrix2, matrix3),
    calculate_jsd(matrix2, matrix4),
    calculate_jsd(matrix3, matrix4)
  )
)
# 打印相似度数据框
print(similarities)

# 创建数据框用于绘制散点图
df <- data.frame(
  value1 = c(as.vector(matrix1), as.vector(matrix1), as.vector(matrix1), 
             as.vector(matrix2), as.vector(matrix2), as.vector(matrix3)),
  value2 = c(as.vector(matrix2), as.vector(matrix3), as.vector(matrix4),
             as.vector(matrix3), as.vector(matrix4), as.vector(matrix4)),
  pair = c(
    rep("matrix1 vs matrix2", length(matrix1)),
    rep("matrix1 vs matrix3", length(matrix1)),
    rep("matrix1 vs matrix4", length(matrix1)),
    rep("matrix2 vs matrix3", length(matrix2)),
    rep("matrix2 vs matrix4", length(matrix2)),
    rep("matrix3 vs matrix4", length(matrix3))
  )
)


# 绘制散点图并添加相似度值和蓝色渐变色
svglite("/Users/gaoyang/Documents/COMMENbT/Rscripts/experiment/散点图线性回归/ARG_Abundance_calculation2.svg", width = 10, height = 6)
p <- ggplot(df, aes(x = value1, y = value2)) +
  geom_point(aes(color = value2), alpha = 0.5) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  facet_wrap(~pair, scales = "free") +
  labs(title = "Scatter Plots of ARG Abundance Calculations",
       x = "Value 1",
       y = "Value 2") +
  theme_minimal() +
  geom_text(data = similarities, aes(x = -Inf, y = Inf, label = sprintf("r = %.2f\ncosine = %.2f\neuclidean = %.2f\nmanhattan = %.2f\njaccard = %.2f\nbray_curtis = %.2f\njsd = %.2f",
                                                                        similarity, cosine_similarity, euclidean_distance, manhattan_distance, jaccard_similarity, bray_curtis_dissimilarity,jsd_sqrt)),
            hjust = -0.1, vjust = 1.1, size = 3, color = "red")

print(p)
dev.off()




# 计算所有组合之间的相似度
similarities <- data.frame(
  pair = c("matrix1 vs matrix2", "matrix1 vs matrix3", "matrix1 vs matrix4",
           "matrix2 vs matrix3", "matrix2 vs matrix4", "matrix3 vs matrix4"),
  similarity = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                           list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                             res <- calculate_similarity(m[[1]], m[[2]])
                             return(res$cor)
                           }),
  similarity_p = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                             list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                               res <- calculate_similarity(m[[1]], m[[2]])
                               return(res$p_value)
                             }),
  cosine_similarity = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                                  list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                                    res <- calculate_cosine_similarity(m[[1]], m[[2]])
                                    return(res$cosine)
                                  }),
  cosine_p = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                         list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                           res <- calculate_cosine_similarity(m[[1]], m[[2]])
                           return(res$p_value)
                         }),
  euclidean_distance = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                                   list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                                     return(calculate_euclidean_distance(m[[1]], m[[2]]))
                                   }),
  manhattan_distance = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                                   list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                                     return(calculate_manhattan_distance(m[[1]], m[[2]]))
                                   }),
  jaccard_similarity = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                                   list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                                     return(calculate_jaccard_similarity(m[[1]], m[[2]]))
                                   }),
  bray_curtis_dissimilarity = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                                          list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                                            return(calculate_bray_curtis_dissimilarity(m[[1]], m[[2]]))
                                          }),
  jsd = sapply(list(list(matrix1, matrix2), list(matrix1, matrix3), list(matrix1, matrix4), 
                    list(matrix2, matrix3), list(matrix2, matrix4), list(matrix3, matrix4)), function(m) {
                      return(calculate_jsd(m[[1]], m[[2]]))
                    })
)

# 转换数据框为长格式
plot_data <- similarities %>%
  gather(key = "metric", value = "value", -pair, -similarity_p, -cosine_p)

# 为每个统计值绘制箱图
metrics <- c("similarity", "cosine_similarity", "euclidean_distance", 
             "manhattan_distance", "jaccard_similarity", 
             "bray_curtis_dissimilarity", "jsd")

for (metric in metrics) {
  print(
    ggplot(plot_data %>% filter(metric == !!metric), aes(x = pair, y = value, color = pair)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Boxplot of", metric, "between matrices"),
           x = "Matrix Pair", y = metric)
  )
}

for (metric in metrics) {
  print(
        ggplot(plot_data, aes(x = pair, y = value, fill = pair)) +
          geom_bar(stat = "identity", position = "dodge") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "Bar plots of different metrics between matrices",
               x = "Matrix Pair", y = "Value") +
          facet_wrap(~ metric, scales = "free_y")
  )
}
