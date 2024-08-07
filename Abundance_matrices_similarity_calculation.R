install.packages("compositions")
install.packages(c("coda", "MASS", "Matrix", "mvtnorm", "pracma"))
install.packages("remotes")
remotes::install_github("jimhester/bayesm")
install.packages("bayesm")
install.packages("philentropy") 

# Load necessary libraries
library(readxl)
library(ggplot2)
library(reshape2)
library(compositions)
library(philentropy)
library(dplyr)

# Read data and convert to numeric matrix
read_and_convert_to_numeric <- function(file_path, sheet) {
  data <- read_excel(file_path, sheet = sheet)
  # Remove the first column (category column)
  data <- data[, -1]
  # Check data types
  print(sapply(data, class))
  # Convert all data to numeric
  data <- data %>% mutate(across(everything(), as.numeric))
  # Check data types again
  print(sapply(data, class))
  return(as.matrix(data))
}

# Read data
matrix1 <- read_and_convert_to_numeric("DSR_ARG_abundance_cellcopy.xlsx", sheet = 1)
matrix2 <- read_and_convert_to_numeric("DSR_ARG_abundance_cellkk2.xlsx", sheet = 1)
matrix3 <- read_and_convert_to_numeric("DSR_ARG_abundance_rpkg.xlsx", sheet = 1)
matrix4 <- read_and_convert_to_numeric("DSR_ARG_abundance_rpkm.xlsx", sheet = 1)

# Define a function to calculate similarity (Pearson correlation coefficient)
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

# Define a function to calculate Jensen-Shannon divergence (JSD)
calculate_jsd <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  # Convert to probability distributions, ensuring each vector sums to 1
  vec1 <- vec1 / sum(vec1)
  vec2 <- vec2 / sum(vec2)
  prob_matrix <- rbind(vec1, vec2)
  jsd_value <- distance(prob_matrix, method = "jensen-shannon")
  return(sqrt(jsd_value))
}

# Calculate similarity for all combinations
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
# Print similarity dataframe
print(similarities)

# Create dataframe for scatter plot
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


# Plot scatter plot with similarity values and blue gradient
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




# Calculate similarity for all combinations
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

# Convert dataframe to long format
plot_data <- similarities %>%
  gather(key = "metric", value = "value", -pair, -similarity_p, -cosine_p)

# Convert dataframe to long format
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
