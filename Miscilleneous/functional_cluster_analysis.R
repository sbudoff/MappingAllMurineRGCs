library(tidyverse)
library(rstatix)
library(FSA)


root <- '/home/sam/FinalRGC_xenium/'
# group_key <- read_csv(paste0(root,'human_derived_cluster-groups.csv'))

# Embedding derived groups
group_key <- read_csv('/home/sam/FinalRGC_xenium/embeddingModelFinal/20241221_k5_HclustPromising_latent64_splits32/clustering/hierarchical_clusters.csv') %>%
  rename(group = cluster,
         cluster = map_index) %>%
  mutate(cluster = cluster + 1,
         group = group + 1)


unification <- read_csv(paste0(root,'RGC Subtype Unification.csv')) %>%
  rename(cluster = Tran2019_Clusters,
         Spot_Size = `Spot_Size_(um)`) %>%
  left_join(group_key, by = 'cluster') %>%
  select(group, Tran2019_Names,Goetz2022_PatchSeq, Goetz2022_Name, Functional_Group, 
  sup_PeakResponseLatency, sup_OnOffIndex, Spot_Size,sup_PeakFiringRate, 
  sup_ResponseDuration, sup_BLFiringRate,sup_SuppressionIndex
  ) %>%
  drop_na() %>%
  mutate(group = factor(group),
         Functional_Group = factor(Functional_Group))


# List of numerical variables to test
numerical_vars <- c("sup_PeakResponseLatency", "sup_OnOffIndex", "Spot_Size",
                  "sup_PeakFiringRate", "sup_ResponseDuration", "sup_BLFiringRate",
                  "sup_SuppressionIndex")

library(ggplot2)
library(rstatix)
library(FSA)
library(dplyr)
library(ggpubr)

# Function to run tests and create visualizations
run_tests <- function(data, var) {
  # Print diagnostic information
  cat("\n\nChecking variable:", var, "\n")
  cat("Number of non-NA values:", sum(!is.na(data[[var]])), "\n")
  cat("Values by group:\n")
  print(table(data$group, is.na(data[[var]])))
  
  # Create plot
  p <- ggplot(data, aes(x = group, y = !!sym(var))) +
    geom_boxplot() +
    theme_classic() +
    ggtitle(var) +
    ylab(var) +
    xlab("Group")
  
  # Try to run Kruskal-Wallis test with error handling
  tryCatch({
    kw_result <- kruskal.test(formula(paste(var, "~ group")), data = data)
    
    # Check if p-value exists and is not NA
    if(!is.na(kw_result$p.value) && kw_result$p.value < 0.05) {
      dunn_result <- dunnTest(formula(paste(var, "~ group")), 
                              data = data,
                              method="bh")
      
      # Add significance brackets if any pairwise comparisons are significant
      sig_comparisons <- dunn_result$res %>%
        filter(P.adj < 0.05) %>%
        mutate(comparison = strsplit(Comparison, " - "))
      
      if(nrow(sig_comparisons) > 0) {
        p <- p + stat_compare_means(comparisons = sig_comparisons$comparison,
                                    method = "wilcox.test",
                                    label = "p.adj")
      }
      
      return(list(variable = var,
                  kw_pvalue = kw_result$p.value,
                  significant = TRUE,
                  dunn = dunn_result$res,
                  plot = p))
    } else {
      return(list(variable = var,
                  kw_pvalue = ifelse(is.na(kw_result$p.value), NA, kw_result$p.value),
                  significant = FALSE,
                  dunn = NULL,
                  plot = p))
    }
  }, error = function(e) {
    cat("Error in analysis:", conditionMessage(e), "\n")
    return(list(variable = var,
                kw_pvalue = NA,
                significant = FALSE,
                dunn = NULL,
                error = conditionMessage(e),
                plot = p))
  })
}

# Run tests and create plots
results <- lapply(numerical_vars, function(var) run_tests(unification, var))

# Print results and arrange plots
plots <- list()
for(res in results) {
  cat("\n\n=== Results for", res$variable, "===\n")
  if(!is.na(res$kw_pvalue)) {
    cat("Kruskal-Wallis p-value:", format.pval(res$kw_pvalue, digits = 3), "\n")
    
    if(res$significant) {
      cat("Significant differences found. Dunn's test results:\n")
      print(res$dunn)
      
      medians <- unification %>%
        group_by(group) %>%
        summarise(median = median(!!sym(res$variable), na.rm = TRUE)) %>%
        arrange(median)
      
      cat("\nMedians by group:\n")
      print(medians)
    } else {
      cat("No significant differences between groups\n")
    }
  } else {
    if(!is.null(res$error)) {
      cat("Analysis failed with error:", res$error, "\n")
    } else {
      cat("Analysis could not be completed (possibly due to insufficient data)\n")
    }
  }
  
  plots[[res$variable]] <- res$plot
}

# Arrange all plots in a grid
ggarrange(plotlist = plots, ncol = 3, nrow = 3)



#########################################################################


scaled_data <- unification %>%
  select(numerical_vars) %>%
  scale() %>%
  data.frame()%>%
  mutate(group = unification$group)

library(lsa)
library(ggpubr)

# Function to perform both similarity measures bootstrap
bootstrap_similarities <- function(data, numerical_vars, n_iterations = 10000) {
  groups <- unique(data$group[data$group != 1])
  results_cos <- data.frame(matrix(ncol = length(groups), nrow = n_iterations))
  results_euc <- data.frame(matrix(ncol = length(groups), nrow = n_iterations))
  colnames(results_cos) <- paste0("Group_", groups)
  colnames(results_euc) <- paste0("Group_", groups)
  
  group1_data <- data[data$group == 1, numerical_vars]
  
  for(i in 1:n_iterations) {
    group1_sample <- as.numeric(group1_data[sample(nrow(group1_data), 1), ])
    
    for(g in groups) {
      group_data <- data[data$group == g, numerical_vars]
      group_sample <- as.numeric(group_data[sample(nrow(group_data), 1), ])
      
      results_cos[i, paste0("Group_", g)] <- cosine(group1_sample, group_sample)
      results_euc[i, paste0("Group_", g)] <- sqrt(sum((group1_sample - group_sample)^2))
    }
  }
  
  return(list(cosine = results_cos, euclidean = results_euc))
}

# Run bootstrap
set.seed(18)
bootstrap_results <- bootstrap_similarities(scaled_data, numerical_vars)

# Analyze both metrics
analyze_metric <- function(results, metric_name) {
  results_long <- results %>%
    pivot_longer(cols = everything(),
                 names_to = "group",
                 values_to = "similarity")
  
  aov_result <- aov(similarity ~ group, data = results_long)
  tukey_result <- TukeyHSD(aov_result)
  
  # Dynamic significance plotting
  p_values <- tukey_result$group[, "p adj"]
  names(p_values) <- rownames(tukey_result$group)
  y_max <- max(results_long$similarity)
  spacing <- (max(results_long$similarity) - min(results_long$similarity)) * 0.1
  
  comparisons <- list()
  y_positions <- c()
  for(i in seq_along(p_values)) {
    pair <- strsplit(names(p_values)[i], "-")[[1]]
    comparisons[[i]] <- pair
    y_positions[i] <- y_max + (spacing * i)
  }
  
  plot <- ggplot(results_long, aes(x = group, y = similarity)) +
    geom_boxplot() +
    theme_classic() +
    labs(x = "Group", y = metric_name) +
    stat_compare_means(comparisons = comparisons,
                       label = "p.adj",
                       method = "t.test",
                       label.y = y_positions) +
    scale_y_continuous(limits = c(min(results_long$similarity), 
                                  max(y_positions) + spacing))
  
  return(plot)
}

# Create and combine plots
p1 <- analyze_metric(bootstrap_results$cosine, "Cosine Similarity to Group 1")
p2 <- analyze_metric(bootstrap_results$euclidean, "Euclidean Distance to Group 1")

ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))


#############################################################################################
library(tidyverse)
library(lsa)
library(ggpubr)
library(reshape2)
library(pheatmap)

# Function to perform pairwise bootstrap comparisons
bootstrap_pairwise_similarities <- function(data, numerical_vars, n_iterations = 10000) {
  groups <- unique(data$group)
  n_groups <- length(groups)
  
  # Initialize matrices to store results
  cos_matrix <- matrix(NA, nrow = n_groups, ncol = n_groups)
  euc_matrix <- matrix(NA, nrow = n_groups, ncol = n_groups)
  rownames(cos_matrix) <- colnames(cos_matrix) <- paste0("Group_", groups)
  rownames(euc_matrix) <- colnames(euc_matrix) <- paste0("Group_", groups)
  
  # Lists to store all bootstrap results
  all_cos_results <- list()
  all_euc_results <- list()
  all_plots <- list()
  
  for(base_group in groups) {
    compare_groups <- groups[groups != base_group]
    results_cos <- data.frame(matrix(ncol = length(compare_groups), nrow = n_iterations))
    results_euc <- data.frame(matrix(ncol = length(compare_groups), nrow = n_iterations))
    colnames(results_cos) <- paste0("Group_", compare_groups)
    colnames(results_euc) <- paste0("Group_", compare_groups)
    
    base_group_data <- data[data$group == base_group, numerical_vars]
    
    for(i in 1:n_iterations) {
      base_sample <- as.numeric(base_group_data[sample(nrow(base_group_data), 1), ])
      
      for(g in compare_groups) {
        group_data <- data[data$group == g, numerical_vars]
        group_sample <- as.numeric(group_data[sample(nrow(group_data), 1), ])
        
        results_cos[i, paste0("Group_", g)] <- cosine(base_sample, group_sample)
        results_euc[i, paste0("Group_", g)] <- sqrt(sum((base_sample - group_sample)^2))
      }
    }
    
    # Store results
    all_cos_results[[paste0("Group_", base_group)]] <- results_cos
    all_euc_results[[paste0("Group_", base_group)]] <- results_euc
    
    # Calculate means for matrices
    for(g in compare_groups) {
      cos_matrix[paste0("Group_", base_group), paste0("Group_", g)] <- 
        mean(results_cos[[paste0("Group_", g)]], na.rm = TRUE)
      euc_matrix[paste0("Group_", base_group), paste0("Group_", g)] <- 
        mean(results_euc[[paste0("Group_", g)]], na.rm = TRUE)
    }
    
    # Create plots for this base group
    plots <- analyze_metric_pair(results_cos, results_euc, 
                                 paste("Comparisons with Group", base_group))
    all_plots[[paste0("Group_", base_group)]] <- plots
  }
  
  # Create heatmaps
  cos_heatmap <- create_lower_triangle_heatmap(cos_matrix, "Cosine Similarity")
  euc_heatmap <- create_lower_triangle_heatmap(euc_matrix, "Euclidean Distance")
  
  return(list(
    cosine_matrix = cos_matrix,
    euclidean_matrix = euc_matrix,
    all_cos_results = all_cos_results,
    all_euc_results = all_euc_results,
    plots = all_plots,
    heatmaps = list(cosine = cos_heatmap, euclidean = euc_heatmap)
  ))
}

# Function to analyze metrics and create plots
analyze_metric_pair <- function(cos_results, euc_results, title) {
  # Process cosine similarities
  cos_long <- cos_results %>%
    pivot_longer(cols = everything(),
                 names_to = "group",
                 values_to = "similarity")
  
  # Process euclidean distances
  euc_long <- euc_results %>%
    pivot_longer(cols = everything(),
                 names_to = "group",
                 values_to = "distance")
  
  # Create cosine plot
  p1 <- ggplot(cos_long, aes(x = group, y = similarity)) +
    geom_boxplot() +
    theme_classic() +
    labs(x = "Group", y = "Cosine Similarity", title = paste(title, "- Cosine Similarity")) +
    theme(plot.title = element_text(size = 10))
  
  # Create euclidean plot
  p2 <- ggplot(euc_long, aes(x = group, y = distance)) +
    geom_boxplot() +
    theme_classic() +
    labs(x = "Group", y = "Euclidean Distance", title = paste(title, "- Euclidean Distance")) +
    theme(plot.title = element_text(size = 10))
  
  return(ggarrange(p1, p2, ncol = 2))
}

# Function to create lower triangle heatmap
create_lower_triangle_heatmap <- function(matrix_data, title) {
  # Create lower triangle matrix
  matrix_data[upper.tri(matrix_data)] <- NA
  
  # Convert to long format for plotting
  melted_matrix <- melt(matrix_data, na.rm = TRUE)
  
  # Create heatmap with appropriate scale based on metric type
  p <- ggplot(melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    theme_minimal() +
    labs(title = title,
         x = "Base Group",
         y = "Comparison Group",
         fill = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add appropriate scale based on metric type
  if (grepl("Cosine", title)) {
    p <- p + scale_fill_viridis_c(limits = c(-1, 1))
  } else {
    p <- p + scale_fill_viridis_c()
  }
  
  return(p)
}

# Run the analysis
results <- bootstrap_pairwise_similarities(scaled_data, numerical_vars)

# Save all plots in root directory
for(group_name in names(results$plots)) {
  ggsave(
    filename = paste0(root, "bootstrap_comparison_", group_name, ".pdf"),
    plot = results$plots[[group_name]],
    width = 10,
    height = 6
  )
}

# Save heatmaps in root directory
ggsave(
  filename = paste0(root, "similarity_heatmaps.pdf"),
  plot = ggarrange(results$heatmaps$cosine, results$heatmaps$euclidean, 
                   ncol = 2, labels = c("A", "B")),
  width = 12,
  height = 6
)

# Save numerical results in root directory
write.csv(results$cosine_matrix, paste0(root, "cosine_similarity_matrix.csv"))
write.csv(results$euclidean_matrix, paste0(root, "euclidean_distance_matrix.csv"))

#############################################################################################
# Run PCA
pca_result <- prcomp(scaled_data[numerical_vars])

# Create PCA plot with ggplot2
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  group = scaled_data$group
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  theme_classic() +
  stat_ellipse(level = 0.95) +  # Using default level = 0.95
  # stat_ellipse(aes(linetype = group), level = 0.5) +  # 50% confidence
  # stat_ellipse(aes(linetype = group), level = 0.95, alpha = 0.5) +  # 95% confidence
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)"))


# Convex hull
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_polygon(data = pca_df %>% group_by(group) %>% slice(chull(PC1, PC2)), 
               aes(fill = group), alpha = 0.2) +
  theme_classic()

# Print variance explained and loadings
print(summary(pca_result))
print(pca_result$rotation[,1:3])

