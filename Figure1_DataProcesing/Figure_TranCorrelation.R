library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)

tran_path <- '/home/sam/scRNAseq/SCP509/cluster/RGC_Atlas_coordinates.txt'

tran <- read_delim(tran_path, skip = 1)%>%
  select(group...4)  %>%
  rename(scRNAseq = group...4) %>%
  mutate(scRNAseq = factor(paste0("T", as.numeric(sub("_.*", "", scRNAseq)))))

mouse <-read.csv("/home/sam/FinalRGC_xenium/volcan_expression_matrix.csv")   %>%
  select(Prediction)  %>%
  rename(Xenium = Prediction) %>%
  mutate(Xenium = factor(paste0("T", as.numeric(sub("_.*", "", Xenium)))))

# Calculate relative frequencies
scRNA_freq <- tran %>%
  mutate(scRNAseq = factor(scRNAseq, levels = paste0("T", 45:1))) %>%
  count(scRNAseq) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(scRNAseq))

xenium_freq <- mouse %>%
  mutate(Xenium = factor(Xenium, levels = paste0("T", 1:45))) %>%
  count(Xenium) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Xenium)

# Create a combined dataset
combined_freq <- data.frame(
  scRNA_freq = scRNA_freq$freq,
  xenium_freq = xenium_freq$freq,
  scRNA_type = scRNA_freq$scRNAseq,
  xenium_type = xenium_freq$Xenium
)

# Calculate R-squared
model <- lm(xenium_freq ~ scRNA_freq, data = combined_freq)
r_squared <- round(summary(model)$r.squared, 3)

# Create main scatter plot
main_plot <- ggplot(combined_freq, aes(x = scRNA_freq, y = xenium_freq)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text", x = max(combined_freq$scRNA_freq) * 0.8, 
           y = max(combined_freq$xenium_freq) * 0.9,
           label = paste("R² =", r_squared)) +
  scale_x_continuous(
    breaks = scRNA_freq$freq,
    labels = scRNA_freq$scRNAseq,
    limits = c(0, max(scRNA_freq$freq) * 1.1)
  ) +
  scale_y_continuous(
    breaks = xenium_freq$freq,
    labels = xenium_freq$Xenium,
    limits = c(0, max(xenium_freq$freq) * 1.1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1),
    # Stagger labels by adjusting margin
    axis.text = element_text(margin = margin(t = 10, r = 10)),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "scRNAseq Types", y = "Xenium Types",
       title = "Correlation of Cell Type Frequencies")

# Create right histogram
right_hist <- ggplot(xenium_freq, aes(x = freq, y = Xenium)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Relative Frequency", y = "Cell Type")

# Create bottom histogram with flipped order
bottom_hist <- ggplot(xenium_freq, aes(x = scRNA_freq$scRNAseq, y = -scRNA_freq$freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = function(x) abs(x)) +
  scale_x_discrete(limits = rev(levels(scRNA_freq$scRNAseq))) +  # Reverse order
  labs(x = "Cell Type", y = "Relative Frequency")

# Combine plots using grid.arrange
layout_matrix <- rbind(
  c(1, 1, 1, 2),
  c(1, 1, 1, 2),
  c(1, 1, 1, 2),
  c(3, 3, 3, NA)
)

combined_plot <- grid.arrange(
  main_plot, right_hist, bottom_hist,
  layout_matrix = layout_matrix,
  heights = c(3, 3, 3, 1)
)

# Export combined frequencies
write.csv(combined_freq, "/home/sam/FinalRGC_xenium/Tran_combined_frequencies.csv", row.names = FALSE)

# Create summary statistics dataframe
model_summary <- summary(model)
regression_stats <- data.frame(
  R_squared = model_summary$r.squared,
  Adjusted_R_squared = model_summary$adj.r.squared,
  F_statistic = model_summary$fstatistic[1],
  P_value = pf(model_summary$fstatistic[1], 
               model_summary$fstatistic[2], 
               model_summary$fstatistic[3], 
               lower.tail = FALSE),
  Intercept = model_summary$coefficients[1,1],
  Intercept_SE = model_summary$coefficients[1,2],
  Slope = model_summary$coefficients[2,1],
  Slope_SE = model_summary$coefficients[2,2]
)

# Export regression statistics
write.csv(regression_stats, "/home/sam/FinalRGC_xenium/Tran_regression_statistics.csv", row.names = FALSE)

# Print summary to console
print("Model Summary:")
print(model_summary)

print("\nRegression Statistics Saved:")
print(regression_stats)