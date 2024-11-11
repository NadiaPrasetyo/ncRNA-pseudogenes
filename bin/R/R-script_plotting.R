# Load necessary libraries
library(readr)
library(ggplot2)

# Load data
RNU2_cons <- read_csv("../../data/RNU2_cons.csv")

# View the data (optional, remove if running in a non-interactive environment)
# View(RNU2_cons)

# Remove rows with NA values in the Score column
RNU2_cons_clean <- RNU2_cons[!is.na(RNU2_cons$Score), ]

# 1. Extract the numeric part from the Gene column
RNU2_cons_clean$Gene_number <- as.numeric(gsub("RNU2-(\\d+)P?", "\\1", RNU2_cons_clean$Gene))

# 2. Reorder the Gene factor based on the extracted number
RNU2_cons_clean$Gene <- factor(RNU2_cons_clean$Gene, levels = unique(RNU2_cons_clean$Gene[order(RNU2_cons_clean$Gene_number)]))

# Create the plot
plot_RNU2 <- ggplot(RNU2_cons_clean, aes(x = Score, y = Gene)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  labs(title = "Distribution of Scores by Gene (Ordered)",
       x = "Score",
       y = "Gene") +
  theme_minimal()

# Save the plot to a file
ggsave("../../results/RNU2_score_plot_vertical_labels_large.pdf", plot = plot_RNU2, width = 12, height = 8)
