# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)

# Load data
RNU5_cons <- read_csv("../../data/phyloP100_summary/RNU5_cons_phyloP100_.csv")

# View the data (optional, remove if running in a non-interactive environment)
# View(RNU5_cons)

# Remove rows with NA values in the Score column
RNU5_cons_clean <- RNU5_cons[!is.na(RNU5_cons$Score), ]

# Step 1: Separate gene into letter and number parts
RNU5_cons_clean$Gene_letter <- gsub("RNU5([A-Z]+)-.*", "\\1", RNU5_cons_clean$Gene)
RNU5_cons_clean$Gene_number <- gsub("RNU5[A-Z]+-(\\d+P?)", "\\1", RNU5_cons_clean$Gene)

# Step 2: Extract numeric part and "P" status for sorting
RNU5_cons_clean$Numeric_part <- as.numeric(gsub("P", "", RNU5_cons_clean$Gene_number)) # numeric portion only
RNU5_cons_clean$Has_P <- grepl("P$", RNU5_cons_clean$Gene_number) # logical column for "P" status

# Step 3: Sort by Gene_letter, then Numeric_part, then Has_P (to prioritize non-P over P)
RNU5_cons_clean <- RNU5_cons_clean[order(RNU5_cons_clean$Gene_letter,
                                         RNU5_cons_clean$Numeric_part,
                                         RNU5_cons_clean$Has_P), ]

# Step 4: Optional - Drop intermediate columns if needed
RNU5_cons_clean$Numeric_part <- NULL
RNU5_cons_clean$Has_P <- NULL

# Calculate summary metrics: median and max conservation for each gene
summary_metrics_RNU5 <- RNU5_cons_clean %>%
  group_by(Gene) %>%
  summarise(
    Median_Conservation = median(Score, na.rm = TRUE),
    Max_Conservation = max(Score, na.rm = TRUE)
  )

# Optionally, print the summary metrics to the console
print(summary_metrics_RNU5)

# Create the plot (no summary metrics on plot)
plot_RNU5 <- ggplot(RNU5_cons_clean, aes(x = Score, y = Gene)) +
  geom_boxplot(fill = "olivedrab2", color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  labs(title = "Distribution of Scores by Gene (PhyloP100way)",
       x = "Score",
       y = "Gene") +
  theme_minimal()

# Save the plot to a file
ggsave("../../results/RNU5__phylo100P_score_plot.pdf", plot = plot_RNU5, width = 12, height = 8)

# Optionally: Save summary metrics to a CSV file for RNU5
write_csv(summary_metrics_RNU5, "../../data/phyloP100_summary/RNU5_phyloP100_summary_metrics.csv")
