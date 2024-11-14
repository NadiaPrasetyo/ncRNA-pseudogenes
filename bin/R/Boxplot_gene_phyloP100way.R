# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)

# Define the dataset names and corresponding colors
#datasets <- c("RNU1", "RNU2", "RNU4", "RNU6")
#colors <- c("paleturquoise2", "rosybrown1", "lemonchiffon1", "thistle1")  # Customize the colors as needed

datasets <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA")
colors <- c("palevioletred1", "darkturquoise", "tan3", "slategray", "springgreen4")

# Loop through each dataset
for (i in seq_along(datasets)) {
  # Get dataset and color for the current iteration
  dataset <- datasets[i]
  plot_color <- colors[i]
  
  # Construct the file path for the dataset and the plot
  data_path <- paste0("../../data/phyloP100_summary/", dataset, "_cons_phyloP100_.csv")
  plot_path <- paste0("../../results/", dataset, "_phyloP100_score_plot.pdf")
  
  # Load data
  dataset_cons <- read_csv(data_path)
  
  # Remove rows with NA values in the Score column
  dataset_cons_clean <- dataset_cons[!is.na(dataset_cons$Score), ]
  
  # 1. Extract the numeric part from the Gene column
  # Combined extraction for conventions: RNU#-#(P), dataset#(P), and VTRNA#-#P
  dataset_cons_clean$Gene_number <- as.numeric(
    ifelse(
      grepl("VTRNA-(\\d+)P", dataset_cons_clean$Gene),
      gsub("VTRNA-(\\d+)P", "\\1", dataset_cons_clean$Gene),
      ifelse(
        grepl(paste0(dataset, "-(\\d+)P?"), dataset_cons_clean$Gene),
        gsub(paste0(dataset, "-(\\d+)P?"), "\\1", dataset_cons_clean$Gene),
        gsub(paste0(dataset, "(\\d+)P?"), "\\1", dataset_cons_clean$Gene)
      )
    )
  )
  
  # 2. Reorder the Gene factor based on the extracted number
  dataset_cons_clean$Gene <- factor(dataset_cons_clean$Gene, levels = unique(dataset_cons_clean$Gene[order(dataset_cons_clean$Gene_number)]))
  
  # Calculate summary metrics: median and max conservation for each gene
  summary_metrics <- dataset_cons_clean %>%
    group_by(Gene) %>%
    summarise(
      Median_Conservation = median(Score, na.rm = TRUE),
      Max_Conservation = max(Score, na.rm = TRUE)
    )
  
  # Optionally, print the summary metrics to the console
  print(summary_metrics)
  
  # Create the plot with a unique color for each dataset
  plot <- ggplot(dataset_cons_clean, aes(x = Score, y = Gene)) +
    geom_boxplot(fill = plot_color, color = "black") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    labs(title = paste("Distribution of Scores by Gene (phyloP100way) for", dataset),
         x = "Score",
         y = "Gene") +
    theme_minimal()
  
  # Save the plot to a file
  ggsave(plot_path, plot = plot, width = 12, height = 10)
  
  # Optionally: Save summary metrics to a CSV file for each dataset
  summary_file <- paste0("../../data/phyloP100_summary/", dataset, "_phyloP100_summary_metrics.csv")
  write_csv(summary_metrics, summary_file)
}
