# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)

# Define the dataset names and corresponding colors
#datasets <- c("RNU1", "RNU2", "RNU4", "RNU6")
#colors <- c("lightblue", "lightpink", "lightyellow", "plum2")  # Customize the colors as needed

#datasets <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA")
#colors <- c("hotpink1", "cyan2", "darkorange", "azure3", "springgreen3")

datasets <- c("RNY","TRNA","RN7SL","RNU7","RN7SK")
colors <- c("dodgerblue", "darkmagenta", "lightsalmon", "red", "darkgreen")

# Loop through each dataset
for (i in seq_along(datasets)) {
  # Get dataset and color for the current iteration
  dataset <- datasets[i]
  plot_color <- colors[i]
  
  # Construct the file path for the dataset and the plot
  data_path <- paste0("../../data/phastCons30_summary/", dataset, "_cons.csv")
  plot_path <- paste0("../../results/", dataset, "_phastCons30_score_plot.pdf")
  
  # Load data
  dataset_cons <- read_csv(data_path)
  
  # Remove rows with NA values in the Score column
  dataset_cons_clean <- dataset_cons[!is.na(dataset_cons$Score), ]
  
  # 1. Extract the numeric part from the Gene column
  # Combined extraction for conventions: VTRNA#-#, RNU#-#(P), dataset#(P), RNY#(P)#, RN7SL#(P), RN7SK(P)#
  dataset_cons_clean$Gene_number <- as.numeric(
    ifelse(
      grepl("VTRNA\\d+-(\\d+)P?", dataset_cons_clean$Gene),
      gsub("VTRNA\\d+-(\\d+)P?", "\\1", dataset_cons_clean$Gene),
      ifelse(
        grepl(paste0(dataset, "-(\\d+)P?"), dataset_cons_clean$Gene),
        gsub(paste0(dataset, "-(\\d+)P?"), "\\1", dataset_cons_clean$Gene),
        ifelse(
          grepl("RNY(\\d+)(P?)(\\d*)", dataset_cons_clean$Gene),
          # Extract second number if present, else default to 1
          ifelse(
            nchar(gsub("RNY(\\d+)(P?)(\\d*)", "\\3", dataset_cons_clean$Gene)) > 0,
            gsub("RNY(\\d+)(P?)(\\d*)", "\\3", dataset_cons_clean$Gene),
            "1"
          ),
          ifelse(
            grepl("RN7SL(\\d+)(P?)", dataset_cons_clean$Gene),
            gsub("RN7SL(\\d+)(P?)", "\\1", dataset_cons_clean$Gene),
            ifelse(
              grepl("RN7SK(P?)(\\d+)", dataset_cons_clean$Gene),
              gsub("RN7SK(P?)(\\d+)", "\\2", dataset_cons_clean$Gene),
              # Default case for other formats
              gsub(paste0(dataset, "(\\d+)P?"), "\\1", dataset_cons_clean$Gene)
            )
          )
        )
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
    labs(title = paste("Distribution of Scores by Gene (phastCons30way) for", dataset),
         x = "Score",
         y = "Gene") +
    theme_minimal()
  
  # Save the plot to a file
  ggsave(plot_path, plot = plot, width = 12, height = 50, limitsize = FALSE)
  
  # Save summary metrics to a CSV file for each dataset
  summary_file <- paste0("../../data/phastCons30_summary/", dataset, "_phastCons30_summary_metrics.csv")
  write_csv(summary_metrics, summary_file)
}
