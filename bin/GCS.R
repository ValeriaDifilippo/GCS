## To run the script Rscript script.R /path/to/working/directory
library(ggplot2)
library(dplyr)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure a working directory is provided as an argument
if (length(args) < 1) {
  stop("Please provide the working directory as an argument.")
}

# First argument: working directory
working_dir <- args[1]

# Set the working directory
setwd(working_dir)
cat("Working directory set to:", getwd(), "\n")

# List .txt files from the "data" subfolder in the working directory
filenames <- list.files(path = file.path(working_dir, "input_files"), recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)

for (f in filenames) {
  print(f)
  
  # Extract sample_ID by removing the .txt extension from the file name
  sample_ID <- sub(pattern = "\\.txt$", replacement = "", basename(f))
  
  # Create the directory path
  dir_path <- file.path(working_dir, "results", sample_ID)
  
  # Create the directory if it doesn't already exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("Directory created:", dir_path, "\n")
  } else {
    cat("Directory already exists:", dir_path, "\n")
  }
  
  # Create the plots subfolder inside the sample-specific directory
  plots_dir <- file.path(dir_path, "plots")
  
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
    cat("Plots directory created:", plots_dir, "\n")
  } else {
    cat("Plots directory already exists:", plots_dir, "\n")
  }
  
  # Read the file (Assuming tab-separated values and the file has a header)
  data <- read.table(f, header = TRUE, sep = "\t", as.is = TRUE)
  data_clean <- na.omit(data)
  
  # Ensure the column `LogR` exists in the data
  if (!"LogR" %in% colnames(data_clean)) {
    warning("Column 'LogR' not found in file:", f)
    next
  }
  
  # Remove the outliers by calculating the quantile (Q1, Q3)
  Q1 <- quantile(data_clean$LogR, 0.25, na.rm = TRUE)
  Q3 <- quantile(data_clean$LogR, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data_no_outliers_iqr <- data_clean[data_clean >= lower_bound & data_clean <= upper_bound]
  
  # Calculate mean and standard deviation without outliers
  mean_no_outliers_iqr <- mean(data_no_outliers_iqr, na.rm = TRUE)
  sd_no_outliers_iqr <- sd(data_no_outliers_iqr, na.rm = TRUE)
  
  # Calculate the standard deviation normalized for all the genome
  sd1 <- mean_no_outliers_iqr + (sd_no_outliers_iqr * 1.96)
  sd2 <- mean_no_outliers_iqr - (sd_no_outliers_iqr * 1.96)
  sdnorm <- (sd1 +- sd2) / 2
  
  print(sdnorm)
  
  # Create a data frame with results for the whole genome
  data_frame <- data.frame(Case = sample_ID, GCS_WholeGenome = sdnorm)
  
  # Write the results to a file in the results directory
  result_file <- file.path(dir_path, paste0(sample_ID, "_whole_genome_results.txt"))
  write.table(data_frame, file = result_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  cat("Results written to:", result_file, "\n")
  
  # Split the data by chromosome
  sptdf <- split(data_clean, data_clean$Chromosome)
  
  # Initialize an empty data frame to collect results for all chromosomes
  all_chr_results <- data.frame()
  
  for (chr in names(sptdf)) {
    chr_data <- sptdf[[chr]]
    cat("Calculating the GCS for chromosome",chr,". ")
    
    # Ensure chromosome names are handled properly
    chrom <- gsub("[/\\*?<>|:]", "_", chr)
    
    # Remove the outliers by calculating the quantile (Q1, Q3)
    Q1_chr <- quantile(chr_data$LogR, 0.25, na.rm = TRUE)
    Q3_chr <- quantile(chr_data$LogR, 0.75, na.rm = TRUE)
    IQR_chr <- Q3_chr - Q1_chr
    lower_bound_chr <- Q1_chr - 1.5 * IQR_chr
    upper_bound_chr <- Q3_chr + 1.5 * IQR_chr
    data_no_outliers_iqr_chr <- chr_data[chr_data >= lower_bound_chr & chr_data <= upper_bound_chr]
    
    # Calculate mean and standard deviation without outliers
    mean_no_outliers_iqr_chr <- mean(data_no_outliers_iqr_chr, na.rm = TRUE)
    sd_no_outliers_iqr_chr <- sd(data_no_outliers_iqr_chr, na.rm = TRUE)
    
    # Calculate the standard deviation normalized for all the genome
    sd1_chr <- mean_no_outliers_iqr_chr + (sd_no_outliers_iqr_chr * 1.96)
    sd2_chr <- mean_no_outliers_iqr_chr - (sd_no_outliers_iqr_chr * 1.96)
    sdnorm_chr <- (sd1_chr +- sd2_chr) / 2
   
    data_frame_chr <- data.frame(Case = sample_ID, GCS = sdnorm_chr, Chr = chrom, SD_pos = sd1_chr, SD_neg = sd2_chr)
    
    # Append results to the cumulative data frame
    all_chr_results <- rbind(all_chr_results, data_frame_chr)
    
    # Plot for each chromosome
    p <- ggplot(chr_data, aes(x = Position, y = LogR)) +
      geom_point(color = "skyblue3", size = 0.1) +  # Set dot color and size
      geom_hline(yintercept = 0, linetype = "dashed", color = "snow4", show.legend = TRUE) +
      #geom_hline(yintercept = sdnorm, linetype = "dashed", color = "red", show.legend = TRUE) +
      geom_hline(yintercept = sd1_chr, linetype = "dashed", color = "royalblue4", show.legend = TRUE) +
      geom_hline(yintercept = sd2_chr, linetype = "dashed", color = "royalblue4", show.legend = TRUE) +
      labs(title = paste("Chromosome", chrom, " GCS= ",sdnorm_chr), x = "Position", y = "LogR") +
      theme_minimal() +
      theme(panel.grid = element_blank(),  # Remove grid lines
            panel.background = element_rect(fill = "white"),  # White background
            legend.position = "bottom")  # Position legend at the bottom
    
    # Save plot to file
    plot_file <- file.path(plots_dir, paste0(sample_ID, "_chr_", chrom, "_plot.png"))
    ggsave(filename = plot_file, plot = p, width = 10, height = 6)
    
    cat("Plot saved to:", plot_file, "\n")
  }
  
  # Write all chromosome results to a single file
  result_file_all_chr <- file.path(dir_path, paste0(sample_ID, "_all_chromosomes_results.txt"))
  write.table(all_chr_results, file = result_file_all_chr, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  cat("Chromosome-specific results written to:", result_file_all_chr, "\n")
  
  # Calculate the mean of the SD values across all chromosomes. The mean represents the GCS.
  GCS <- mean(all_chr_results$GCS, na.rm = TRUE)
  
  # Create a data frame for GCS
  data_frame_GCS <- data.frame(Case = sample_ID, GCS = GCS)
  result_GCS <- file.path(dir_path, paste0(sample_ID, "_GCS.txt"))
  write.table(data_frame_GCS, file = result_GCS, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  cat("Global Chromosome SD written to:", result_GCS, "\n")
}
