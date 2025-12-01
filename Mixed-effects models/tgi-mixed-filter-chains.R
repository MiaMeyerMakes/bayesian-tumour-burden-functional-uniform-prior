library(ggplot2)
library(tidyr)
library(dplyr)

###############################################################################
numfiltered_tgi <- 0
files_to_process <- c(
  "Nov-03-13h41-DataV2-mixed-effects-growth-posterior-sample-uniform-NARROW-15.csv",
  "Nov-03-13h41-DataV2-mixed-effects-shrink-posterior-sample-uniform-NARROW-15.csv",
  "Nov-07-00h52-DataV2-mixed-effects-growth-posterior-sample-uniform-NARROW-50.csv",
  "Nov-07-00h52-DataV2-mixed-effects-shrink-posterior-sample-uniform-NARROW-50.csv",
  "Nov-07-11h41-DataV2-mixed-effects-growth-posterior-sample-uniform-NARROW-100.csv",
  "Nov-07-11h41-DataV2-mixed-effects-shrink-posterior-sample-uniform-NARROW-100.csv",

  "Nov-09-17h17-DataV2-mixed-effects-growth-posterior-sample-uniform-MEDIUM-15.csv",
  "Nov-09-17h17-DataV2-mixed-effects-shrink-posterior-sample-uniform-MEDIUM-15.csv",
  "Nov-10-01h29-DataV2-mixed-effects-growth-posterior-sample-uniform-MEDIUM-50.csv",
  "Nov-10-01h29-DataV2-mixed-effects-shrink-posterior-sample-uniform-MEDIUM-50.csv",
  "Nov-10-15h31-DataV2-mixed-effects-growth-posterior-sample-uniform-MEDIUM-100.csv",
  "Nov-10-15h31-DataV2-mixed-effects-shrink-posterior-sample-uniform-MEDIUM-100.csv",

  "Sep-15-22h37-DataV2-mixed-effects-growth-posterior-sample-uniform-WIDE-15.csv",
  "Sep-15-22h37-DataV2-mixed-effects-shrink-posterior-sample-uniform-WIDE-15.csv",
  "Nov-07-20h28-DataV2-mixed-effects-growth-posterior-sample-uniform-WIDE-50.csv",
  "Nov-07-20h28-DataV2-mixed-effects-shrink-posterior-sample-uniform-WIDE-50.csv",
  "Nov-08-17h17-DataV2-mixed-effects-growth-posterior-sample-uniform-WIDE-100.csv",
  "Nov-08-17h17-DataV2-mixed-effects-shrink-posterior-sample-uniform-WIDE-100.csv",
  
  "Nov-18-01h28-DataV2-mixed-effects-growth-posterior-sample-jeffreys15.csv",
  "Nov-18-01h28-DataV2-mixed-effects-shrink-posterior-sample-jeffreys15.csv",
  "Nov-18-18h01-DataV2-mixed-effects-growth-posterior-sample-jeffreys50.csv",
  "Nov-18-18h01-DataV2-mixed-effects-shrink-posterior-sample-jeffreys50.csv",
  "Nov-19-14h17-DataV2-mixed-effects-growth-posterior-sample-jeffreys100.csv",
  "Nov-19-14h17-DataV2-mixed-effects-shrink-posterior-sample-jeffreys100.csv",
  
  "Sep-18-01h15-DataV2-mixed-effects-growth-posterior-sample-FUP-15.csv",
  "Sep-18-01h15-DataV2-mixed-effects-shrink-posterior-sample-FUP-15.csv",
  "Sep-18-12h14-DataV2-mixed-effects-growth-posterior-sample-FUP-50.csv",
  "Sep-18-12h14-DataV2-mixed-effects-shrink-posterior-sample-FUP-50.csv",
  "Sep-19-06h09-DataV2-mixed-effects-growth-posterior-sample-FUP-100.csv",
  "Sep-19-06h09-DataV2-mixed-effects-shrink-posterior-sample-FUP-100.csv"
)

# 2. Set the standard deviation threshold.
#    Columns with a standard deviation below this value will be removed.
variance_threshold <- 2e-05

# --- Processing Loop ---

# Loop through each file name in the list.
for (file_path in files_to_process) {
  
  # Check if the file actually exists to avoid errors.
  if (!file.exists(file_path)) {
    warning(paste("File not found, skipping:", file_path))
    next # Skip to the next file in the loop
  }
  
  cat(paste("Processing:", file_path, "\n"))
  
  # 1. Read the dataset from the CSV file.
  #    check.names = FALSE prevents R from changing column names that might be
  #    syntactically invalid (e.g., starting with a number).
  original_data <- read.csv(file_path, header = TRUE, check.names = FALSE)
  cat("Stationary chains: ",sum(apply(original_data, 2, sd)<variance_threshold),"\n")
  
  # 2. Calculate the standard deviation for every column.
  #    The 'apply' function runs the 'sd' function on each column (MARGIN = 2).
  sds <- apply(original_data, 2, sd, na.rm = TRUE)
  
  # 3. Identify which columns to KEEP.
  #    This creates a logical vector (TRUE/FALSE) where TRUE corresponds to
  #    columns with a standard deviation greater than or equal to the threshold.
  columns_to_keep <- sds >= variance_threshold
  cat(paste("Minimum sd: ",min(apply(original_data,2,sd)),"\n"))
  cat(paste("Maximum sd: ",max(apply(original_data,2,sd)),"\n"))
  cat(paste("Mean sd: ",mean(apply(original_data,2,sd)),"\n"))
  
  # Report the number of columns being removed.
  num_removed <- sum(!columns_to_keep)
  numfiltered_tgi <- numfiltered_tgi + num_removed
  cat(paste("  - Found and removed", num_removed, "columns with low variance.\n"))
  
  # 4. Filter the original dataset to only include the desired columns.
  filtered_data <- original_data[, columns_to_keep]
  
  # 5. Create the new file name.
  #    This uses 'sub' to replace the ".csv" extension with "_filtered.csv".
  new_file_path <- sub("\\.csv$", "_filtered.csv", file_path)
  
  # 6. Save the new, filtered dataset to a CSV file.
  #    'row.names = FALSE' prevents R from writing an extra column for row numbers.
  write.csv(filtered_data, new_file_path, row.names = FALSE)
  
  cat(paste("  - Saved filtered data to:", new_file_path, "\n\n"))
}

cat("--- Filtering complete for all files. ---\n")
numfiltered_tgi
