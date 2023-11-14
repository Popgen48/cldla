library("mbend")

# Check the number of command-line arguments
if (length(commandArgs(trailingOnly = TRUE)) < 2) {
  cat("Usage: Rscript script.R <input_file> <output_file> <method>([lrs]/hj)\nTo know more about which method to use, refer to the description of \'https://cran.r-project.org/web/packages/mbend/mbend.pdf\'\n")
  quit("no")
}

# Retrieve command-line arguments
input_file <- commandArgs(trailingOnly = TRUE)[[1]]
output_file <- commandArgs(trailingOnly = TRUE)[[2]]

# Default method is 'lrs'
m <- "lrs"

# If a method argument is provided, use it instead
if (length(commandArgs(trailingOnly = TRUE)) >= 3) {
  m <- commandArgs(trailingOnly = TRUE)[[3]]
}

# read matrix from file
matrix_data <- read.table(input_file, header = FALSE, skip=1)

# get the last record for size
last_entry <- tail(matrix_data, n = 1)
n = last_entry[[1]]

mat <- matrix(0, nrow = n, ncol = n)

for (i in 1:nrow(matrix_data)) {
  row_idx <- matrix_data[[1]][i]
  col_idx <- matrix_data[[2]][i]
  value <- matrix_data[[3]][i]
  mat[row_idx, col_idx] <- value
  mat[col_idx, row_idx] <- value  # Set the value at [col, row] to account for symmetry
}

# Perform bending using the specified method
if (m != "lrs" && m != "hj") {
  cat("Invalid bending method. Supported methods: lrs(Schaeffer) and hj(Jorjani)\nMore about them at \'https://cran.r-project.org/web/packages/mbend/mbend.pdf\'\n")
  quit("no")
}

bended = bend(mat, method=m)

# Open the file for writing
file_conn <- file(output_file, "w")

# Write the number of elements (assuming it's the same as the size of the matrix)
# cat(n, "\n", file = file_conn)

# Iterate through the matrix and write its content in the "row col value" format
for (row_idx in 1:n) {
    for (col_idx in 1:row_idx) {
        value <- mat[row_idx, col_idx]
        cat(row_idx, col_idx, value, "\n", file = file_conn)
  }
}

# Close the file
close(file_conn)