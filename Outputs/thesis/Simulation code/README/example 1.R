# Load necessary package
library(GPArotation)

# Step 1: Create a non-sparse loading matrix
# Let's define a 6x2 matrix with arbitrary values (non-sparse)
non_sparse_matrix <- matrix(c(0.4, 0.3, 0.6, 0.7, 0.1, 0.4,0.8, 0.5, 0.3, 0.1, 0.8, 0.5), nrow = 6, ncol = 2, byrow = TRUE)

# Display the original non-sparse matrix
print("Original Non-Sparse Loading Matrix:")
print(non_sparse_matrix)

# Step 2: Apply Varimax Rotation
rotated_result <- Varimax(non_sparse_matrix)

# Step 3: Extract and display the rotated (sparser) loading matrix
rotated_matrix <- rotated_result$loadings

rotated_matrix[abs(rotated_matrix) < 0.05] = 0
print(rotated_matrix)
