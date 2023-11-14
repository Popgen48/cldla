import numpy as np
import sys

# Function to read a matrix from a text file
def readMatrixFromFile(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        num_rows = int(lines[0])
        matrix = np.zeros((num_rows, num_rows), dtype=float)
        for line in lines[1:]:
            row, col, value = map(float, line.split())
            matrix[int(row) - 1, int(col) - 1] = value
            matrix[int(col) - 1, int(row) - 1] = value
        print(f"Matrix size: {num_rows} × {num_rows}")
        return matrix

# Function to write a matrix to a text file
def writeMatrixToFile(matrix, filename):
    num_rows = matrix.shape[0]
    with open(filename, 'w') as file:
        # file.write(f"{num_rows}\n")
        for row in range(num_rows):
            for col in range(row + 1):
                file.write(f"{row + 1} {col + 1} {matrix[row, col]}\n")

# Read the input and output file names from command-line arguments
if len(sys.argv) != 3:
    print("Usage: python3 bend_and_inverse.py <input_file> <output_file>")
    sys.exit(1)

input_filename = sys.argv[1]
output_filename = sys.argv[2]

# Read the input matrix
A = readMatrixFromFile(input_filename)

# Perform eigenvalue decomposition
eigenvalues, eigenvectors = np.linalg.eigh(A)

# Check for negative eigenvalues
has_negative_eigenvalues = any(eigenvalues < 0)

if has_negative_eigenvalues:
    print("The matrix has negative eigenvalues. Proceeding to bending...")
    
    # Step 1: Calculate 's = (sum of all neg eigenvals) × 2'
    s = 2.0 * np.sum(eigenvalues[eigenvalues < 0])

    # Step 2: Calculate 't = (s × s) × 100 + 1'
    t = (s**2) * 100.0 + 1.0

    # Step 3: Replace negative eigenvalues (n's) using 'n_new = p × (s − n) × (s − n)/t' where p is the smallest positive eigenvalue
    eigenvalues[eigenvalues < 0] = np.min(eigenvalues[eigenvalues > 0]) * (s - eigenvalues[eigenvalues < 0]) * (s - eigenvalues[eigenvalues < 0]) / t

    # Print the modified eigenvalues
    # print("Modified Eigenvalues:")
    # print(eigenvalues)

    # Construct the new matrix B
    B = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T

    # Save the new matrix B
    # writeMatrixToFile(B, output_filename)
    # print(f"Matrix bending finished. Output saved to '{output_filename}'.")
    
    # Calculate the generalized inverse 
    B_gen_inverse = np.linalg.pinv(B)

    # Print the generalized inverse
    # print("Generalized Inverse (Pseudo-Inverse) of B:")
    # print(B_gen_inverse)

    # Save the generalized inverse to a .giv file
    writeMatrixToFile(B_gen_inverse, output_filename)
    print(f"Generalized inverse saved to {output_filename}.")
    
else:
    print("No negative eigenvalues found. Matrix can be used as it is for inversion.")
    
    # No negative eigenvalues found; remove the first line and save
    # with open(input_filename, 'r') as infile:
    #     first_line = infile.readline()
    #     with open(output_filename, 'w') as outfile:
    #         outfile.write(infile.read())
    # print(f"Output saved to '{output_filename}'.")
    
    # Calculate the generalized inverse 
    A_gen_inverse = np.linalg.pinv(A)

    # Print the generalized inverse
    # print("Generalized Inverse (Pseudo-Inverse) of B:")
    # print(B_gen_inverse)

    # Save the generalized inverse to a .giv file
    writeMatrixToFile(A_gen_inverse, output_filename)
    print(f"Generalized inverse saved to {output_filename}.")
    
    
