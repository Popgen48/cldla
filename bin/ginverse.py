# python3 ginverse.py <input_file> <output_file>

from scipy import linalg
import numpy
import sys

matrix = numpy.empty((0, 0))

with open(sys.argv[1], 'r') as input_file:
    for line in input_file:
        row = float(line.split()[0]) - 1
        col = float(line.split()[1]) - 1
        value = float(line.split()[2])
        max_row = int(row) + 1
        max_col = int(col) + 1
        
        # dynamically resize the matrix
        if max_row > matrix.shape[0]:
            matrix = numpy.pad(matrix, ((0, max_row - matrix.shape[0]), (0, 0)), mode='constant')
        if max_col > matrix.shape[1]:
            matrix = numpy.pad(matrix, ((0, 0), (0, max_col - matrix.shape[1])), mode='constant')
        
        matrix[int(row), int(col)] = value   

# matrix is symmetrical
for i in range(matrix.shape[0]):
    for j in range(i+1, matrix.shape[1]):
        matrix[i][j] = matrix[j][i]

# calculate the g-inverse (https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pinv.html)      
inverse = linalg.pinv(matrix)

with open(sys.argv[2], 'w') as output_file:
    for i in range(inverse.shape[0]):
        for j in range(i+1):
            output_file.write(f'{i+1}\t{j+1}\t{inverse[i][j]}\n')
            

