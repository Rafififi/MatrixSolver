import matplotlib.pyplot as plt
import numpy as np
from scipy.io import mmread

# Read the .mtx file
matrix = mmread('file.mtx')

# Get the row and column indices and the values of the non-zero elements
row, col = matrix.nonzero()
data = matrix.data

# Create a scatter plot of the non-zero elements
plt.figure(figsize=(10, 10))
plt.scatter(row, col, c=data, s=1, cmap='hot')
plt.gca().invert_yaxis()
plt.show()