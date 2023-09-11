import numpy as np

# Create a sample numpy array
data = np.arange(25).reshape(5, 5)
print(data)
# Exclude the first two diagonal elements
np.fill_diagonal(data[2:, 2:], data[2:, 2:] * 0.5)

print(data)
