import matplotlib.pyplot as plt
import numpy as np

# Create some sample data (negative values)
x = np.linspace(-10, -1, 100)
y = -np.log10(-x)  # Working with the absolute values for logarithmic scale

# Create a plot
fig, ax = plt.subplots()
ax.plot(x, y)

# Set the y-axis scale to logarithmic
ax.set_yscale('log')

# Adjust the y-axis tick labels to reflect the original signs
yticks = [-10**i for i in range(1, 11)]
yticklabels = [f'-10^{i}' for i in range(1, 11)]
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)

# Show the plot
plt.show()