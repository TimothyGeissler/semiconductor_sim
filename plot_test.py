import matplotlib.pyplot as plt
import numpy as np

# Generate sample data for the plots
x = np.linspace(0, 5)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = x
y4 = x**2
y5 = np.sqrt(x)

# Create a figure with a single set of axes
fig, ax = plt.subplots(figsize=(10, 6))

# Overlay each plot on the same axes
ax.plot(x, y1, label='Plot 1')
ax.plot(x, y2, label='Plot 2')
ax.plot(x, y3, label='Plot 3')
ax.plot(x, y4, label='Plot 4')
ax.plot(x, y5, label='Plot 5')

# Add legend
ax.legend()

# Set labels and title
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_title('Overlay of Plots on the Same Axes')

# Show the plot
plt.show()
