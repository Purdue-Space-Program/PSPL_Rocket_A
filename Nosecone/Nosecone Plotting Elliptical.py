import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

R = 3  # Outer radius
ratio = 5  # Ratio of length to radius
L = R * ratio  # Length of nosecone
base_h = 77.42854035000 # Base height of rocket up till nosecone

x = np.linspace(0, L, 20)  # Initialize x values
y = R * np.sqrt(1 - (x**2) / (L**2))  # Calculate y values based on ellipse equation
z = np.zeros_like(x)  # z values are zero for 2D plot

z_temp = z # Store original z values
z = x + base_h # Shift z values to start from base height
x = z_temp # Set x values to original z values

# Create a DataFrame
df = pd.DataFrame({'x': x, 'y': y, 'z': z})

# Save to CSV
df.to_csv('nosecone_elliptical.csv', index=False, header=False)
