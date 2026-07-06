import pandas as pd
import pwlf
import numpy as np
import matplotlib.pyplot as plt

# Load your allele frequency data (time and frequency columns)
data = pd.read_csv('allele_frequency_data.csv')
time = data['Time'].values  # X-axis (e.g., years or generations)
frequency = data['Allele_Frequency'].values  # Y-axis (allele frequencies)

# Initialize the piecewise linear model
model = pwlf.PiecewiseLinFit(time, frequency)

# Let the model find the optimal breakpoint (1 breakpoint implies 2 segments)
n_breakpoints = 1  # You can adjust this to increase the number of segments
breakpoints = model.fit(n_breakpoints)

# Generate predictions for visualization
time_hat = np.linspace(min(time), max(time), num=100)
frequency_hat = model.predict(time_hat)

# Plot the results
plt.scatter(time, frequency, color='blue', label='Observed Data')
plt.plot(time_hat, frequency_hat, color='red', label='Piecewise Fit')
plt.title('Broken Stick Model for Allele Frequency with Automatic Breakpoint')
plt.xlabel('Time (years)')
plt.ylabel('Allele Frequency')

# Mark the detected breakpoint
for bp in breakpoints[1:-1]:  # Skip the start and end points (they're not actual breakpoints)
    plt.axvline(x=bp, color='green', linestyle='--', label=f'Breakpoint at {bp:.2f}')

plt.legend()
plt.show()

# Print the found breakpoints
print("Breakpoints (knots) detected:", breakpoints)
