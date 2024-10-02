import matplotlib.pyplot as plt
import numpy as np

# Sample data (replace these with actual data)
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
precipitation = [5.93, 3.86, 4.06, 3.00, 2.11, 1.57, 0.68, 0.82, 1.29, 3.70, 6.68, 5.53]  # in inches
temperature = [41, 43, 46, 50, 55, 60, 65, 65, 60, 52, 45, 41]  # in Fahrenheit

fig, ax1 = plt.subplots()

# Bar chart for precipitation
ax1.bar(months, precipitation, color='blue', alpha=0.7)
ax1.set_xlabel('Month')
ax1.set_ylabel('Precipitation (in)', color='blue')
ax1.set_ylim(0, max(precipitation) + 1)
ax1.tick_params(axis='y', labelcolor='blue')

# Line chart for temperature
ax2 = ax1.twinx()
ax2.plot(months, temperature, color='red', marker='o')
ax2.set_ylabel('Avg Temp (°F)', color='red')
ax2.set_ylim(0, max(temperature) + 10)
ax2.tick_params(axis='y', labelcolor='red')

# Titles and grid
ax1.set_title('Tacoma, WA')
ax1.grid(True, which='both', axis='y', linestyle='--', alpha=0.5)

plt.show()