import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load the NetCDF file
file_path = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Temp2022SSW1.nc")
data = xr.open_dataset(file_path, decode_times=True)

# Select data for 10 hPa pressure level and latitude range 90N to 60N
pressure_level = 10.0
latitude_range = slice(90, 60)  # Selecting latitudes from 90N to 60N
temperature_data = data.sel(pressure_level=pressure_level, latitude=latitude_range).t

# Convert temperature from Kelvin to Celsius if needed (for this example it's left in Kelvin)
# temperature_data = temperature_data - 273.15

# Calculate daily mean temperature over the latitude and longitude range
daily_mean_temp = temperature_data.mean(dim=['latitude', 'longitude'])

# Calculate statistics for shading (min, max, percentiles)
mean_temp = daily_mean_temp.mean(dim='valid_time')
max_temp = daily_mean_temp.max(dim='valid_time')
min_temp = daily_mean_temp.min(dim='valid_time')
percentiles = daily_mean_temp.quantile([0.1, 0.3, 0.7, 0.9], dim='valid_time')

# Create a range of days for plotting
days_of_year = np.arange(len(daily_mean_temp.valid_time))

# Plot configuration
plt.figure(figsize=(14, 7))
plt.plot(days_of_year, daily_mean_temp, color='green', label='Daily Mean Temperature')
plt.plot(days_of_year, [mean_temp] * len(days_of_year), color='yellow', label='Mean Temperature', linewidth=2)
plt.fill_between(days_of_year, percentiles.sel(quantile=0.1), percentiles.sel(quantile=0.9), color='gray', alpha=0.3, label='10%-90% Range')
plt.fill_between(days_of_year, percentiles.sel(quantile=0.3), percentiles.sel(quantile=0.7), color='gray', alpha=0.6, label='30%-70% Range')
plt.fill_between(days_of_year, min_temp, max_temp, color='gray', alpha=0.1, label='Min-Max Range')

# Plot details
plt.title(f'Temperature [K] at {pressure_level} hPa (90N to 60N)')
plt.xlabel('Day of Year')
plt.ylabel('Temperature (K)')
plt.legend()
plt.grid(True)
plt.show()
