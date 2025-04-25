import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

# Step 1: Load the dataset
file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Snow_2223.nc"
dataset = xr.open_dataset(file_path)

# Step 2: Verify the 'valid_time' coordinate and use it directly if it’s already in datetime format
times_converted = dataset['valid_time'].values  # Get the time values directly

# Step 3: Select the data between 2022-07-01 and 2023-07-01
start_date = "2022-07-01"
end_date = "2023-07-01"
start_datetime = np.datetime64(start_date)
end_datetime = np.datetime64(end_date)

# Slice the dataset to the desired time range
swe_data = dataset['sf'].sel(valid_time=slice(start_datetime, end_datetime))

# Step 4: Select the grid point for a specific latitude and longitude
latitude = 37.5  # Latitude of interest
longitude = 241.25  # Longitude of interest (Western Hemisphere)

# Nearest neighbor selection for the grid point
swe_at_point = swe_data.sel(latitude=latitude, longitude=longitude, method=None)

# Step 5: Convert SWE to Snowfall (using a ratio of 1:10)
snowfall_at_point = swe_at_point / 0.1  # Assuming 1mm of SWE equals 10mm of snowfall

# Step 6: Calculate cumulative snowfall
cumulative_snowfall = snowfall_at_point.cumsum(dim='valid_time')

# Step 7: Plot the cumulative snowfall time series at the specified grid point
plt.figure(figsize=(10, 6))
plt.plot(snowfall_at_point['valid_time'], cumulative_snowfall, marker='o', linestyle='-', color='b')
plt.title(f'Cumulative Snowfall at ({latitude}°N, {longitude-360}°W) from {start_date} to {end_date}', fontsize=14)
plt.xlabel('Time')
plt.ylabel('Cumulative Snowfall (m)')
plt.grid(True)

# Show plot
plt.show()
