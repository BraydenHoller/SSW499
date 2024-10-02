# import netCDF4 as nc
# import xarray as xr
# import matplotlib.pyplot as plt
# from datetime import datetime

# # Step 1: Load the dataset
# file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Snow_2223.nc"  
# dataset = xr.open_dataset(file_path)
# times = dataset.variables['valid_time'][:]  # Assuming this is in hourly intervals

# # # Step 2: Convert UNIX time (valid_time) to datetime if needed
# # time_units = "seconds since 1970-01-01"
# # time_calendar = "proleptic_gregorian"
# # dataset['valid_time'] = xr.coding.times.decode_cf_datetime(dataset['valid_time'], units=time_units, calendar=time_calendar)

# # # Step 3: Select the data between 2022-07-01 and 2023-07-01
# # start_date = "2022-07-01"
# # end_date = "2023-07-01"
# # snowfall_data = dataset['sf'].sel(valid_time=slice(start_date, end_date))

# time_units = dataset.variables['valid_time'].units
# time_calendar = dataset.variables['valid_time'].calendar
# times_converted = nc.num2date(times, units=time_units, calendar=time_calendar)

# # Convert cftime.DatetimeProlepticGregorian to native datetime.datetime
# times_converted = [datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in times_converted]

# start_date = "2022-07-01"
# end_date = "2023-07-01"

# # Filter time steps based on the desired date range
# start_datetime = datetime.strptime(start_date, "%Y-%m-%d")
# end_datetime = datetime.strptime(end_date, "%Y-%m-%d")

# snowfall_data = dataset['sf'].sel(valid_time=slice(start_date, end_date))

# # Step 4: Select the grid point for a specific latitude and longitude
# latitude = 61.25  # Update to the specific latitude
# longitude = -149.5  # Update to the specific longitude

# # Nearest neighbor selection for the grid point
# snowfall_at_point = snowfall_data.sel(latitude=latitude, longitude=longitude, method='nearest')

# # Step 5: Plot the snowfall time series at the specified grid point
# plt.figure(figsize=(10, 6))
# plt.plot(snowfall_at_point['valid_time'], snowfall_at_point, marker='o', linestyle='-', color='b')
# plt.title(f'Snowfall at ({latitude}°N, {149.5}°W) from {start_date} to {end_date}', fontsize=14)
# plt.xlabel('Time')
# plt.ylabel('Snowfall (m)')
# plt.grid(True)

# # Show plot
# plt.show()
# import xarray as xr
# import matplotlib.pyplot as plt
# import numpy as np
# from datetime import datetime

# # Step 1: Load the dataset
# file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Snow_2223.nc"
# dataset = xr.open_dataset(file_path)

# # Step 2: Verify the 'valid_time' coordinate and use it directly if it’s already in datetime format
# times_converted = dataset['valid_time'].values  # Get the time values directly

# # Step 3: Select the data between 2022-07-01 and 2023-07-01
# start_date = "2022-07-01"
# end_date = "2023-07-01"
# start_datetime = np.datetime64(start_date)
# end_datetime = np.datetime64(end_date)

# # Slice the dataset to the desired time range
# snowfall_data = dataset['sf'].sel(valid_time=slice(start_datetime, end_datetime))

# # Step 4: Select the grid point for a specific latitude and longitude
# latitude = 61.25  # Latitude of interest
# longitude = -149.5  # Longitude of interest (Western Hemisphere)

# # Nearest neighbor selection for the grid point
# snowfall_at_point = snowfall_data.sel(latitude=latitude, longitude=longitude, method='nearest')

# # Step 5: Calculate cumulative snowfall
# cumulative_snowfall = snowfall_at_point.cumsum(dim='valid_time')

# # Step 6: Plot the cumulative snowfall time series at the specified grid point
# plt.figure(figsize=(10, 6))
# plt.plot(snowfall_at_point['valid_time'], 100* cumulative_snowfall, marker='o', linestyle='-', color='b') # Including 10:1 SWE ratio
# plt.title(f'Cumulative Snowfall at ({latitude}°N, {abs(longitude)}°W) from {start_date} to {end_date}', fontsize=14)
# plt.xlabel('Time')
# plt.ylabel('Cumulative Snowfall (m)')
# plt.grid(True)

# # Show plot
# plt.show()

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

# Step 1: Load the dataset
file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Snow_2223.nc"
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
latitude = 61.25  # Latitude of interest
longitude = -149.5  # Longitude of interest (Western Hemisphere)

# Nearest neighbor selection for the grid point
swe_at_point = swe_data.sel(latitude=latitude, longitude=longitude, method='nearest')

# Step 5: Convert SWE to Snowfall (using a ratio of 1:10)
snowfall_at_point = swe_at_point / 0.1  # Assuming 1mm of SWE equals 10mm of snowfall

# Step 6: Calculate cumulative snowfall
cumulative_snowfall = snowfall_at_point.cumsum(dim='valid_time')

# Step 7: Plot the cumulative snowfall time series at the specified grid point
plt.figure(figsize=(10, 6))
plt.plot(snowfall_at_point['valid_time'], cumulative_snowfall, marker='o', linestyle='-', color='b')
plt.title(f'Cumulative Snowfall at ({latitude}°N, {abs(longitude)}°W) from {start_date} to {end_date}', fontsize=14)
plt.xlabel('Time')
plt.ylabel('Cumulative Snowfall (m)')
plt.grid(True)

# Show plot
plt.show()
