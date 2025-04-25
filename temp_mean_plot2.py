import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import cftime

# Load the NetCDF data file (replace with your file path)
file_path = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\SSWC_v1.0_tempClimMean_ERAi_s19790101_e20140630_c20160701.nc")
data = xr.open_dataset(file_path, decode_times=False)

file_path2 = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Temp202223SSW.nc")
data2 = xr.open_dataset(file_path2, decode_times=False)

# Specify latitude and longitude for the location of interest
latitude = 48.5  # Replace with your latitude of interest
longitude = -101.5  # Replace with your longitude of interest
pressure_level = 1000  # Replace with desired pressure level (e.g., 500 hPa)

# Select the temperature data for the specified location and pressure level
temperature_data = data.sel(lat=latitude, lon=longitude, pres=pressure_level, method='nearest').tempClimMean
temperature_data2 = data2.sel(latitude=latitude, longitude=longitude, pressure_level=pressure_level, method='nearest').t

# Manually set the 'valid_time' coordinate using cftime_range
temperature_data2['valid_time'] = xr.cftime_range(start='2022-07-01', periods=temperature_data2.sizes['valid_time'], freq='H', calendar='proleptic_gregorian')
temperature_data2 = temperature_data2.sel(valid_time=slice('2022-07-01', '2023-07-01'))

# Check the structure of the time variable
print("Original time variable structure:")
print(temperature_data['timeClim'])
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt
# from pathlib import Path
# import pandas as pd
# import cftime

# # Load the NetCDF data file (replace with your file path)
# file_path = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\SSWC_v1.0_tempClimMean_ERAi_s19790101_e20140630_c20160701.nc")
# data = xr.open_dataset(file_path, decode_times=False)

# file_path2 = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Temp202223SSW.nc")
# data2 = xr.open_dataset(file_path2, decode_times=False)

# # Specify latitude and longitude for the location of interest
# latitude = 48.5  # Replace with your latitude of interest
# longitude = -101.5  # Replace with your longitude of interest
# pressure_level = 1000  # Replace with desired pressure level (e.g., 500 hPa)

# # Select the temperature data for the specified location and pressure level
# temperature_data = data.sel(lat=latitude, lon=longitude, pres=pressure_level, method='nearest').tempClimMean
# temperature_data2 = data2.sel(latitude=latitude, longitude=longitude, pressure_level=pressure_level, method='nearest').t

# # Manually set the 'valid_time' coordinate using cftime_range
# temperature_data2['valid_time'] = xr.cftime_range(start='2022-07-01', periods=temperature_data2.sizes['valid_time'], freq='H', calendar='proleptic_gregorian')
# temperature_data2 = temperature_data2.sel(valid_time=slice('2022-07-01', '2023-07-01'))

# # Define the reference date (2023-01-01)
# reference_date = cftime.DatetimeProlepticGregorian(2023, 1, 1)

# # Calculate float days from the reference date
# # Dates before 2023-01-01 will be negative, and dates after will be positive
# float_days_from_ref = [
#     (t - reference_date).days + (t - reference_date).seconds / 86400 for t in temperature_data2['valid_time'].values
# ]

# # Replace the 'valid_time' values with the converted float days
# temperature_data2 = temperature_data2.assign_coords(valid_time=float_days_from_ref)

# # Check the updated structure of the time variable
# print("Updated time variable with float values (negative for before 2023-01-01, positive for after):")
# print(temperature_data2['valid_time'])

# # Check for missing values
# print("Number of missing values before handling:", temperature_data.isnull().sum().values)
# print("Number of missing values before handling:", temperature_data2.isnull().sum().values)

# # Strategy: Interpolation (or choose another strategy as needed)
# temperature_data_clean = temperature_data.interpolate_na(dim='timeClim', method='linear')
# temperature_data_clean2 = temperature_data2.interpolate_na(dim='valid_time', method='linear')

# # Check for missing values after handling
# print("Number of missing values after handling:", temperature_data_clean.isnull().sum().values)
# print("Number of missing values after handling:", temperature_data_clean2.isnull().sum().values)

# # Check the dimensions of temperature_data_clean
# print("Dimensions of 'temperature_data_clean':", temperature_data_clean.dims)
# print("Dimensions of 'temperature_data_clean':", temperature_data_clean2.dims)

# # Check if 'lat' and 'lon' are still present, and calculate mean if they are
# if 'lat' in temperature_data_clean.dims and 'lon' in temperature_data_clean.dims:
#     daily_mean_temp = temperature_data_clean.mean(dim=['lat', 'lon'])
# else:
#     # If 'lat' and 'lon' are not present, assume reduction has already happened
#     daily_mean_temp = temperature_data_clean

# # Check if 'lat' and 'lon' are still present, and calculate mean if they are
# if 'latitude' in temperature_data_clean2.dims and 'longitude' in temperature_data_clean2.dims:
#     daily_mean_temp2 = temperature_data_clean2.mean(dim=['latitude', 'longitude'])
# else:
#     # If 'lat' and 'lon' are not present, assume reduction has already happened
#     daily_mean_temp2 = temperature_data_clean2    

# # Calculate percentiles over the timeClim dimension directly
# q10 = daily_mean_temp.quantile(0.1, dim='timeClim')
# q30 = daily_mean_temp.quantile(0.3, dim='timeClim')
# q70 = daily_mean_temp.quantile(0.7, dim='timeClim')
# q90 = daily_mean_temp.quantile(0.9, dim='timeClim')

# # Resample temperature_data2 to daily means
# temperature_data2_daily = temperature_data_clean2.resample(valid_time='1D').mean()

# # Plotting the results after resampling to daily means
# plt.figure(figsize=(12, 6))

# # Plot the daily mean temperature from temperature_data
# plt.plot(daily_mean_temp['timeClim'], daily_mean_temp, label='Mean Daily Temperature (1979-2014)', color='yellow')

# # Plot the resampled daily mean temperature from temperature_data2
# plt.plot(temperature_data2_daily['valid_time'], temperature_data2_daily, label='Mean Daily Temperature (2022-2023)', color='red')

# # Add percentiles for temperature_data
# plt.fill_between(daily_mean_temp['timeClim'], q10, q90, color='lightgray', alpha=0.5, label='10th-90th Percentile (1979-2014)')
# plt.fill_between(daily_mean_temp['timeClim'], q30, q70, color='lightblue', alpha=0.5, label='30th-70th Percentile (1979-2014)')

# # Plot details
# plt.title(f'Mean Daily Temperature at {latitude}°N, {longitude}°W and {pressure_level} hPa')
# plt.xlabel('Day of Year')
# plt.ylabel('Temperature (K)')  # Adjust units if necessary
# plt.legend()
# plt.grid(True)
# plt.show()