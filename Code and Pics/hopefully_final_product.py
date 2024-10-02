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
temperature_data2['valid_time'] = xr.cftime_range(start='2022-01-01', periods=temperature_data2.sizes['valid_time'], freq='H', calendar='proleptic_gregorian')
temperature_data2 = temperature_data2.sel(valid_time=slice('2022-07-01', '2023-07-01'))

# Check the structure of the time variable
print("Original time variable structure:")
print(temperature_data2['valid_time'])

# Manually set calendar and units for the time variable if not already present
temperature_data2['valid_time'].attrs['units'] = 'days since 1980-01-01'
temperature_data2['valid_time'].attrs['calendar'] = 'julian'

# Check the updated attributes
time_attrs = temperature_data2['valid_time'].attrs
print("Updated valid_time variable attributes:", time_attrs)

# Define a new reference date for conversion (e.g., 2000-01-01)
new_reference_date = cftime.DatetimeProlepticGregorian(2023, 1, 1)

# Convert valid_time values to the same calendar as the reference date
converted_valid_time = [
    t if t.calendar == new_reference_date.calendar else cftime.DatetimeProlepticGregorian(t.year, t.month, t.day, t.hour, t.minute, t.second) 
    for t in temperature_data2['valid_time'].values
]

# Calculate float days from the new reference date
float_days_from_new_ref = [
    (t - new_reference_date).days + (t - new_reference_date).seconds / 86400 for t in converted_valid_time
]

# Replace the 'valid_time' values with the converted float days
temperature_data2 = temperature_data2.assign_coords(valid_time=float_days_from_new_ref)

# Check the updated structure of the time variable
print("Updated time variable with float days:")
print(temperature_data2['valid_time'])

# Check for missing values
print("Number of missing values before handling:", temperature_data.isnull().sum().values)
print("Number of missing values before handling:", temperature_data2.isnull().sum().values)

# Strategy: Interpolation (or choose another strategy as needed)
temperature_data_clean = temperature_data.interpolate_na(dim='timeClim', method='linear')
temperature_data_clean2 = temperature_data2.interpolate_na(dim='valid_time', method='linear')

# Check for missing values after handling
print("Number of missing values after handling:", temperature_data_clean.isnull().sum().values)
print("Number of missing values after handling:", temperature_data_clean2.isnull().sum().values)

# Check the dimensions of temperature_data_clean
print("Dimensions of 'temperature_data_clean':", temperature_data_clean.dims)
print("Dimensions of 'temperature_data_clean':", temperature_data_clean2.dims)

# Check if 'lat' and 'lon' are still present, and calculate mean if they are
if 'lat' in temperature_data_clean.dims and 'lon' in temperature_data_clean.dims:
    daily_mean_temp = temperature_data_clean.mean(dim=['lat', 'lon'])
else:
    # If 'lat' and 'lon' are not present, assume reduction has already happened
    daily_mean_temp = temperature_data_clean

# Check if 'lat' and 'lon' are still present, and calculate mean if they are
if 'latitude' in temperature_data_clean2.dims and 'longitude' in temperature_data_clean2.dims:
    daily_mean_temp2 = temperature_data_clean2.mean(dim=['latitude', 'longitude'])
else:
    # If 'lat' and 'lon' are not present, assume reduction has already happened
    daily_mean_temp2 = temperature_data_clean2    

# Calculate percentiles over the timeClim dimension directly
# Calculate percentiles for each day along the curve using apply_ufunc
def calculate_percentiles(array, q):
    return np.percentile(array, q, axis=0)

percentiles = xr.apply_ufunc(
    calculate_percentiles,
    daily_mean_temp,
    input_core_dims=[['timeClim']],
    output_core_dims=[['percentile']],
    vectorize=True,
    kwargs={'q': [10, 30, 70, 90]},
    dask='parallelized',
    output_dtypes=[float]
)

# Extract individual percentiles from the result
q10, q30, q70, q90 = percentiles.sel(percentile=0), percentiles.sel(percentile=1), percentiles.sel(percentile=2), percentiles.sel(percentile=3)

count = 0
total = 0
dims=("time")
coords={
    "time": pd.RangeIndex(-184.0, 182.0, 1, dtype= int)
}
empty_array = xr.DataArray(np.zeros(366), coords=coords, dims=dims)
empty_list = []
for t in daily_mean_temp2.values:
    total = total + t
    count = count + 1
    if count == 24:
        avg = total / 24
        empty_list.append(avg)
        count = 0
        total = 0
    else:
        pass

empty_array += empty_list

# Plotting the results
plt.figure(figsize=(12, 6))
plt.plot(daily_mean_temp['timeClim'], daily_mean_temp, label='Mean Daily Temperature (Climatology)', color='yellow', linewidth=2)
plt.plot(empty_array['time'], empty_array, label='Mean Daily Temperature (2022-2023)', color='red', linewidth=2)
plt.fill_between(daily_mean_temp['timeClim'], q10, q90, color='lightgray', alpha=0.5, label='10th-90th Percentile (Climatology)')
plt.fill_between(daily_mean_temp['timeClim'], q30, q70, color='lightblue', alpha=0.5, label='30th-70th Percentile (Climatology)')

# Plot details
plt.title(f'Mean Daily Temperature at {latitude}°N, 101.5°W and {pressure_level} hPa')
plt.xlabel('Day of Year Measured from Jan 01')
plt.ylabel('Temperature (K)')  # Adjust units if necessary

# Enhanced legend
plt.legend(title='Temperature Data', title_fontsize='13', fontsize='11', shadow=True, facecolor='white', edgecolor='black', loc='upper left', frameon=True)
plt.grid(True)
plt.show()
