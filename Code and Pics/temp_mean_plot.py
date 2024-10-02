import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load the NetCDF data file (replace with your file path)
file_path = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\SSWC_v1.0_tempClimMean_ERAi_s19790101_e20140630_c20160701.nc")
data = xr.open_dataset(file_path, decode_times=False)

# Specify latitude and longitude for the location of interest
latitude = 48.5  # Replace with your latitude of interest
longitude = -101.5  # Replace with your longitude of interest
pressure_level = 10  # Replace with desired pressure level (e.g., 500 hPa)

# Select the temperature data for the specified location and pressure level
temperature_data = data.sel(lat=latitude, lon=longitude, pres=pressure_level, method='nearest').tempClimMean

# Print the selected data's dimensions and coordinates
print("Selected temperature data dimensions:", temperature_data.dims)
print("Selected temperature data coordinates:", temperature_data.coords)

# Ensure 'timeClim' is treated as a coordinate if not present
if 'timeClim' not in temperature_data.coords:
    temperature_data = temperature_data.assign_coords({'timeClim': data['timeClim']})

# Check if 'timeClim' is now in the coordinates
print("Updated temperature data coordinates:", temperature_data.coords)

# Check the shape and dimensions of the data before applying groupby
print("Shape of temperature data before groupby:", temperature_data.shape)

# Group by 'timeClim' and calculate the mean
try:
    daily_mean_temp = temperature_data.groupby('timeClim').mean()
    print("Successfully calculated daily mean temperature.")
except Exception as e:
    print("Error during groupby operation:", e)

# If successful, proceed with percentile calculations and plotting
if 'daily_mean_temp' in locals():
    try:
        # Manually compute percentiles since groupby with reduce caused issues
        # Calculate percentiles along the time dimension
        q10 = temperature_data.quantile(0.1, dim='timeClim')
        q30 = temperature_data.quantile(0.3, dim='timeClim')
        q70 = temperature_data.quantile(0.7, dim='timeClim')
        q90 = temperature_data.quantile(0.9, dim='timeClim')

        # Generate the plot
        plt.figure(figsize=(12, 6))
        plt.plot(daily_mean_temp['timeClim'], daily_mean_temp, label='Mean Daily Temperature', color='blue')
        plt.fill_between(daily_mean_temp['timeClim'], q10, q90, color='lightgray', alpha=0.5, label='10th-90th Percentile')
        plt.fill_between(daily_mean_temp['timeClim'], q30, q70, color='lightblue', alpha=0.5, label='30th-70th Percentile')

        # Plot details
        plt.title(f'Mean Daily Temperature at {latitude}°N, {longitude}°E and {pressure_level} hPa')
        plt.xlabel('Day of Year')
        plt.ylabel('Temperature (K)')  # Adjust units if necessary
        plt.legend()
        plt.grid(True)
        plt.show()

    except Exception as e:
        print("Error during percentile calculation or plotting:", e)