import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd 

# Load the NetCDF data file without decoding times
file_path = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Temp2022SSW1.nc")
data = xr.open_dataset(file_path, decode_times=False)

# Select the temperature data at 10 hPa and for latitudes between 90N to 60N
pressure_level = 10
latitude_range = slice(90, 60)
temperature_data = data.sel(pres=pressure_level, lat=latitude_range).geopClimMean

# Convert the time coordinates to datetime and filter by the required range
temperature_data['timeClim'] = xr.cftime_range(start='1970-01-01', periods=temperature_data.sizes['timeClim'], freq='D')
temperature_data = temperature_data.sel(timeClim=slice('2022-07-01', '2023-07-01'))

# Further reduce the size by subsetting the longitude dimension if necessary
temperature_data = temperature_data.isel(lon=slice(None, None, 10))  # Select every 10th longitude point

# Calculate daily mean temperature over the latitude and longitude range
daily_mean_temp = temperature_data.mean(dim=['lat', 'lon'])

# Calculate statistics for shading (min, max, percentiles)
mean_temp = daily_mean_temp.mean(dim='timeClim').compute()
max_temp = daily_mean_temp.max(dim='timeClim').compute()
min_temp = daily_mean_temp.min(dim='timeClim').compute()
percentiles = daily_mean_temp.quantile([0.1, 0.3, 0.7, 0.9], dim='timeClim').compute()

# Convert the filtered time data to day of year format for plotting
days_of_year = np.arange(len(daily_mean_temp.timeClim))

# Plot configuration
plt.figure(figsize=(14, 7))
plt.plot(days_of_year, daily_mean_temp.compute(), color='green', label='Daily Mean Temperature')
plt.plot(days_of_year, [mean_temp] * len(days_of_year), color='yellow', label='Mean Temperature', linewidth=2)
plt.fill_between(days_of_year, percentiles.sel(quantile=0.1), percentiles.sel(quantile=0.9), color='gray', alpha=0.3, label='10%-90% Range')
plt.fill_between(days_of_year, percentiles.sel(quantile=0.3), percentiles.sel(quantile=0.7), color='gray', alpha=0.6, label='30%-70% Range')
plt.fill_between(days_of_year, min_temp, max_temp, color='gray', alpha=0.1, label='Min-Max Range')

# Set x-axis range and labels for date visualization
plt.xlim(0, len(days_of_year)-1)
plt.xticks(np.arange(0, len(days_of_year), 30), labels=pd.date_range('2022-07-01', '2023-07-01', freq='M').strftime('%b %Y'))

# Plot details
plt.title(f'Temperature [K] at {pressure_level} hPa (90N to 60N)')
plt.xlabel('Date')
plt.ylabel('Temperature (K)')
plt.legend()
plt.grid(True)
plt.show()
