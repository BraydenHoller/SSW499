# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 09:59:47 2024

@author: Nathaniel.George
"""

'''This section will show you what time steps and levels are in the file which you will use below in order to create a plot'''

# Import necessary libraries
import netCDF4 as nc
from netCDF4 import num2date  # For converting time units to readable format

# Open the NetCDF file
dataset = nc.Dataset(r"C:\Users\Nathaniel.George\OneDrive - afacademy.af.edu\Desktop\SSW\SSWC_v1.0_geopClimMean_ERA40_s19570901_e20020630_c20160701 - Copy.nc")

# Extract and print all available pressure levels
pressure_levels = dataset.variables['pres'][:]  # Pressure levels
print("Pressure Levels (hPa):")
print(pressure_levels)

# Extract and print all available time steps in human-readable format
time_var = dataset.variables['timeClim'][:]  # Time values
time_units = dataset.variables['timeClim'].units  # Time units (e.g., "days since 1970-01-01")
time_readable = num2date(time_var, time_units)  # Convert time to datetime objects

print("\nTime Steps (Converted to Date):")
for t in time_readable:
    print(t.strftime('%Y-%m-%d %H:%M:%S'))  # Print each time step in human-readable format

# Close the dataset
dataset.close()



#%%

'''This section of the code will allow you to plot the netCDF file. You will need to know the variables when exploring the file in order to plot. Additionally you will need to make adjustments to the time step and levels below'''

# Import necessary libraries
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs  # For geographic projections
import cartopy.feature as cfeature
from netCDF4 import num2date  # For converting time units to readable format

'''You need to select the time step and pressure level you want. The portion of the code at the top will go through and tell you what time steps are available and the different pressure levels. The title of the code will also display this'''
# Define indices for time step and pressure level
time_step_index = 100  # Set your desired time step index here
pressure_level_index = 15  # Set your desired pressure level index here

'''You can add contours to the plot to help get a better visualization. BUT the code will do its best to do the intervals equally based on how many intervals you want to see. You will have to play around with this as it can make the map very messy'''
# Define contour intervals
interval = 200  # Set your desired interval here
num_intervals = 10  # Number of intervals (used for dynamic intervals)

# Define latitude and longitude range or indices
lat_range = slice(-30, 60)  # Adjust as needed (e.g., latitudes from index 10 to 50)
lon_range = slice(-180, 180)  # Adjust as needed (e.g., longitudes from index 10 to 50)

# Open the NetCDF file
dataset = nc.Dataset(r"C:\Users\Nathaniel.George\OneDrive - afacademy.af.edu\Desktop\SSW\SSWC_v1.0_geopClimMean_ERA40_s19570901_e20020630_c20160701 - Copy.nc")

# Extract latitude, longitude, and geopotential climatological mean
lat = dataset.variables['lat'][:]  # Latitude values
lon = dataset.variables['lon'][:]  # Longitude values
geopClimMean = dataset.variables['geopClimMean'][:, :, :, :]  # All time steps and pressure levels

# Extract time variable and convert to readable format
time_var = dataset.variables['timeClim'][:]  # Time values
time_units = dataset.variables['timeClim'].units  # Time units (e.g., "days since 1970-01-01")
time_readable = num2date(time_var, time_units)  # Convert time to datetime objects

# Extract pressure levels and their attributes
pres_var = dataset.variables['pres']
pressure_levels = pres_var[:]  # Pressure levels
pressure_units = pres_var.units if hasattr(pres_var, 'units') else 'unknown'  # Get units if available

# Use the specified time step and pressure level for plotting
plot_time = time_readable[time_step_index].strftime('%Y-%m-%d %H:%M:%S')  # Convert time to string
selected_pressure = pressure_levels[pressure_level_index]  # Get the pressure level at the specified index

# Extract the subset of data for the specified time and pressure level
subset_geopClimMean = geopClimMean[time_step_index, pressure_level_index, lat_range, lon_range]
lat_subset = lat[lat_range]
lon_subset = lon[lon_range]

# Create a meshgrid for latitude and longitude
lon2d, lat2d = np.meshgrid(lon_subset, lat_subset)

# Plotting the data
plt.figure(figsize=(12, 8))

# Create a Cartopy projection for the map (PlateCarree for global projection)
ax = plt.axes(projection=ccrs.PlateCarree())

# Add map features (coastlines, borders, etc.)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')

# Plot the geopotential climatological mean data on the map
# Use pcolormesh for colored grid cells
pcolor = plt.pcolormesh(lon2d, lat2d, subset_geopClimMean, transform=ccrs.PlateCarree(), cmap='coolwarm')

# Define contour intervals
min_value = np.min(subset_geopClimMean)
max_value = np.max(subset_geopClimMean)
dynamic_interval = (max_value - min_value) / num_intervals  # Dynamic interval calculation

# Generate contour levels based on calculated intervals
levels = np.arange(min_value, max_value + dynamic_interval, dynamic_interval)
contour = plt.contour(lon2d, lat2d, subset_geopClimMean, levels=levels, colors='black', transform=ccrs.PlateCarree())

# Add contour labels if needed
plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f')

# Set the map extent based on the specified latitude and longitude range
ax.set_extent([lon_subset.min(), lon_subset.max(), lat_subset.min(), lat_subset.max()])

# Add title and labels, including date, time, and pressure level information
plt.title(f"Geopotential Climatological Mean\nDate: {plot_time} | Pressure Level: {selected_pressure} {pressure_units}")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

# Add a horizontal colorbar below the plot
cbar = plt.colorbar(pcolor, orientation='horizontal', pad=0.1)
cbar.set_label("Geopotential Climatological Mean")

# Show the plot
plt.show()

# Close the dataset
dataset.close()



