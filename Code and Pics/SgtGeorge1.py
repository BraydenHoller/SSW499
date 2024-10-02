#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 20:42:16 2024

@author: nathanielgeorge
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def calculate_mean_temp(file_path, lat_min, lat_max, pressure_level):
    # Open the NetCDF file
    dataset = nc.Dataset(file_path)
    
    # Extract variables
    lats = dataset.variables['lat'][:]
    pres = dataset.variables['pres'][:]
    times = dataset.variables['timeClim'][:]
    temps = dataset.variables['tempClimMean']
    
    # Find indices for the specified latitudes and pressure level
    lat_indices = np.where((lats >= lat_min) & (lats <= lat_max))[0]
    pres_index = np.where(pres == pressure_level)[0][0]

    # Initialize an array to store mean temperatures
    mean_temps = []

    # Iterate through each time step
    for t in range(len(times)):
        # Extract temperature data for the current time step, pressure level, and latitude range
        temp_slice = temps[t, pres_index, lat_indices, :]

        # Print the shape of temp_slice to see how many points are used
        print(f"Time step {t}: temp_slice shape = {temp_slice.shape}")

        # Calculate the mean temperature across the latitude and longitude dimensions
        mean_temp = np.mean(temp_slice)
        mean_temps.append(mean_temp)

    # Convert the times to standard datetime objects
    time_units = dataset.variables['timeClim'].units
    time_calendar = dataset.variables['timeClim'].calendar if 'calendar' in dataset.variables['timeClim'].ncattrs() else 'standard'
    times_converted = nc.num2date(times, units=time_units, calendar=time_calendar)

    # Convert cftime objects to standard datetime objects
    times_converted = [datetime(t.year, t.month, t.day) for t in times_converted]

    # Close the dataset
    dataset.close()
    
    return times_converted, mean_temps

# Define the file path, latitude range, and pressure level
file_path = r"/Users/nathanielgeorge/Documents/SSWC_v1.0_tempClimMean_MERRA2_s19800101_e20140630_c20160701.nc"  # Replace with your file path
lat_min, lat_max = -90, 90
pressure_level = 500  # Replace with your desired pressure level in hPa

# Calculate the mean temperature for each time step
times, mean_temps = calculate_mean_temp(file_path, lat_min, lat_max, pressure_level)

#%%
import netCDF4 as nc

def count_grid_points(file_path):
    # Open the NetCDF file
    dataset = nc.Dataset(file_path)
    
    # Extract dimensions
    lat_dim = dataset.dimensions['lat'].size
    lon_dim = dataset.dimensions['lon'].size
    pres_dim = dataset.dimensions['pres'].size
    time_dim = dataset.dimensions['timeClim'].size
    
    # Calculate the total number of grid points
    total_grid_points = lat_dim * lon_dim
    total_with_pres = total_grid_points * pres_dim
    total_with_time = total_with_pres * time_dim
    
    # Print out the information
    print(f"Latitude points: {lat_dim}")
    print(f"Longitude points: {lon_dim}")
    print(f"Pressure levels: {pres_dim}")
    print(f"Time steps: {time_dim}")
    print(f"Total number of grid points (lat x lon): {total_grid_points}")
    print(f"Total number of grid points (lat x lon x pres): {total_with_pres}")
    print(f"Total number of grid points (lat x lon x pres x time): {total_with_time}")
    
    # Close the dataset
    dataset.close()

# Define the file path
file_path = r"/Users/nathanielgeorge/Documents/SSWC_v1.0_tempClimMean_MERRA2_s19800101_e20140630_c20160701.nc"  # Replace with your file path

# Count the grid points
count_grid_points(file_path)