# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 09:59:47 2024

@author: Nathaniel.George
"""

# Importing the necessary libraries to read a NetCDF File
import netCDF4 as nc  # Library for working with NetCDF files
import numpy as np    # Useful for numerical operations
import datetime       # To convert time units to a readable format


# Open the NetCDF file
# We need to use certain libraries to open NetCDF files. Below is what can be used to read a netcdf file
dataset = nc.Dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Temp2022SSW1.nc")

# Print basic metadata of the file
# General information about the file, like its global attributes (Lat, Lon, ect..)
print("NetCDF File Metadata:")
print(dataset)

# List all dimensions
print("\nDimensions:")
for dim in dataset.dimensions.items():
    print(f"{dim[0]}: {dim[1].size}")

# List all variables
# Determines what variables are in the file like temperature, pressure, wind, etc.
print("\nVariables:")
for var in dataset.variables.keys():
    print(var)

# Access the variable 'time' (if present) and inspect it
# We need to be able to determine why time format of the file so we can determine how to use it
if 'time' in dataset.variables:
    time_var = dataset.variables['time']
    print("\nTime Variable Details:")
    print(time_var)
    
    # Convert time values into a form that can be read
    try:
        time_units = time_var.units  # Time units (e.g., "days since 1970-01-01")
        time_values = time_var[:]  # The actual time data
        
        # Convert time to a readable format (assuming units are in a common format)
        time_readable = nc.num2date(time_values, time_units)
        print("\nReadable Time Data:")
        for t in time_readable:
            print(t)
    except Exception as e:
        print(f"Could not convert time: {e}")

# Close the NetCDF file
dataset.close()
