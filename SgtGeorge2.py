import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def calculate_mean_temp(file_path, lat_min, lat_max, pressure_level):
    # Open the NetCDF file
    dataset = nc.Dataset(file_path)
    
    # Extract variables
    lats = dataset.variables['latitude'][:]
    pres = dataset.variables['pressure_level'][:]
    times = dataset.variables['valid_time'][:]
    temps = dataset.variables['t']
    
    # Find indices for the specified latitudes and pressure level
    lat_indices = np.where((lats >= lat_min) & (lats <= lat_max))[0]
    pres_index = np.where(pres == pressure_level)[0][0]

    # Initialize an array to store mean temperatures
    mean_temps = []

    # Iterate through each time step
    for t in range(len(times)):
        # Extract temperature data for the current time step, pressure level, and latitude range
        temp_slice = temps[t, pres_index, lat_indices, :]
        
        # Calculate the mean temperature across the latitude and longitude dimensions
        mean_temp = np.mean(temp_slice)
        mean_temps.append(mean_temp)

    # Convert the times to standard datetime objects
    time_units = dataset.variables['valid_time'].units
    time_calendar = dataset.variables['valid_time'].calendar if 'calendar' in dataset.variables['valid_time'].ncattrs() else 'standard'
    times_converted = nc.num2date(times, units=time_units, calendar=time_calendar)

    # Convert cftime objects to standard datetime objects
    times_converted = [datetime(t.year, t.month, t.day) for t in times_converted]

    # Close the dataset
    dataset.close()
    
    return times_converted, mean_temps

def plot_mean_temp(times, mean_temps, lat_min, lat_max, pressure_level):
    # Plot the mean temperature over time
    plt.figure(figsize=(10, 6))
    plt.plot(times, mean_temps, marker='o')
    plt.title(f'Mean Temperature ({lat_min}-{lat_max}°N, {pressure_level} hPa)')
    plt.xlabel('Time')
    plt.ylabel('Temperature (K)')
    plt.grid(True)
    plt.show()

# Define the file path, latitude range, and pressure level
file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Temp2233.nc"  # Replace with your file path
lat_min, lat_max = 60, 90
pressure_level = 1  # Replace with your desired pressure level in hPa

# Calculate the mean temperature for each time step
times, mean_temps = calculate_mean_temp(file_path, lat_min, lat_max, pressure_level)

# Plot the results
plot_mean_temp(times, mean_temps, lat_min, lat_max, pressure_level)
