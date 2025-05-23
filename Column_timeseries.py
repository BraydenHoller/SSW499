import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def calculate_mean_wind(file_path, lat_min, lat_max, pressure_levels, start_date, end_date):
    try:
        # Open the NetCDF file
        print("Opening NetCDF file...")
        dataset = nc.Dataset(file_path)
        print("File opened successfully.")
        
        # Extract variables
        lats = dataset.variables['latitude'][:]
        pres = dataset.variables['pressure_level'][:]
        times = dataset.variables['valid_time'][:]  # Assuming this is in hourly intervals
        U_winds = dataset.variables['u']  # Don't load the entire U winds array
        V_winds = dataset.variables['v']  # Don't load the entire V winds array
        
        # Find indices for the specified latitudes
        lat_indices = np.where((lats >= lat_min) & (lats <= lat_max))[0]

        # Convert all time steps to datetime objects
        time_units = dataset.variables['valid_time'].units
        time_calendar = dataset.variables['valid_time'].calendar
        times_converted = nc.num2date(times, units=time_units, calendar=time_calendar)

        # Convert cftime.DatetimeProlepticGregorian to native datetime.datetime
        times_converted = [datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in times_converted]

        # Filter time steps based on the desired date range
        start_datetime = datetime.strptime(start_date, "%Y-%m-%d")
        end_datetime = datetime.strptime(end_date, "%Y-%m-%d")
        time_indices = [i for i, t in enumerate(times_converted) if start_datetime <= t < end_datetime]

        if not time_indices:
            print(f"No data found between {start_date} and {end_date}.")
            return [], {}

        # Initialize a dictionary to store mean winds for each pressure level
        mean_winds_by_level = {}

        # Iterate through the selected pressure levels
        for pressure_level in pressure_levels:
            print(f"Processing pressure level: {pressure_level} hPa")
            pres_index = np.where(pres == pressure_level)[0][0]  # Find the index for the given pressure level

            mean_winds = []

            # Iterate through the filtered time steps
            for t in time_indices:
                # Extract wind data for the current time step, pressure level, and latitude range
                U_wind_slice = U_winds[t, pres_index, lat_indices, :]  # Lazy slice
                V_wind_slice = V_winds[t, pres_index, lat_indices, :]  # Lazy slice
                mean_wind = np.sqrt(U_wind_slice**2 + V_wind_slice**2).mean()
                mean_winds.append(mean_wind)

            # Store the results for this pressure level
            mean_winds_by_level[pressure_level] = mean_winds

        # Select only the filtered times
        selected_times = [times_converted[t] for t in time_indices]

        print("Time conversion successful.")
        
        # Close the dataset
        dataset.close()
        
        return selected_times, mean_winds_by_level
    except Exception as e:
        print(f"An error occurred: {e}")
        return [], {}

def plot_mean_wind(times, mean_winds_by_level, lat_min, lat_max):
    if len(times) == 0 or len(mean_winds_by_level) == 0:
        print("No data to plot.")
        return
    
    # Plot the mean winds over time for each pressure level
    plt.figure(figsize=(10, 6))

    for pressure_level, mean_winds in mean_winds_by_level.items():
        plt.plot(times, mean_winds, marker='o', label=f'{pressure_level} hPa')

    # Customize plot
    plt.title(f'Mean Winds ({lat_min}-{lat_max}°N)')
    plt.xlabel('Time')
    plt.ylabel('Wind Speed (m/s)')
    plt.grid(True)
    plt.legend()
    plt.show()

# Define the file path, latitude range, and pressure levels
file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\U_V_Wind_2223.nc"  # Replace with your file path
lat_min, lat_max = 60, 90
pressure_levels = [1, 10, 300, 500, 700, 850, 1000]  # List of desired pressure levels in hPa

# Specify the date range
start_date = "2022-07-01"
end_date = "2023-07-01"

# Calculate the mean wind for each day in the specified date range and for the specified pressure levels
times, mean_winds_by_level = calculate_mean_wind(file_path, lat_min, lat_max, pressure_levels, start_date, end_date)

# Plot the results
plot_mean_wind(times, mean_winds_by_level, lat_min, lat_max)
