# import netCDF4 as nc
# import numpy as np
# import matplotlib.pyplot as plt
# from datetime import datetime

# def calculate_mean_temp(file_path, lat_min, lat_max, pressure_level, start_date, end_date):
#     try:
#         # Open the NetCDF file
#         print("Opening NetCDF file...")
#         dataset = nc.Dataset(file_path)
#         print("File opened successfully.")
        
#         # Extract variables
#         lats = dataset.variables['latitude'][:]
#         pres = dataset.variables['pressure_level'][:]
#         times = dataset.variables['valid_time'][:]  # Assuming this is in hourly intervals
#         temps = dataset.variables['t']  # Don't load the entire temperature array
        
#         # Handle missing values if any
#         if '_FillValue' in dataset.variables['t'].ncattrs():
#             fill_value = dataset.variables['t']._FillValue
#             temps.set_auto_maskandscale(True)  # Automatically mask and scale missing values

#         print(f"Data extracted: {len(times)} time steps.")
        
#         # Find indices for the specified latitudes and pressure level
#         lat_indices = np.where((lats >= lat_min) & (lats <= lat_max))[0]
#         pres_index = np.where(pres == pressure_level)[0][0]

#         # Convert all time steps to datetime objects
#         time_units = dataset.variables['valid_time'].units
#         time_calendar = dataset.variables['valid_time'].calendar
#         times_converted = nc.num2date(times, units=time_units, calendar=time_calendar)

#         # Convert cftime.DatetimeProlepticGregorian to native datetime.datetime
#         times_converted = [datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in times_converted]

#         # Filter time steps based on the desired date range
#         start_datetime = datetime.strptime(start_date, "%Y-%m-%d")
#         end_datetime = datetime.strptime(end_date, "%Y-%m-%d")
#         time_indices = [i for i, t in enumerate(times_converted) if start_datetime <= t < end_datetime]

#         if not time_indices:
#             print(f"No data found between {start_date} and {end_date}.")
#             return [], []

#         print(f"Processing {len(time_indices)} time steps within the date range {start_date} to {end_date}.")

#         # Initialize an array to store mean temperatures
#         mean_temps = []
        
#         # Iterate through the filtered time steps (one time step per day)
#         for i in range(0, len(time_indices)):
#             t = time_indices[i]
#             print(f"Processing time step {t} (Date: {times_converted[t]})")
            
#             # Extract temperature data for the current time step, pressure level, and latitude range
#             temp_slice = temps[t, pres_index, lat_indices, :]  # Lazy slice
#             mean_temp = np.mean(temp_slice)
#             mean_temps.append(mean_temp)

#         # Select only the filtered times
#         selected_times = [times_converted[t] for t in time_indices]

#         print("Time conversion successful.")
        
#         # Close the dataset
#         dataset.close()
        
#         return selected_times, mean_temps
#     except Exception as e:
#         print(f"An error occurred: {e}")
#         return [], []

# def plot_mean_temp(times, mean_temps, lat_min, lat_max, pressure_level):
#     if len(times) == 0 or len(mean_temps) == 0:
#         print("No data to plot.")
#         return
    
#     # Plot the mean temperature over time
#     plt.figure(figsize=(10, 6))
#     plt.plot(times, mean_temps, marker='o')
#     plt.title(f'Mean Temperature ({lat_min}-{lat_max}°N, {pressure_level} hPa)')
#     plt.xlabel('Time')
#     plt.ylabel('Temperature (K)')
#     plt.grid(True)
#     plt.show()

# # Define the file path, latitude range, and pressure level
# file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Temp2233.nc"  # Replace with your file path
# lat_min, lat_max = 60, 90
# pressure_level = 10  # Replace with your desired pressure level in hPa

# # Specify the date range
# start_date = "2022-07-01"
# end_date = "2023-07-01"

# # Calculate the mean temperature for each day in the specified date range
# times1, mean_temps1 = calculate_mean_temp(file_path, lat_min, lat_max, 1, start_date, end_date)
# times2, mean_temps2 = calculate_mean_temp(file_path, lat_min, lat_max, 10, start_date, end_date)
# times3, mean_temps3 = calculate_mean_temp(file_path, lat_min, lat_max, 300, start_date, end_date)
# times4, mean_temps4 = calculate_mean_temp(file_path, lat_min, lat_max, 500, start_date, end_date)
# times5, mean_temps5 = calculate_mean_temp(file_path, lat_min, lat_max, 700, start_date, end_date)
# times6, mean_temps6 = calculate_mean_temp(file_path, lat_min, lat_max, 850, start_date, end_date)
# times7, mean_temps7 = calculate_mean_temp(file_path, lat_min, lat_max, 1000, start_date, end_date)
# # Plot the results
# plot_mean_temp(times1, mean_temps1, lat_min, lat_max, 1)
# plot_mean_temp(times2, mean_temps2, lat_min, lat_max, 10)
# plot_mean_temp(times3, mean_temps3, lat_min, lat_max, 300)
# plot_mean_temp(times4, mean_temps4, lat_min, lat_max, 500)
# plot_mean_temp(times5, mean_temps5, lat_min, lat_max, 700)
# plot_mean_temp(times6, mean_temps6, lat_min, lat_max, 850)
# plot_mean_temp(times7, mean_temps7, lat_min, lat_max, 1000)
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def calculate_mean_temp(file_path, lat_min, lat_max, pressure_levels, start_date, end_date):
    try:
        # Open the NetCDF file
        print("Opening NetCDF file...")
        dataset = nc.Dataset(file_path)
        print("File opened successfully.")
        
        # Extract variables
        lats = dataset.variables['latitude'][:]
        pres = dataset.variables['pressure_level'][:]
        times = dataset.variables['valid_time'][:]  # Assuming this is in hourly intervals
        temps = dataset.variables['t']  # Don't load the entire U temps array
        
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

        # Initialize a dictionary to store mean temps for each pressure level
        mean_temps_by_level = {}

        # Iterate through the selected pressure levels
        for pressure_level in pressure_levels:
            print(f"Processing pressure level: {pressure_level} hPa")
            pres_index = np.where(pres == pressure_level)[0][0]  # Find the index for the given pressure level

            mean_temps = []

            # Iterate through the filtered time steps
            for t in time_indices:
                # Extract temp data for the current time step, pressure level, and latitude range
                temp_slice = temps[t, pres_index, lat_indices, :]  # Lazy slice
                mean_temp = np.mean(temp_slice)
                mean_temps.append(mean_temp)

            # Store the results for this pressure level
            mean_temps_by_level[pressure_level] = mean_temps

        # Select only the filtered times
        selected_times = [times_converted[t] for t in time_indices]

        print("Time conversion successful.")
        
        # Close the dataset
        dataset.close()
        
        return selected_times, mean_temps_by_level
    except Exception as e:
        print(f"An error occurred: {e}")
        return [], {}

def plot_mean_temp(times, mean_temps_by_level, lat_min, lat_max):
    if len(times) == 0 or len(mean_temps_by_level) == 0:
        print("No data to plot.")
        return
    
    # Plot the mean temps over time for each pressure level
    plt.figure(figsize=(10, 6))

    for pressure_level, mean_temps in mean_temps_by_level.items():
        plt.plot(times, mean_temps, marker='o', label=f'{pressure_level} hPa')

    # Customize plot
    plt.title(f'Mean temps ({lat_min}-{lat_max}°N)')
    plt.xlabel('Time')
    plt.ylabel('Temperature (K)')
    plt.grid(True)
    plt.legend()
    plt.show()

# Define the file path, latitude range, and pressure levels
file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Nov_Jan_2025_T_U_V.nc"  # Replace with your file path
lat_min, lat_max = 60, 90
pressure_levels = [1, 10]  # List of desired pressure levels in hPa , 300, 500, 700, 850, 1000

# Specify the date range
start_date = "2024-11-01"
end_date = "2025-01-29"

# Calculate the mean temp for each day in the specified date range and for the specified pressure levels
times, mean_temps_by_level = calculate_mean_temp(file_path, lat_min, lat_max, pressure_levels, start_date, end_date)

# Plot the results
plot_mean_temp(times, mean_temps_by_level, lat_min, lat_max)