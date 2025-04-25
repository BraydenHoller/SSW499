# # Import necessary libraries
# import xarray as xr
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import ttest_ind
# from scipy.interpolate import interp1d

# # Load your dataset
# ds = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Since_2010_Temps_1000_500.nc")

# # Adjust longitudes from 0-360 to -180 to 180 if necessary
# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(
#             longitude=(((ds.longitude + 180) % 360) - 180)
#         )
#         ds = ds.sortby('longitude')
#     return ds

# ds = adjust_longitude(ds)

# # Ensure temperature is in Celsius
# ds['t_celsius'] = ds['t'] - 273.15  # Convert from Kelvin to Celsius if necessary

# # Define the cities with their coordinates
# cities = {
#     'Washington D.C.': {'lat': 38.9072, 'lon': -77.0369},
#     'Beijing': {'lat': 39.9042, 'lon': 116.4074},
#     'Moscow': {'lat': 55.7558, 'lon': 37.6173}
# }

# # Define SSW event dates (from 2010 onwards)
# ssw_dates = pd.to_datetime([
#     "2010-02-10", "2010-03-24", "2013-01-07",
#     "2018-02-12", "2019-01-02", "2021-01-05", "2023-02-16"
# ])

# # Function to label dates as 'SSW' or 'Non-SSW' based on proximity to SSW events
# def label_ssw_periods(dates, ssw_dates, window=15):
#     labels = []
#     for date in dates:
#         if any(abs((date - ssw_date).days) <= window for ssw_date in ssw_dates):
#             labels.append('SSW')
#         else:
#             labels.append('Non-SSW')
#     return labels

# # Get the current date for plotting limits
# current_date = pd.to_datetime('today')

# # Extract the temperature data at 1000 mb
# temperature = ds['t_celsius'].sel(pressure_level=1000.0)

# # Define the time range from 2010 onwards
# start_date = pd.to_datetime('2010-01-01')
# end_date = current_date

# # Reindex the temperature data to have a consistent daily time series
# temperature = temperature.sel(valid_time=slice(start_date, end_date)).resample(valid_time='1D').nearest()

# # Prepare a dictionary to store the 0°C latitude data for each city
# city_zero_deg_latitudes = {}

# for city, coords in cities.items():
#     print(f"Processing {city}...")
#     lon = coords['lon']
    
#     # Adjust longitude to the dataset's longitude coordinates
#     lon = ((lon + 180) % 360) - 180  # Ensure longitude is between -180 and 180
    
#     # Find the nearest longitude index in the dataset
#     lon_idx = ds.longitude.sel(longitude=lon, method='nearest')
    
#     # Extract temperature profiles along the city's longitude over time
#     temp_profile = temperature.sel(longitude=lon_idx)
    
#     # Initialize a list to store the latitude where temperature crosses 0°C
#     zero_deg_latitudes = []
#     dates = temp_profile.valid_time.values
#     for i in range(len(dates)):
#         temp_slice = temp_profile.isel(valid_time=i)
#         temp_values = temp_slice.values
#         latitudes = temp_slice.latitude.values
        
#         # Check if temperatures cross 0°C
#         if np.any((temp_values[:-1] * temp_values[1:]) <= 0):
#             # Interpolate to find the latitude where temperature is 0°C
#             f_interp = interp1d(temp_values, latitudes, kind='linear', bounds_error=False, fill_value=np.nan)
#             zero_deg_lat = f_interp(0)
#             zero_deg_latitudes.append(zero_deg_lat)
#         else:
#             # If the 0°C isotherm doesn't cross the longitude, append NaN
#             zero_deg_latitudes.append(np.nan)
    
#     # Store the results in the dictionary
#     city_zero_deg_latitudes[city] = pd.Series(zero_deg_latitudes, index=dates)

# # Create a DataFrame to hold the results for all cities
# df_zero_deg = pd.DataFrame(city_zero_deg_latitudes)

# # Label each date as 'SSW' or 'Non-SSW'
# dates = pd.to_datetime(df_zero_deg.index)
# labels = label_ssw_periods(dates, ssw_dates)
# df_zero_deg['SSW_Label'] = labels

# # Statistical analysis and plotting
# for city in cities.keys():
#     # Extract the latitudes
#     latitudes = df_zero_deg[city]
    
#     # Separate SSW and Non-SSW periods
#     ssw_mask = df_zero_deg['SSW_Label'] == 'SSW'
#     ssw_latitudes = latitudes[ssw_mask]
#     non_ssw_latitudes = latitudes[~ssw_mask]
    
#     # Convert to numeric data types, coerce errors to NaN
#     ssw_latitudes = pd.to_numeric(ssw_latitudes, errors='coerce')
#     non_ssw_latitudes = pd.to_numeric(non_ssw_latitudes, errors='coerce')
    
#     # Remove NaN values
#     ssw_latitudes = ssw_latitudes.dropna()
#     non_ssw_latitudes = non_ssw_latitudes.dropna()
    
#     # Convert to numpy arrays with float64 data type
#     ssw_latitudes = ssw_latitudes.values.astype(np.float64)
#     non_ssw_latitudes = non_ssw_latitudes.values.astype(np.float64)
    
#     # Check for empty arrays
#     if len(ssw_latitudes) == 0 or len(non_ssw_latitudes) == 0:
#         print(f"\n{city}:")
#         print(f"Not enough data to perform t-test.")
#         continue
    
#     # Calculate means
#     ssw_mean_lat = ssw_latitudes.mean()
#     non_ssw_mean_lat = non_ssw_latitudes.mean()
    
#     # Perform t-test
#     t_stat, p_value = ttest_ind(ssw_latitudes, non_ssw_latitudes, equal_var=False)
    
#     print(f"\n{city}:")
#     print(f"Mean Latitude of 0°C Isotherm during SSW: {ssw_mean_lat:.2f}°N")
#     print(f"Mean Latitude of 0°C Isotherm during Non-SSW: {non_ssw_mean_lat:.2f}°N")
#     print(f"T-test p-value: {p_value:.4f}")


# 500mb, -25C code block
# Import necessary libraries
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from scipy.interpolate import interp1d

# Load your dataset (ensure it includes 500 mb temperature data)
ds = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Since_2010_Temps_1000_500.nc")

# Adjust longitudes from 0-360 to -180 to 180 if necessary
def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(
            longitude=(((ds.longitude + 180) % 360) - 180)
        )
        ds = ds.sortby('longitude')
    return ds

ds = adjust_longitude(ds)

# Ensure temperature is in Celsius
ds['t_celsius'] = ds['t'] - 273.15  # Convert from Kelvin to Celsius if necessary

# Define the cities with their coordinates
cities = {
    'Washington D.C.': {'lat': 38.9072, 'lon': -77.0369},
    'Beijing': {'lat': 39.9042, 'lon': 116.4074},
    'Moscow': {'lat': 55.7558, 'lon': 37.6173}
}

# Define SSW event dates (from 2010 onwards)
ssw_dates = pd.to_datetime([
    "2010-02-10", "2010-03-24", "2013-01-07",
    "2018-02-12", "2019-01-02", "2021-01-05", "2023-02-16"
])

# Function to label dates as 'SSW' or 'Non-SSW' based on proximity to SSW events
def label_ssw_periods(dates, ssw_dates, window=15):
    labels = []
    for date in dates:
        if any(abs((date - ssw_date).days) <= window for ssw_date in ssw_dates):
            labels.append('SSW')
        else:
            labels.append('Non-SSW')
    return labels

# Get the current date for plotting limits
current_date = pd.to_datetime('today')

# **Extract the temperature data at 500 mb**
temperature = ds['t_celsius'].sel(pressure_level=500.0)

# Define the time range from 2010 onwards
start_date = pd.to_datetime('2010-01-01')
end_date = current_date

# Reindex the temperature data to have a consistent daily time series
temperature = temperature.sel(valid_time=slice(start_date, end_date)).resample(valid_time='1D').nearest()

# Prepare a dictionary to store the **-25°C** latitude data for each city
city_minus25_deg_latitudes = {}

for city, coords in cities.items():
    print(f"Processing {city}...")
    lon = coords['lon']
    
    # Adjust longitude to the dataset's longitude coordinates
    lon = ((lon + 180) % 360) - 180  # Ensure longitude is between -180 and 180
    
    # Find the nearest longitude index in the dataset
    lon_idx = ds.longitude.sel(longitude=lon, method='nearest')
    
    # Extract temperature profiles along the city's longitude over time
    temp_profile = temperature.sel(longitude=lon_idx)
    
    # Initialize a list to store the latitude where temperature crosses **-25°C**
    minus25_deg_latitudes = []
    dates = temp_profile.valid_time.values
    for i in range(len(dates)):
        temp_slice = temp_profile.isel(valid_time=i)
        temp_values = temp_slice.values
        latitudes = temp_slice.latitude.values
        
        # **Shift temperatures to find zero crossing at -25°C**
        temp_values_shifted = temp_values - (-25.0)
        
        # Check if temperatures cross -25°C
        if np.any((temp_values_shifted[:-1] * temp_values_shifted[1:]) <= 0):
            # Interpolate to find the latitude where temperature is -25°C
            try:
                f_interp = interp1d(temp_values_shifted, latitudes, kind='linear', bounds_error=False, fill_value=np.nan)
                minus25_deg_lat = f_interp(0)
                minus25_deg_latitudes.append(minus25_deg_lat)
            except ValueError:
                # Handle cases where interpolation fails
                minus25_deg_latitudes.append(np.nan)
        else:
            # If the -25°C isotherm doesn't cross the longitude, append NaN
            minus25_deg_latitudes.append(np.nan)
    
    # Store the results in the dictionary
    city_minus25_deg_latitudes[city] = pd.Series(minus25_deg_latitudes, index=dates)

# Create a DataFrame to hold the results for all cities
df_minus25_deg = pd.DataFrame(city_minus25_deg_latitudes)

# Label each date as 'SSW' or 'Non-SSW'
dates = pd.to_datetime(df_minus25_deg.index)
labels = label_ssw_periods(dates, ssw_dates)
df_minus25_deg['SSW_Label'] = labels

# Statistical analysis and plotting
for city in cities.keys():
    # Extract the latitudes
    latitudes = df_minus25_deg[city]
    
    # Separate SSW and Non-SSW periods
    ssw_mask = df_minus25_deg['SSW_Label'] == 'SSW'
    ssw_latitudes = latitudes[ssw_mask]
    non_ssw_latitudes = latitudes[~ssw_mask]
    
    # Convert to numeric data types, coerce errors to NaN
    ssw_latitudes = pd.to_numeric(ssw_latitudes, errors='coerce')
    non_ssw_latitudes = pd.to_numeric(non_ssw_latitudes, errors='coerce')
    
    # Remove NaN values
    ssw_latitudes = ssw_latitudes.dropna()
    non_ssw_latitudes = non_ssw_latitudes.dropna()
    
    # Convert to numpy arrays with float64 data type
    ssw_latitudes = ssw_latitudes.values.astype(np.float64)
    non_ssw_latitudes = non_ssw_latitudes.values.astype(np.float64)
    
    # Check for empty arrays
    if len(ssw_latitudes) == 0 or len(non_ssw_latitudes) == 0:
        print(f"\n{city}:")
        print(f"Not enough data to perform t-test.")
        continue
    
    # Calculate means
    ssw_mean_lat = ssw_latitudes.mean()
    non_ssw_mean_lat = non_ssw_latitudes.mean()
    
    # Perform t-test
    t_stat, p_value = ttest_ind(ssw_latitudes, non_ssw_latitudes, equal_var=False)
    
    print(f"\n{city}:")
    print(f"Mean Latitude of -25°C Isotherm during SSW: {ssw_mean_lat:.2f}°N")
    print(f"Mean Latitude of -25°C Isotherm during Non-SSW: {non_ssw_mean_lat:.2f}°N")
    print(f"T-test p-value: {p_value:.4f}")