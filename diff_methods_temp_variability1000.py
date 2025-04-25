# Import necessary libraries
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

# Load your dataset
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

# Extract temperature time series data for each city at the 1000 mb pressure level
city_temp_data = {}

for city, coords in cities.items():
    temp_series = ds['t_celsius'].sel(
        latitude=coords['lat'], longitude=coords['lon'], method='nearest'
    ).sel(pressure_level=1000.0)
    city_temp_data[city] = temp_series.to_series()

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

# Plotting temperature time series for each city
for city, temp_series in city_temp_data.items():
    temp_series = temp_series.dropna()
    dates = temp_series.index
    labels = label_ssw_periods(dates, ssw_dates)
    
    # Filter dates from 2010 onwards
    mask = dates >= pd.to_datetime('2010-01-01')
    dates = dates[mask]
    temp_series = temp_series[mask]
    labels = [label for idx, label in enumerate(labels) if mask[idx]]
    
    # Create a DataFrame for plotting
    df = pd.DataFrame({'Temperature': temp_series.values}, index=dates)
    df['SSW_Label'] = labels
    
    plt.figure(figsize=(14, 6))
    for label, group in df.groupby('SSW_Label'):
        plt.plot(group.index, group['Temperature'], label=label)
    
    plt.title(f'Temperature at 1000 mb in {city} (2010 - {current_date.year})')
    plt.xlabel('Date')
    plt.ylabel('Temperature (°C)')
    plt.xlim(pd.to_datetime('2010-01-01'), current_date)
    plt.legend()
    plt.show()

# Statistical comparison between SSW and Non-SSW periods
for city, temp_series in city_temp_data.items():
    temp_series = temp_series.dropna()
    dates = temp_series.index
    labels = label_ssw_periods(dates, ssw_dates)
    
    # Filter dates from 2010 onwards
    mask = dates >= pd.to_datetime('2010-01-01')
    dates = dates[mask]
    temp_series = temp_series[mask]
    labels = [label for idx, label in enumerate(labels) if mask[idx]]
    
    df = pd.DataFrame({'Temperature': temp_series.values}, index=dates)
    df['SSW_Label'] = labels
    
    ssw_temps = df[df['SSW_Label'] == 'SSW']['Temperature']
    non_ssw_temps = df[df['SSW_Label'] == 'Non-SSW']['Temperature']
    
    print(f"\n{city}:")
    print(f"Mean Temperature during SSW: {ssw_temps.mean():.2f}°C")
    print(f"Mean Temperature during Non-SSW: {non_ssw_temps.mean():.2f}°C")
    
    # Perform t-test
    t_stat, p_value = ttest_ind(ssw_temps, non_ssw_temps, equal_var=False)
    print(f"T-test p-value: {p_value:.4f}")
