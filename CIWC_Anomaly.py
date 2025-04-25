import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.stats import pearsonr
from matplotlib.lines import Line2D

##############################################
# 1. Load the Datasets                       #
##############################################

# Replace these file paths with the locations of your datasets.
precip_path = 'path/to/precipitation_data.nc'
cloud_path = 'path/to/cloud_data.nc'

# Load precipitation and cloud water data.
# It is assumed that the datasets have a time coordinate named 'time'
# and spatial coordinates 'latitude' and 'longitude'.
ds_precip = xr.open_dataset(precip_path)
ds_cloud = xr.open_dataset(cloud_path)

# For the cloud dataset, we assume there are two variables:
# 'ciwc' for Specific Cloud Ice Water Content and 'clwc' for Specific Cloud Liquid Water Content.
# Adjust the variable names if they differ in your data.

##############################################
# 2. Compute Anomalies                       #
##############################################

# Compute the overall time mean for each variable as a simple baseline.
precip_mean = ds_precip['precipitation'].mean(dim='time')
ciwc_mean = ds_cloud['ciwc'].mean(dim='time')
clwc_mean = ds_cloud['clwc'].mean(dim='time')

# Compute anomaly fields (i.e. deviation from the overall mean)
precip_anom = ds_precip['precipitation'] - precip_mean
ciwc_anom = ds_cloud['ciwc'] - ciwc_mean
clwc_anom = ds_cloud['clwc'] - clwc_mean

##############################################
# 3. Define SSW Event Dates and Time Window    #
##############################################

# Use your provided SSW event dates.
ssw_event_dates = [
    np.datetime64("2010-02-10", 'D'),
    np.datetime64("2010-03-24", 'D'),
    np.datetime64("2013-01-06", 'D'),
    np.datetime64("2018-02-12", 'D'),
    np.datetime64("2019-01-01", 'D'),
    np.datetime64("2021-01-05", 'D'),
    np.datetime64("2023-02-16", 'D')
]

# Define a time window around each event, e.g. 7 days before and after the event.
window_days = 7

##############################################
# 4. Build Composite Anomalies for SSW Events  #
##############################################

# Create lists to store the composites.
precip_composites = []
ciwc_composites = []
clwc_composites = []
time_coords = None  # to save the relative time coordinate

# Loop over each SSW event.
for event_date in ssw_event_dates:
    # Define the start and end of the analysis window.
    start_time = event_date - np.timedelta64(window_days, 'D')
    end_time = event_date + np.timedelta64(window_days, 'D')
    
    # Select data from the anomaly fields within this window.
    precip_event = precip_anom.sel(time=slice(start_time, end_time))
    ciwc_event = ciwc_anom.sel(time=slice(start_time, end_time))
    clwc_event = clwc_anom.sel(time=slice(start_time, end_time))
    
    # Create a relative time coordinate (days from event)
    rel_time = (precip_event['time'] - event_date) / np.timedelta64(1, 'D')
    precip_event = precip_event.assign_coords(rel_time=rel_time)
    ciwc_event = ciwc_event.assign_coords(rel_time=rel_time)
    clwc_event = clwc_event.assign_coords(rel_time=rel_time)
    
    # Append to our composite lists.
    precip_composites.append(precip_event)
    ciwc_composites.append(ciwc_event)
    clwc_composites.append(clwc_event)
    
    # Save the relative time axis (assumed to be identical across events)
    if time_coords is None:
        time_coords = rel_time

# Concatenate along a new 'event' dimension and average over events.
precip_comp = xr.concat(precip_composites, dim='event').mean(dim='event')
ciwc_comp = xr.concat(ciwc_composites, dim='event').mean(dim='event')
clwc_comp = xr.concat(clwc_composites, dim='event').mean(dim='event')

##############################################
# 5. Plot Composite Anomaly Maps             #
##############################################

# Use a Plate Carree projection for geographic plots.
proj = ccrs.PlateCarree()
fig, axs = plt.subplots(1, 3, subplot_kw={'projection': proj}, figsize=(18, 6))

# Helper function to add map features.
def setup_map(ax, title):
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.set_title(title)

# Plot composite precipitation anomaly (using day 0, the event day).
day_idx = np.where(time_coords == 0)[0][0]
im0 = precip_comp.isel(rel_time=day_idx).plot(ax=axs[0], transform=proj, add_colorbar=False)
setup_map(axs[0], "Precipitation Anomaly (Day 0)")

# Plot composite cloud ice water anomaly.
im1 = ciwc_comp.isel(rel_time=day_idx).plot(ax=axs[1], transform=proj, add_colorbar=False)
setup_map(axs[1], "CIWC Anomaly (Day 0)")

# Plot composite cloud liquid water anomaly.
im2 = clwc_comp.isel(rel_time=day_idx).plot(ax=axs[2], transform=proj, add_colorbar=False)
setup_map(axs[2], "CLWC Anomaly (Day 0)")

# Add a common colorbar.
fig.colorbar(im0, ax=axs, orientation='horizontal', fraction=0.05)
plt.tight_layout()
plt.show()

##############################################
# 6. Regional Time Series & Correlation Analysis
##############################################

# Define a region of interest (update these limits as needed).
region = {
    'lat_min': 40.0,
    'lat_max': 50.0,
    'lon_min': -100.0,
    'lon_max': -90.0
}

# Compute area-averaged time series for the full anomaly datasets over the region.
precip_region = precip_anom.sel(latitude=slice(region['lat_min'], region['lat_max']),
                                longitude=slice(region['lon_min'], region['lon_max'])).mean(dim=['latitude','longitude'])
ciwc_region = ciwc_anom.sel(latitude=slice(region['lat_min'], region['lat_max']),
                            longitude=slice(region['lon_min'], region['lon_max'])).mean(dim=['latitude','longitude'])
clwc_region = clwc_anom.sel(latitude=slice(region['lat_min'], region['lat_max']),
                            longitude=slice(region['lon_min'], region['lon_max'])).mean(dim=['latitude','longitude'])

# Plot the area-averaged time series.
plt.figure(figsize=(10, 6))
plt.plot(precip_region['time'], precip_region, label='Precipitation Anomaly')
plt.plot(ciwc_region['time'], ciwc_region, label='CIWC Anomaly')
plt.plot(clwc_region['time'], clwc_region, label='CLWC Anomaly')
plt.xlabel('Time')
plt.ylabel('Anomaly')
plt.legend()
plt.title('Area-Averaged Anomaly Time Series')
plt.show()

# Compute Pearson correlations between precipitation and the cloud variables.
precip_vals = precip_region.values
ciwc_vals = ciwc_region.values
clwc_vals = clwc_region.values

# Remove any NaN values.
mask = np.isfinite(precip_vals) & np.isfinite(ciwc_vals) & np.isfinite(clwc_vals)
r_precip_ciwc, p_precip_ciwc = pearsonr(precip_vals[mask], ciwc_vals[mask])
r_precip_clwc, p_precip_clwc = pearsonr(precip_vals[mask], clwc_vals[mask])

print("Correlation between precipitation and CIWC anomalies: r = {:.3f}, p = {:.3e}".format(r_precip_ciwc, p_precip_ciwc))
print("Correlation between precipitation and CLWC anomalies: r = {:.3f}, p = {:.3e}".format(r_precip_clwc, p_precip_clwc))
