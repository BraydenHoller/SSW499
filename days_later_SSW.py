import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D

# Load the new geopotential height dataset
# Replace the path with the path to your new dataset
ds = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSWs_stereo.nc")

# Adjust longitudes from 0-360 to -180 to 180 if necessary
if ds.longitude.max() > 180:
    ds = ds.assign_coords(
        longitude=(((ds.longitude + 180) % 360) - 180)
    )
    ds = ds.sortby('longitude')

# Define outset dates of events with explicit units
event_dates = [
    np.datetime64("2010-02-10", 'D'),
    np.datetime64("2010-03-24", 'D'),
    np.datetime64("2013-01-07", 'D'),
    np.datetime64("2018-02-12", 'D'),
    np.datetime64("2019-01-02", 'D'),
    np.datetime64("2021-01-05", 'D'),
    np.datetime64("2022-02-16", 'D')
]

# Event labels for legend
event_labels = [
    "2010-02-10",
    "2010-03-24",
    "2013-01-07",
    "2018-02-12",
    "2019-01-02",
    "2021-01-05",
    "2023-02-16"
]

# Define days since event outset
days_later = np.timedelta64(3, 'D')

# Define the latitude and longitude boundaries
lat_min, lat_max = 20, 90  # From 20 degrees latitude up to the North Pole
lon_min, lon_max = -180, 180  # Global longitude coverage

# List of colors for each event
colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'cyan']

# Prepare the figure and axes
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# Add map features
ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

def plot(ax, date, color, label):
    # Calculate the event date
    event_date = date + days_later

    # Select the valid_time using method='nearest'
    geopotential = ds['z'].sel(valid_time=event_date, method='nearest')

    # Then, select pressure levels and slice latitude and longitude
    geopotential = geopotential.sel(
        pressure_level=[500.0, 1000.0],
        latitude=slice(lat_max, lat_min),
        longitude=slice(lon_min, lon_max)
    )

    # Check if data is available
    if geopotential.isnull().all():
        print(f"No data available for event date: {event_date}")
        return

    # Convert geopotential to geopotential height
    geopotential_height = geopotential / 9.80665

    # Extract heights for 500 mb and 1000 mb
    hgt_500mb = geopotential_height.sel(pressure_level=500.0)
    hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)

    # Calculate the thickness
    thickness = hgt_500mb - hgt_1000mb

    # Plot the 5400 m thickness contour
    cs = ax.contour(
        thickness.longitude,
        thickness.latitude,
        thickness,
        levels=[5400],
        colors=color,
        transform=ccrs.PlateCarree(),
        linewidths=2
    )

    # Check if contour was created
    if len(cs.allsegs[0]) == 0:
        print(f"No contour at 5400 m for event date: {event_date}")
        return

    # For legend purposes, create a custom Line2D object
    line = Line2D([0], [0], color=color, lw=2, label=label)
    return line

# Collect legend elements
legend_elements = []

# Loop over events, colors, and labels, plotting onto the same axes
for date, color, label in zip(event_dates, colors, event_labels):
    line = plot(ax, date, color, label)
    if line:
        legend_elements.append(line)

# Add a title
ax.set_title(f'5400 m Thickness Contours {days_later.astype(int)} Days After SSW Events')

# Add the legend
ax.legend(handles=legend_elements, loc='lower left')

# Display the plot
plt.show()