import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D

# Load the two geopotential height datasets
# Replace the paths with the paths to your datasets
ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSWs_stereo.nc")
ds_non_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\non_SSWs_stereo.nc")

# Adjust longitudes from 0-360 to -180 to 180 if necessary for both datasets
def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(
            longitude=(((ds.longitude + 180) % 360) - 180)
        )
        ds = ds.sortby('longitude')
    return ds

ds_ssw = adjust_longitude(ds_ssw)
ds_non_ssw = adjust_longitude(ds_non_ssw)

# Define event dates and labels for both datasets
# SSW events
ssw_event_dates = [
    np.datetime64("2010-02-10", 'D'),
    np.datetime64("2010-03-24", 'D'),
    np.datetime64("2013-01-07", 'D'),
    np.datetime64("2018-02-12", 'D'),
    np.datetime64("2019-01-02", 'D'),
    np.datetime64("2021-01-05", 'D'),
    np.datetime64("2022-02-16", 'D')
]
ssw_event_labels = [
    "SSW 2010-02-10",
    "SSW 2010-03-24",
    "SSW 2013-01-07",
    "SSW 2018-02-12",
    "SSW 2019-01-02",
    "SSW 2021-01-05",
    "SSW 2023-02-16"
]

# Non-SSW events
non_ssw_event_dates = [
    np.datetime64("2012-01-05", 'D'),
    np.datetime64("2014-01-05", 'D'),
    np.datetime64("2015-01-05", 'D'),
    np.datetime64("2016-01-05", 'D'),
    np.datetime64("2017-01-05", 'D'),
    np.datetime64("2020-01-05", 'D'),
    np.datetime64("2022-01-05", 'D')
]
non_ssw_event_labels = [
    "Non-SSW 2012-01-05",
    "Non-SSW 2014-01-05",
    "Non-SSW 2015-01-05",
    "Non-SSW 2016-01-05",
    "Non-SSW 2017-01-05",
    "Non-SSW 2020-01-05",
    "Non-SSW 2022-01-05"
]

# Define days since event outset
days_later_ssw = np.timedelta64(7, 'D')      # For SSW events
days_later_non_ssw = np.timedelta64(0, 'D')  # For Non-SSW events

# Define the latitude and longitude boundaries
lat_min, lat_max = 20, 90  # From 20 degrees latitude up to the North Pole
lon_min, lon_max = -180, 180  # Global longitude coverage

# Colors for SSW and Non-SSW events
colors_ssw = ['red', 'red', 'red', 'red', 'red', 'red', 'red']
colors_non_ssw = ['navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy']

# Prepare the figure and axes
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# Add map features
ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

def plot(ax, date, days_later, color, label, ds):
    # Calculate the event date
    event_date = date + days_later

    # Select the valid_time using method='nearest'
    try:
        geopotential = ds['z'].sel(valid_time=event_date, method='nearest')
    except KeyError:
        print(f"No data available for event date: {event_date}")
        return

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

# Plot SSW events
for date, color, label in zip(ssw_event_dates, colors_ssw, ssw_event_labels):
    line = plot(ax, date, days_later_ssw, color, label, ds_ssw)
    if line:
        legend_elements.append(line)

# Plot Non-SSW events
for date, color, label in zip(non_ssw_event_dates, colors_non_ssw, non_ssw_event_labels):
    line = plot(ax, date, days_later_non_ssw, color, label, ds_non_ssw)
    if line:
        legend_elements.append(line)

# Add a title
ax.set_title(f'5400 m Thickness Contours {days_later_ssw.astype(int)} Days After SSW Events')

# Add the legend
ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)

# Adjust layout to make room for the legend
plt.tight_layout()

# Display the plot
plt.show()