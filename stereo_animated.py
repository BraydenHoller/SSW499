import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D

# Load the two geopotential height datasets
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
    np.datetime64("2023-02-16", 'D')
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

# Define the latitude and longitude boundaries
lat_min, lat_max = 20, 90  # From 20 degrees latitude up to the North Pole
lon_min, lon_max = -180, 180  # Global longitude coverage

# Colors for SSW and Non-SSW events
colors_ssw = ['red'] * len(ssw_event_dates)  # All SSW events in red
colors_non_ssw = ['navy'] * len(non_ssw_event_dates)  # All Non-SSW events in navy

# Prepare the figure and axes
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# Add map features
ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

# Collect legend elements
legend_elements = [
    Line2D([0], [0], color='red', lw=2, label='SSW Events'),
    Line2D([0], [0], color='navy', lw=2, label='Non-SSW Events')
]

# Combine SSW and Non-SSW events into one list with their corresponding datasets and colors
all_events = []
for date, label in zip(ssw_event_dates, ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': 'red', 'ds': ds_ssw})

for date, label in zip(non_ssw_event_dates, non_ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': 'navy', 'ds': ds_non_ssw})

# Define the total number of days to animate
total_days = 14  # Animate for 14 days after event dates

# Create a list of days to animate
days_to_animate = np.arange(0, total_days + 1)  # From day 0 to day 14

def update(day_offset):
    ax.clear()
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

    # Title for the current frame
    current_day = day_offset
    ax.set_title(f'5400 m Thickness Contours - Day {current_day} After Event Dates')

    # Plot contours for all events at the current day offset
    for event in all_events:
        event_date = event['date'] + np.timedelta64(day_offset, 'D')
        ds = event['ds']
        color = event['color']

        # Select the valid_time using method='nearest'
        try:
            geopotential = ds['z'].sel(valid_time=event_date, method='nearest')
        except KeyError:
            print(f"No data available for event date: {event_date}")
            continue

        # Then, select pressure levels and slice latitude and longitude
        geopotential = geopotential.sel(
            pressure_level=[500.0, 1000.0],
            latitude=slice(lat_max, lat_min),
            longitude=slice(lon_min, lon_max)
        )

        # Check if data is available
        if geopotential.isnull().all():
            print(f"No data available for event date: {event_date}")
            continue

        # Convert geopotential to geopotential height
        geopotential_height = geopotential / 9.80665

        # Extract heights for 500 mb and 1000 mb
        try:
            hgt_500mb = geopotential_height.sel(pressure_level=500.0)
            hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)
        except KeyError:
            print(f"Pressure levels not found for event date: {event_date}")
            continue

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
            linewidths=1
        )

    # Add the legend
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)

# Create the animation
ani = animation.FuncAnimation(
    fig,
    update,
    frames=days_to_animate,
    interval=1000,
    blit=False
)

# To save the animation as a video file (optional), uncomment the following lines:
# from matplotlib.animation import FFMpegWriter

from matplotlib.animation import FFMpegWriter
import matplotlib
matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
ani.save('thickness_animation.mp4', writer=writer)

# Display the animation
plt.show()