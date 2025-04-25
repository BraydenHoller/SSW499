import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import pyplot as plt

# Set the path to ffmpeg
matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"

# Load the geopotential height dataset
ds = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\5400m_0z.nc")

# Adjust longitudes from 0-360 to -180 to 180
ds = ds.assign_coords(
    longitude=(((ds.longitude + 180) % 360) - 180)
)
ds = ds.sortby('longitude')

# Load the precipitation dataset
precip_ds = xr.open_dataset(r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\precip_540_0z.nc')

# Adjust longitudes in precipitation dataset
precip_ds = precip_ds.assign_coords(
    longitude=(((precip_ds.longitude + 180) % 360) - 180)
)
precip_ds = precip_ds.sortby('longitude')

# Define time range
start_date = "2022-09-01"
end_date = "2023-06-01"
start_datetime = np.datetime64(start_date)
end_datetime = np.datetime64(end_date)

# Get list of valid times within the time range
valid_times1 = ds['valid_time'].sel(valid_time=slice(start_datetime, end_datetime))
valid_times2 = precip_ds['valid_time'].sel(valid_time=slice(start_datetime, end_datetime))

# Prepare the figure and axis
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Create an axes on the right side for the colorbar
divider = make_axes_locatable(ax)
from matplotlib.axes import Axes
cax = divider.append_axes("right", size="5%", pad=0.05, axes_class=Axes)

# Initialize global variables
precip_plot = None
cbar = None

# Define US region boundaries
lat_min, lat_max = 24, 50
lon_min_original, lon_max_original = 235, 294  # Original longitudes in degrees East

# Convert longitude bounds to -180 to 180 degrees
lon_min = lon_min_original - 360  # Converts to -125 degrees
lon_max = lon_max_original - 360  # Converts to -66 degrees

# Subset precipitation data to US region
precipitation_us = precip_ds['lsp'].sel(
    latitude=slice(lat_max, lat_min),
    longitude=slice(lon_min, lon_max)
)

# Align precipitation times with thickness data
precipitation_us = precipitation_us.sel(valid_time=valid_times2)

def update(frame):
    global precip_plot, cbar

    ax.clear()
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='black')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle=':')

    # Get the current time
    current_time = valid_times2[frame].values

    # Select data for the current time step and pressure levels
    geopotential = ds['z'].sel(
        valid_time=current_time,
        pressure_level=[500.0, 1000.0],
        latitude=slice(lat_max, lat_min),
        longitude=slice(lon_min, lon_max)
    )

    # Check if data is available
    if geopotential.longitude.size == 0 or geopotential.latitude.size == 0:
        print(f"No data available for time {current_time}")
        return

    # Convert geopotential to geopotential height
    geopotential_height = geopotential / 9.80665

    # Extract heights for 500mb and 1000mb
    hgt_500mb = geopotential_height.sel(pressure_level=500.0)
    hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)

    # Calculate the thickness
    thickness = hgt_500mb - hgt_1000mb

    # Check if thickness data is valid
    if thickness.longitude.size == 0 or thickness.latitude.size == 0:
        print(f"No thickness data available for time {current_time}")
        return

    # Plot the precipitation data
    precip_data = precipitation_us.sel(valid_time=current_time)

    # Handle potential missing data
    if precip_data.longitude.size == 0 or precip_data.latitude.size == 0:
        print(f"No precipitation data available for time {current_time}")
    else:
        # Remove the previous precipitation plot if it exists
        if precip_plot is not None:
            precip_plot.remove()

        # Create the pcolormesh
        precip_plot = ax.pcolormesh(
            precip_data.longitude,
            precip_data.latitude,
            precip_data,
            cmap='Blues',
            shading='auto',
            transform=ccrs.PlateCarree()
        )

        # Update the colorbar
        if cbar is None:
            cbar = fig.colorbar(precip_plot, cax=cax, orientation='vertical')
            cbar.set_label('Precipitation (mm)')  # Replace with correct units
        else:
            cbar.update_normal(precip_plot)

    # Plot the 5400 m thickness contour
    cs = ax.contour(
        thickness.longitude,
        thickness.latitude,
        thickness,
        levels=[5400],
        colors='red',
        transform=ccrs.PlateCarree()
    )

    # Check if contour was created successfully
    if len(cs.allsegs[0]) == 0:
        print(f"No contour found for time {current_time}")
        return

    # Add labels to the contour
    ax.clabel(cs, fmt='%d', inline=True, fontsize=10)

    # Add a title with the current time
    ax.set_title(f'5400 m Thickness and Precipitation over CONUS\nTime: {np.datetime_as_string(current_time, unit="h")}')


from matplotlib.animation import FFMpegWriter

writer = FFMpegWriter(fps=5, metadata=dict(artist='Brayden Holler'), bitrate=1800)

# Create the animation
ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(valid_times1),
    interval=200,
    blit=False
)

# Save the animation
ani.save('5400m_thickness_precipitation_animation2.mp4', writer=writer)

#plot a set # days after events start