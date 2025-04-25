import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib

matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"

# Load the dataset
ds = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Snow_2223.nc")
snowfall = ds['sf']

start_date = "2022-07-01"
end_date = "2023-07-01"
start_datetime = np.datetime64(start_date)
end_datetime = np.datetime64(end_date)

# Slice the dataset to the desired time range
snowfall_timeline = snowfall.sel(valid_time=slice(start_datetime, end_datetime))

# Define US region boundaries
lat_min, lat_max = 24, 50
lon_min, lon_max = 235, 294

# Subset the data to the US region
snowfall_us = snowfall_timeline.sel(
    latitude=slice(lat_max, lat_min),
    longitude=slice(lon_min, lon_max)
)

# Create arrays of full degree values
lat_full_degrees = np.arange(np.ceil(lat_min), np.floor(lat_max) + 1, 1)
lon_full_degrees = np.arange(np.ceil(lon_min), np.floor(lon_max) + 1, 1)

# Select data at these latitude and longitude values
snowfall_thinned = snowfall_us.sel(
    latitude=lat_full_degrees,
    longitude=lon_full_degrees,
    method='nearest'
)

snowfall_thinned_true = snowfall_thinned / 0.1

# Create the animation
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

def update(frame):
    ax.clear()
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='black')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle=':')
    
    data = snowfall_thinned_true.isel(valid_time=frame)
    data.plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='Blues',
        add_colorbar=False,
        vmin=0,
        vmax=snowfall_thinned_true.max()
    )
    ax.set_title(f'Snowfall on {str(data.valid_time.values)[:16]}')

from matplotlib.animation import FFMpegWriter

writer = FFMpegWriter(fps=5, metadata=dict(artist='Brayden Holler'), bitrate=1800)

# Create the animation (using your existing code)
ani = animation.FuncAnimation(
    fig, update, frames=len(snowfall_thinned.valid_time), interval=200, blit=False
)

ani.save('us_snowfall_animation0722-0723.mp4', writer=writer)


# Display the animation in a Jupyter notebook
from IPython.display import Video

Video('us_snowfall_animation0722-0723.mp4', embed=True)