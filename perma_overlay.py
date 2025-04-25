import rasterio
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D

#############################################
# 1. Load and subset the permafrost raster  #
#############################################

# Define the path to your .byte permafrost file
byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\llipa.byte'

# Open the permafrost file with Rasterio and subset to latitudes above 25°N.
with rasterio.open(byte_file_path) as perma_ds:
    # Read the first band as a NumPy array
    perma_data = perma_ds.read(1)
    
    # Determine the pixel size in the y-direction (assumes north-up image)
    pixel_size_y = abs(perma_ds.transform.e)
    nrows = perma_data.shape[0]
    rows = np.arange(nrows)
    
    # Compute each row's center latitude.
    center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    
    # Find indices where the center latitude is above 25°N
    valid_rows = np.where(center_lats > 25)[0]
    if valid_rows.size == 0:
        raise ValueError("No rows found with latitude above 25°N.")
    
    # Since rows are ordered from north to south, take rows from the top down to the last valid row.
    last_valid_row = valid_rows[-1]
    perma_subset = perma_data[:last_valid_row + 1, :]
    
    # Compute the new geographic extent from the raster bounds.
    new_left = perma_ds.bounds.left
    new_right = perma_ds.bounds.right
    new_top = perma_ds.bounds.top
    # For row i, the bottom edge is: top - (i+1)*pixel_size_y.
    new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
    # This extent is in (left, right, bottom, top) order.
    perma_extent = (new_left, new_right, new_bottom, new_top)

#############################################
# 2. Load geopotential height datasets      #
#############################################

# Load the two geopotential height datasets
ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSWs_stereo.nc")
ds_non_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\non_SSWs_stereo.nc")

# Function to adjust longitudes from 0-360 to -180 to 180 if needed
def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(
            longitude=(((ds.longitude + 180) % 360) - 180)
        )
        ds = ds.sortby('longitude')
    return ds

ds_ssw = adjust_longitude(ds_ssw)
ds_non_ssw = adjust_longitude(ds_non_ssw)

#############################################
# 3. Define event dates, labels, and colors   #
#############################################

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

# Colors for SSW and Non-SSW events
colors_ssw = ['red'] * len(ssw_event_dates)
colors_non_ssw = ['navy'] * len(non_ssw_event_dates)

# Combine events into one list (each event dict includes its date, label, color, and dataset)
all_events = []
for date, label in zip(ssw_event_dates, ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': 'red', 'ds': ds_ssw})
for date, label in zip(non_ssw_event_dates, non_ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': 'navy', 'ds': ds_non_ssw})

#############################################
# 4. Prepare the figure and animation axis    #
#############################################

# We'll use a Cartopy axis with North Polar Stereo projection.
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
# Set the extent using the permafrost raster extent.
ax.set_extent(perma_extent, crs=ccrs.PlateCarree())

# Add basic map features
ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

# Prepare legend elements for the events
legend_elements = [
    Line2D([0], [0], color='red', lw=2, label='SSW Events'),
    Line2D([0], [0], color='navy', lw=2, label='Non-SSW Events')
]

#############################################
# 5. Define animation parameters and update   #
#############################################

# Define the total number of days to animate after each event date.
total_days = 14
days_to_animate = np.arange(0, total_days + 1)  # Days 0 to 14

def update(day_offset):
    # Clear the axis and re-set the projection and extent.
    ax.clear()
    ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
    
    # Re-add map features
    ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

    # Plot the permafrost background.
    # We use transform=ccrs.PlateCarree() since the image's extent is in geographic coordinates.
    ax.imshow(
        perma_subset,
        transform=ccrs.PlateCarree(),
        extent=perma_extent,
        origin='upper',
        cmap='viridis',
        alpha=0.5  # Adjust transparency as desired.
    )

    # Add a title that shows the current day offset.
    current_day = day_offset
    ax.set_title(f'5400 m Thickness Contours - Day {current_day} After Event Dates', fontsize=14)

    # Plot the thickness contour for each event at the current day offset.
    for event in all_events:
        event_date = event['date'] + np.timedelta64(day_offset, 'D')
        ds = event['ds']
        color = event['color']

        # Attempt to select the valid_time using method='nearest'
        try:
            geopotential = ds['z'].sel(valid_time=event_date, method='nearest')
        except KeyError:
            print(f"No data available for event date: {event_date}")
            continue

        # For consistency with the permafrost background, we select the same geographic domain.
        try:
            geopotential = geopotential.sel(
                pressure_level=[500.0, 1000.0],
                latitude=slice(new_top, new_bottom),  # from high (top) to low (bottom)
                longitude=slice(new_left, new_right)
            )
        except Exception as e:
            print(f"Error selecting data for event date {event_date}: {e}")
            continue

        if geopotential.isnull().all():
            print(f"No data available (all NaN) for event date: {event_date}")
            continue

        # Convert geopotential to geopotential height (divide by g = 9.80665 m/s²)
        geopotential_height = geopotential / 9.80665

        try:
            hgt_500mb = geopotential_height.sel(pressure_level=500.0)
            hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)
        except KeyError:
            print(f"Pressure levels not found for event date: {event_date}")
            continue

        # Compute the thickness (difference between 500 mb and 1000 mb heights)
        thickness = hgt_500mb - hgt_1000mb

        # Plot the 5400 m thickness contour.
        # The contour is plotted in PlateCarree coordinates.
        cs = ax.contour(
            thickness.longitude,
            thickness.latitude,
            thickness,
            levels=[5400],
            colors=color,
            transform=ccrs.PlateCarree(),
            linewidths=1.5
        )

    # Add the legend for SSW and Non-SSW events.
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)

#############################################
# 6. Create and save the animation           #
#############################################

ani = animation.FuncAnimation(
    fig,
    update,
    frames=days_to_animate,
    interval=1000,   # 1 second per frame; adjust as needed
    blit=False
)

# (Optional) Save the animation to a video file.
from matplotlib.animation import FFMpegWriter
import matplotlib
# Set the path to your ffmpeg executable
matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
ani.save('thickness_permafrost_animation.mp4', writer=writer)

plt.show()