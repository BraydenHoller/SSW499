# import rasterio
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature

# #############################################
# # 1. Load and subset the permafrost raster  #
# #############################################

# # Define the path to your .byte permafrost file
# byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\PermaFrost\llipa.byte'

# with rasterio.open(byte_file_path) as perma_ds:
#     # Read the first band as a NumPy array.
#     perma_data = perma_ds.read(1)
    
#     # Determine pixel size in the y-direction (assumes north-up image)
#     pixel_size_y = abs(perma_ds.transform.e)
#     nrows = perma_data.shape[0]
#     rows = np.arange(nrows)
    
#     # Compute each row's center latitude.
#     center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    
#     # Find indices where the center latitude is above 25°N.
#     valid_rows = np.where(center_lats > 25)[0]
#     if valid_rows.size == 0:
#         raise ValueError("No rows found with latitude above 25°N.")
    
#     # Rows are ordered from north to south; select from the top to the last valid row.
#     last_valid_row = valid_rows[-1]
#     perma_subset = perma_data[:last_valid_row + 1, :]
    
#     # Compute the new geographic extent from the raster bounds.
#     new_left = perma_ds.bounds.left
#     new_right = perma_ds.bounds.right
#     new_top = perma_ds.bounds.top
#     new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
#     perma_extent = (new_left, new_right, new_bottom, new_top)

# #############################################
# # 2. Load geopotential height datasets      #
# #############################################

# # Load the two geopotential height datasets
# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSWs_stereo.nc")
# ds_non_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\non_SSWs_stereo.nc")

# # Adjust longitudes from 0–360 to –180–180 if needed.
# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(
#             longitude=(((ds.longitude + 180) % 360) - 180)
#         )
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)
# ds_non_ssw = adjust_longitude(ds_non_ssw)

# #############################################
# # 3. Define event dates (for both datasets)  #
# #############################################

# # SSW event dates
# ssw_event_dates = [
#     np.datetime64("2010-02-10", 'D'),
#     np.datetime64("2010-03-24", 'D'),
#     np.datetime64("2013-01-07", 'D'),
#     np.datetime64("2018-02-12", 'D'),
#     np.datetime64("2019-01-02", 'D'),
#     np.datetime64("2021-01-05", 'D'),
#     np.datetime64("2023-02-16", 'D')
# ]

# # Non-SSW event dates
# non_ssw_event_dates = [
#     np.datetime64("2012-01-05", 'D'),
#     np.datetime64("2014-01-05", 'D'),
#     np.datetime64("2015-01-05", 'D'),
#     np.datetime64("2016-01-05", 'D'),
#     np.datetime64("2017-01-05", 'D'),
#     np.datetime64("2020-01-05", 'D'),
#     np.datetime64("2022-01-05", 'D')
# ]

# # Combine events into one list.
# # (Each event now includes its date and corresponding dataset.)
# all_events = []
# for date in ssw_event_dates:
#     all_events.append({'date': date, 'ds': ds_ssw})
# for date in non_ssw_event_dates:
#     all_events.append({'date': date, 'ds': ds_non_ssw})

# #############################################
# # 4. Prepare the figure and base map         #
# #############################################

# fig = plt.figure(figsize=(12, 12))
# ax = plt.axes(projection=ccrs.NorthPolarStereo())
# ax.set_extent(perma_extent, crs=ccrs.PlateCarree())

# # Add basic map features.
# ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
# gl.top_labels = False
# gl.right_labels = False

# # Plot the permafrost background image.
# ax.imshow(
#     perma_subset,
#     transform=ccrs.PlateCarree(),
#     extent=perma_extent,
#     origin='upper',
#     cmap='viridis',
#     alpha=0.5
# )

# #############################################
# # 5. Define animation update function       #
# #############################################

# # Define the total number of days (0 to 14 after the event)
# total_days = 14
# days_to_animate = np.arange(0, total_days + 1)

# def update(day_offset):
#     # Clear and reconfigure the axis.
#     ax.clear()
#     ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
#     ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
#     ax.add_feature(cfeature.COASTLINE)
#     ax.add_feature(cfeature.BORDERS, linestyle=':')
#     gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False
    
#     # Re-plot the permafrost background.
#     ax.imshow(
#         perma_subset,
#         transform=ccrs.PlateCarree(),
#         extent=perma_extent,
#         origin='upper',
#         cmap='viridis',
#         alpha=0.5
#     )
    
#     # List to hold thickness fields for each event at the given day offset.
#     thickness_fields = []
    
#     # Loop over all events and compute the thickness field for the current day.
#     for event in all_events:
#         event_date = event['date'] + np.timedelta64(day_offset, 'D')
#         ds = event['ds']
        
#         try:
#             # Select the nearest valid_time.
#             geopotential = ds['z'].sel(valid_time=event_date, method='nearest')
#         except Exception as e:
#             print(f"No data for {event_date}: {e}")
#             continue
        
#         try:
#             geopotential = geopotential.sel(
#                 pressure_level=[500.0, 1000.0],
#                 latitude=slice(new_top, new_bottom),  # high-to-low
#                 longitude=slice(new_left, new_right)
#             )
#         except Exception as e:
#             print(f"Error selecting domain for {event_date}: {e}")
#             continue
        
#         if geopotential.isnull().all():
#             print(f"All values are NaN for {event_date}")
#             continue
        
#         # Convert geopotential to geopotential height (in m).
#         geopotential_height = geopotential / 9.80665
        
#         try:
#             hgt_500mb = geopotential_height.sel(pressure_level=500.0)
#             hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)
#         except Exception as e:
#             print(f"Pressure levels issue for {event_date}: {e}")
#             continue
        
#         # Compute thickness as the difference between 500 mb and 1000 mb.
#         thickness = hgt_500mb - hgt_1000mb
#         thickness_fields.append(thickness)
    
#     # Only proceed if we have at least one valid thickness field.
#     if thickness_fields:
#         # Concatenate along a new "event" dimension and compute the average.
#         thickness_all = xr.concat(thickness_fields, dim='event')
#         avg_thickness = thickness_all.mean(dim='event')
        
#         # Plot the 5400 m thickness contour from the averaged thickness field.
#         cs = ax.contour(
#             avg_thickness.longitude,
#             avg_thickness.latitude,
#             avg_thickness,
#             levels=[5400],
#             colors='black',
#             transform=ccrs.PlateCarree(),
#             linewidths=2
#         )
    
#     ax.set_title(f'Average 5400 m Thickness Contour\nDay {day_offset} After Event Dates', fontsize=14)

# #############################################
# # 6. Create and save the animation           #
# #############################################

# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=days_to_animate,
#     interval=1000,  # 1 second per frame (adjust as needed)
#     blit=False
# )

# # Save the animation using FFMpegWriter.
# from matplotlib.animation import FFMpegWriter
# import matplotlib
# # Set the path to your ffmpeg executable
# matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
# writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
# ani.save('avg_thickness_permafrost_animation.mp4', writer=writer)

# plt.show()

import rasterio
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#############################################
# 1. Load and subset the permafrost raster  #
#############################################

# Define the path to your .byte permafrost file
byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\PermaFrost\llipa.byte'

with rasterio.open(byte_file_path) as perma_ds:
    # Read the first band as a NumPy array.
    perma_data = perma_ds.read(1)
    
    # Determine pixel size in the y-direction (assumes north-up image)
    pixel_size_y = abs(perma_ds.transform.e)
    nrows = perma_data.shape[0]
    rows = np.arange(nrows)
    
    # Compute each row's center latitude.
    center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    
    # Find indices where the center latitude is above 25°N.
    valid_rows = np.where(center_lats > 25)[0]
    if valid_rows.size == 0:
        raise ValueError("No rows found with latitude above 25°N.")
    
    # Rows are ordered from north to south; select from the top to the last valid row.
    last_valid_row = valid_rows[-1]
    perma_subset = perma_data[:last_valid_row + 1, :]
    
    # Compute the new geographic extent from the raster bounds.
    new_left = perma_ds.bounds.left
    new_right = perma_ds.bounds.right
    new_top = perma_ds.bounds.top
    new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
    perma_extent = (new_left, new_right, new_bottom, new_top)

#############################################
# 2. Load geopotential height dataset (SSW) #
#############################################

# Load the SSW geopotential height dataset only
ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSWs_stereo.nc")

# Adjust longitudes from 0–360 to –180–180 if needed.
def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(
            longitude=(((ds.longitude + 180) % 360) - 180)
        )
        ds = ds.sortby('longitude')
    return ds

ds_ssw = adjust_longitude(ds_ssw)

#############################################
# 3. Define SSW event dates                  #
#############################################

# SSW event dates only
ssw_event_dates = [
    np.datetime64("2010-02-10", 'D'),
    np.datetime64("2010-03-24", 'D'),
    np.datetime64("2013-01-07", 'D'),
    np.datetime64("2018-02-12", 'D'),
    np.datetime64("2019-01-02", 'D'),
    np.datetime64("2021-01-05", 'D'),
    np.datetime64("2023-02-16", 'D')
]

# Combine SSW events into one list.
# Each event now includes its date and the SSW dataset.
all_events = []
for date in ssw_event_dates:
    all_events.append({'date': date, 'ds': ds_ssw})

#############################################
# 4. Prepare the figure and base map         #
#############################################

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent(perma_extent, crs=ccrs.PlateCarree())

# Add basic map features.
ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

# Plot the permafrost background image.
ax.imshow(
    perma_subset,
    transform=ccrs.PlateCarree(),
    extent=perma_extent,
    origin='upper',
    cmap='viridis',
    alpha=0.5
)

#############################################
# 5. Define animation update function       #
#############################################

# Define the total number of days (0 to 14 after the event)
total_days = 14
days_to_animate = np.arange(0, total_days + 1)

def update(day_offset):
    # Clear and reconfigure the axis.
    ax.clear()
    ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    
    # Re-plot the permafrost background.
    ax.imshow(
        perma_subset,
        transform=ccrs.PlateCarree(),
        extent=perma_extent,
        origin='upper',
        cmap='viridis',
        alpha=0.5
    )
    
    # List to hold thickness fields for each event at the given day offset.
    thickness_fields = []
    
    # Loop over all SSW events and compute the thickness field for the current day.
    for event in all_events:
        event_date = event['date'] + np.timedelta64(day_offset, 'D')
        ds = event['ds']
        
        try:
            # Select the nearest valid_time.
            geopotential = ds['z'].sel(valid_time=event_date, method='nearest')
        except Exception as e:
            print(f"No data for {event_date}: {e}")
            continue
        
        try:
            geopotential = geopotential.sel(
                pressure_level=[500.0, 1000.0],
                latitude=slice(new_top, new_bottom),  # high-to-low
                longitude=slice(new_left, new_right)
            )
        except Exception as e:
            print(f"Error selecting domain for {event_date}: {e}")
            continue
        
        if geopotential.isnull().all():
            print(f"All values are NaN for {event_date}")
            continue
        
        # Convert geopotential to geopotential height (in m).
        geopotential_height = geopotential / 9.80665
        
        try:
            hgt_500mb = geopotential_height.sel(pressure_level=500.0)
            hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)
        except Exception as e:
            print(f"Pressure levels issue for {event_date}: {e}")
            continue
        
        # Compute thickness as the difference between 500 mb and 1000 mb.
        thickness = hgt_500mb - hgt_1000mb
        thickness_fields.append(thickness)
    
    # Only proceed if we have at least one valid thickness field.
    if thickness_fields:
        # Concatenate along a new "event" dimension and compute the average.
        thickness_all = xr.concat(thickness_fields, dim='event')
        avg_thickness = thickness_all.mean(dim='event')
        
        # Plot the 5400 m thickness contour from the averaged thickness field.
        cs = ax.contour(
            avg_thickness.longitude,
            avg_thickness.latitude,
            avg_thickness,
            levels=[5400],
            colors='black',
            transform=ccrs.PlateCarree(),
            linewidths=2
        )
    
    ax.set_title(f'Average 5400 m Thickness Contour\nDay {day_offset} After SSW Event Dates', fontsize=14)

#############################################
# 6. Create and save the animation           #
#############################################

ani = animation.FuncAnimation(
    fig,
    update,
    frames=days_to_animate,
    interval=1000,  # 1 second per frame (adjust as needed)
    blit=False
)

# Save the animation using FFMpegWriter.
from matplotlib.animation import FFMpegWriter
import matplotlib
# Set the path to your ffmpeg executable
matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
ani.save('avg_thickness_permafrost_animation_SSW_only.mp4', writer=writer)

plt.show()
