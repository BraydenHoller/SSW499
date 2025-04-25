# import rasterio
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from matplotlib.lines import Line2D

# #############################################
# # 1. Load and subset the permafrost raster  #
# #############################################

# byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\PermaFrost\llipa.byte'

# with rasterio.open(byte_file_path) as perma_ds:
#     perma_data = perma_ds.read(1)
    
#     # Determine the pixel size in the y-direction (assumes north-up image)
#     pixel_size_y = abs(perma_ds.transform.e)
#     nrows = perma_data.shape[0]
#     rows = np.arange(nrows)
    
#     # Compute each row's center latitude
#     center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    
#     # Find indices where the center latitude is above 25°N
#     valid_rows = np.where(center_lats > 25)[0]
#     if valid_rows.size == 0:
#         raise ValueError("No rows found with latitude above 25°N.")
    
#     # Subset rows from the top down to the last valid row
#     last_valid_row = valid_rows[-1]
#     perma_subset = perma_data[:last_valid_row + 1, :]
    
#     # Compute the new geographic extent from the raster bounds
#     new_left = perma_ds.bounds.left
#     new_right = perma_ds.bounds.right
#     new_top = perma_ds.bounds.top
#     new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
#     perma_extent = (new_left, new_right, new_bottom, new_top)

# #############################################
# # 2. Load geopotential height datasets      #
# #############################################

# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSWs_stereo.nc")
# ds_non_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\non_SSWs_stereo.nc")

# # Adjust longitudes from 0–360 to –180–180 if needed
# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)
# ds_non_ssw = adjust_longitude(ds_non_ssw)

# #############################################
# # 3. Define event dates and combine         #
# #############################################

# # SSW events
# ssw_event_dates = [
#     np.datetime64("2010-02-10", 'D'),
#     np.datetime64("2010-03-24", 'D'),
#     np.datetime64("2013-01-07", 'D'),
#     np.datetime64("2018-02-12", 'D'),
#     np.datetime64("2019-01-02", 'D'),
#     np.datetime64("2021-01-05", 'D'),
#     np.datetime64("2023-02-16", 'D')
# ]

# # Non-SSW events
# non_ssw_event_dates = [
#     np.datetime64("2012-01-05", 'D'),
#     np.datetime64("2014-01-05", 'D'),
#     np.datetime64("2015-01-05", 'D'),
#     np.datetime64("2016-01-05", 'D'),
#     np.datetime64("2017-01-05", 'D'),
#     np.datetime64("2020-01-05", 'D'),
#     np.datetime64("2022-01-05", 'D')
# ]

# # Combine them all
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

# # Add basic map features
# ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
# gl.top_labels = False
# gl.right_labels = False

# # Plot the permafrost background
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

# total_days = 14
# days_to_animate = np.arange(0, total_days + 1)

# def update(day_offset):
#     ax.clear()
#     ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
    
#     # Re-add map features
#     ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
#     ax.add_feature(cfeature.COASTLINE)
#     ax.add_feature(cfeature.BORDERS, linestyle=':')
#     gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False
    
#     # Permafrost background
#     ax.imshow(
#         perma_subset,
#         transform=ccrs.PlateCarree(),
#         extent=perma_extent,
#         origin='upper',
#         cmap='viridis',
#         alpha=0.5
#     )
    
#     # Gather thickness fields for all events at this day offset
#     thickness_fields = []
    
#     for event in all_events:
#         event_date = event['date'] + np.timedelta64(day_offset, 'D')
#         ds = event['ds']
        
#         # Attempt to select the nearest valid_time
#         try:
#             geopotential = ds['z'].sel(valid_time=event_date, method='nearest')
#         except Exception:
#             # No data for this date
#             continue
        
#         # Subset domain, select pressure levels 500 and 1000
#         try:
#             geopotential = geopotential.sel(
#                 pressure_level=[500.0, 1000.0],
#                 latitude=slice(new_top, new_bottom),
#                 longitude=slice(new_left, new_right)
#             )
#         except Exception:
#             continue
        
#         # Skip if all NaNs
#         if geopotential.isnull().all():
#             continue
        
#         # Convert geopotential to geopotential height (m)
#         geopotential_height = geopotential / 9.80665
        
#         # Extract heights
#         try:
#             hgt_500mb = geopotential_height.sel(pressure_level=500.0)
#             hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)
#         except Exception:
#             continue
        
#         # Thickness
#         thickness = hgt_500mb - hgt_1000mb
#         thickness_fields.append(thickness)
    
#     # Only proceed if we have valid thickness fields
#     if len(thickness_fields) > 0:
#         # Concatenate and compute the 25th, 50th (median), and 75th percentile
#         thickness_all = xr.concat(thickness_fields, dim='event')
#         q25 = thickness_all.quantile(0.25, dim='event')
#         q50 = thickness_all.quantile(0.50, dim='event')
#         q75 = thickness_all.quantile(0.75, dim='event')
        
#         # Plot the 5400-m contour for each percentile
#         cs25 = ax.contour(
#             q25.longitude,
#             q25.latitude,
#             q25,
#             levels=[5400],
#             colors='blue',
#             transform=ccrs.PlateCarree(),
#             linewidths=1.5
#         )
        
#         cs50 = ax.contour(
#             q50.longitude,
#             q50.latitude,
#             q50,
#             levels=[5400],
#             colors='green',
#             transform=ccrs.PlateCarree(),
#             linewidths=1.5
#         )
        
#         cs75 = ax.contour(
#             q75.longitude,
#             q75.latitude,
#             q75,
#             levels=[5400],
#             colors='red',
#             transform=ccrs.PlateCarree(),
#             linewidths=1.5
#         )
    
#     # Title
#     ax.set_title(f'5400-m Thickness IQR\nDay {day_offset} After Event Dates', fontsize=14)

#     # Optional: Add a legend for the lines
#     handles = [
#         Line2D([0], [0], color='blue', label='25th Percentile'),
#         Line2D([0], [0], color='green', label='Median'),
#         Line2D([0], [0], color='red', label='75th Percentile')
#     ]
#     ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)

# #############################################
# # 6. Create and save the animation           #
# #############################################

# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=days_to_animate,
#     interval=1000,  # 1 second per frame
#     blit=False
# )

# # Save the animation with FFMpegWriter
# from matplotlib.animation import FFMpegWriter
# import matplotlib

# # Adjust if your ffmpeg path is different
# matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
# writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)

# ani.save('thickness_IQR_animation.mp4', writer=writer)
# plt.show()

# Computationally Inefficient Code
# import rasterio
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from matplotlib.lines import Line2D

# #############################################
# # 1. Load and subset the permafrost raster  #
# #############################################

# byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\PermaFrost\llipa.byte'

# with rasterio.open(byte_file_path) as perma_ds:
#     perma_data = perma_ds.read(1)
    
#     # Determine the pixel size in the y-direction (assumes north-up image)
#     pixel_size_y = abs(perma_ds.transform.e)
#     nrows = perma_data.shape[0]
#     rows = np.arange(nrows)
    
#     # Compute each row's center latitude
#     center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    
#     # Find indices where the center latitude is above 25°N
#     valid_rows = np.where(center_lats > 25)[0]
#     if valid_rows.size == 0:
#         raise ValueError("No rows found with latitude above 25°N.")
    
#     # Subset rows from the top down to the last valid row
#     last_valid_row = valid_rows[-1]
#     perma_subset = perma_data[:last_valid_row + 1, :]
    
#     # Compute the new geographic extent from the raster bounds
#     new_left = perma_ds.bounds.left
#     new_right = perma_ds.bounds.right
#     new_top = perma_ds.bounds.top
#     new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
#     perma_extent = (new_left, new_right, new_bottom, new_top)

# #############################################
# # 2. Load geopotential height dataset (SSW) #
# #############################################

# # Load the SSW geopotential height dataset only
# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All")

# # Adjust longitudes from 0–360 to –180–180 if needed
# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)

# #############################################
# # 3. Define ONLY SSW event dates            #
# #############################################

# # Since 2010
# # ssw_event_dates = [
# #     np.datetime64("2010-02-10", 'D'),
# #     np.datetime64("2010-03-24", 'D'),
# #     np.datetime64("2013-01-06", 'D'),
# #     np.datetime64("2018-02-12", 'D'),
# #     np.datetime64("2019-01-01", 'D'),
# #     np.datetime64("2021-01-05", 'D'),
# #     np.datetime64("2023-02-16", 'D')
# # ]

# # # 2000-2010
# # ssw_event_dates = [
# #     np.datetime64("2000-03-20", 'D'),
# #     np.datetime64("2001-02-11", 'D'),
# #     np.datetime64("2001-12-30", 'D'),
# #     np.datetime64("2002-02-17", 'D'),
# #     np.datetime64("2003-01-18", 'D'),
# #     np.datetime64("2004-01-05", 'D'),
# #     np.datetime64("2006-01-21", 'D'),
# #     np.datetime64("2007-02-24", 'D'),
# #     np.datetime64("2008-02-22", 'D'),
# #     np.datetime64("2009-01-24", 'D')
# # ]

# # 1990-2000
# ssw_event_dates = [
#     np.datetime64("1998-12-15", 'D'),
#     np.datetime64("1999-02-26", 'D')
# ]

# # # 1980-1990
# # ssw_event_dates = [
# #     np.datetime64("1980-02-29", 'D'),
# #     np.datetime64("1981-03-04", 'D'),
# #     np.datetime64("1981-12-04", 'D'),
# #     np.datetime64("1984-02-24", 'D'),
# #     np.datetime64("1985-01-01", 'D'),
# #     np.datetime64("1987-01-23", 'D'),
# #     np.datetime64("1987-12-08", 'D'),
# #     np.datetime64("1988-03-14", 'D'),
# #     np.datetime64("1989-02-21", 'D')
# # ]

# # # 1970-1980
# # ssw_event_dates = [
# #     np.datetime64("1970-01-02", 'D'),
# #     np.datetime64("1971-01-18", 'D'),
# #     np.datetime64("1971-03-20", 'D'),
# #     np.datetime64("1973-01-31", 'D'),
# #     np.datetime64("1977-01-09", 'D'),
# #     np.datetime64("1979-02-22", 'D')
# # ]

# # # 1958-1970
# # ssw_event_dates = [
# #     np.datetime64("1958-02-08", 'D'),
# #     np.datetime64("1960-01-17", 'D'),
# #     np.datetime64("1963-01-27", 'D'),
# #     np.datetime64("1965-12-16", 'D'),
# #     np.datetime64("1966-02-22", 'D'),
# #     np.datetime64("1968-01-07", 'D'),
# #     np.datetime64("1968-11-28", 'D'),
# #     np.datetime64("1969-03-13", 'D')
# # ]

# # Combine them into one list. Each item includes the date and the SSW dataset.
# all_events = []
# for date in ssw_event_dates:
#     all_events.append({'date': date, 'ds': ds_ssw})

# #############################################
# # 4. Prepare the figure and base map         #
# #############################################

# fig = plt.figure(figsize=(12, 12))
# ax = plt.axes(projection=ccrs.NorthPolarStereo())
# ax.set_extent(perma_extent, crs=ccrs.PlateCarree())

# # Add basic map features
# ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
# gl.top_labels = False
# gl.right_labels = False

# # Plot the permafrost background
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

# total_days = 14
# days_to_animate = np.arange(0, total_days + 1)

# def update(day_offset):
#     ax.clear()
#     ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
    
#     # Re-add map features
#     ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
#     ax.add_feature(cfeature.COASTLINE)
#     ax.add_feature(cfeature.BORDERS, linestyle=':')
#     gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False
    
#     # Permafrost background
#     ax.imshow(
#         perma_subset,
#         transform=ccrs.PlateCarree(),
#         extent=perma_extent,
#         origin='upper',
#         cmap='viridis',
#         alpha=0.5
#     )
    
#     # Gather thickness fields for all SSW events at this day offset
#     thickness_fields = []
    
#     for event in all_events:
#         event_date = event['date'] + np.timedelta64(day_offset, 'D')
#         ds = event['ds']
        
#         # Attempt to select the nearest valid_time
#         try:
#             geopotential = ds['z'].sel(valid_time=event_date, method='nearest')
#         except Exception:
#             # No data for this date
#             continue
        
#         # Subset domain, select pressure levels 500 and 1000
#         try:
#             geopotential = geopotential.sel(
#                 pressure_level=[500.0, 1000.0],
#                 latitude=slice(new_top, new_bottom),
#                 longitude=slice(new_left, new_right)
#             )
#         except Exception:
#             continue
        
#         # Skip if all NaNs
#         if geopotential.isnull().all():
#             continue
        
#         # Convert geopotential to geopotential height (m)
#         geopotential_height = geopotential / 9.80665
        
#         # Extract heights
#         try:
#             hgt_500mb = geopotential_height.sel(pressure_level=500.0)
#             hgt_1000mb = geopotential_height.sel(pressure_level=1000.0)
#         except Exception:
#             continue
        
#         # Thickness
#         thickness = hgt_500mb - hgt_1000mb
#         thickness_fields.append(thickness)
    
#     # Only proceed if we have valid thickness fields
#     if len(thickness_fields) > 0:
#         # Concatenate and compute the 25th, 50th (median), and 75th percentile
#         thickness_all = xr.concat(thickness_fields, dim='event')
#         q25 = thickness_all.quantile(0.25, dim='event')
#         q50 = thickness_all.quantile(0.50, dim='event')
#         q75 = thickness_all.quantile(0.75, dim='event')
        
#         # Plot the 5400-m contour for each percentile
#         cs25 = ax.contour(
#             q25.longitude,
#             q25.latitude,
#             q25,
#             levels=[5400],
#             colors='blue',
#             transform=ccrs.PlateCarree(),
#             linewidths=1.5
#         )
        
#         cs50 = ax.contour(
#             q50.longitude,
#             q50.latitude,
#             q50,
#             levels=[5400],
#             colors='green',
#             transform=ccrs.PlateCarree(),
#             linewidths=1.5
#         )
        
#         cs75 = ax.contour(
#             q75.longitude,
#             q75.latitude,
#             q75,
#             levels=[5400],
#             colors='red',
#             transform=ccrs.PlateCarree(),
#             linewidths=1.5
#         )
    
#     # Title
#     ax.set_title(f'SSW-Only 5400-m Thickness IQR\nDay {day_offset} After Event Dates', fontsize=14)

#     # Optional: Add a legend for the lines
#     handles = [
#         Line2D([0], [0], color='blue', label='25th Percentile'),
#         Line2D([0], [0], color='green', label='Median'),
#         Line2D([0], [0], color='red', label='75th Percentile')
#     ]
#     ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)

# #############################################
# # 6. Create and save the animation           #
# #############################################

# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=days_to_animate,
#     interval=1000,  # 1 second per frame
#     blit=False
# )

# # Save the animation with FFMpegWriter
# from matplotlib.animation import FFMpegWriter
# import matplotlib

# matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
# writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)

# ani.save('thickness_IQR_animation_SSW_2000-2010.mp4', writer=writer)
# plt.show()

# More Efficient GPT-O3 Mini-High
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

byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\PermaFrost\llipa.byte'

with rasterio.open(byte_file_path) as perma_ds:
    perma_data = perma_ds.read(1)
    
    # Determine the pixel size in the y-direction (assumes north-up image)
    pixel_size_y = abs(perma_ds.transform.e)
    nrows = perma_data.shape[0]
    rows = np.arange(nrows)
    
    # Compute each row's center latitude
    center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    
    # Find indices where the center latitude is above 25°N
    valid_rows = np.where(center_lats > 25)[0]
    if valid_rows.size == 0:
        raise ValueError("No rows found with latitude above 25°N.")
    
    # Subset rows from the top down to the last valid row
    last_valid_row = valid_rows[-1]
    perma_subset = perma_data[:last_valid_row + 1, :]
    
    # Compute the new geographic extent from the raster bounds
    new_left = perma_ds.bounds.left
    new_right = perma_ds.bounds.right
    new_top = perma_ds.bounds.top
    new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
    perma_extent = (new_left, new_right, new_bottom, new_top)

#############################################
# 2. Load geopotential height dataset (SSW) #
#############################################

ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All")

def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
        ds = ds.sortby('longitude')
    return ds

ds_ssw = adjust_longitude(ds_ssw)

#############################################
# 3. Define SSW event dates (2000-2010)       #
#############################################

# Since 2010
ssw_event_dates = [
    np.datetime64("2010-02-10", 'D'),
    np.datetime64("2010-03-24", 'D'),
    np.datetime64("2013-01-06", 'D'),
    np.datetime64("2018-02-12", 'D'),
    np.datetime64("2019-01-01", 'D'),
    np.datetime64("2021-01-05", 'D'),
    np.datetime64("2023-02-16", 'D')
]

# # 2000-2010
# ssw_event_dates = [
#     np.datetime64("2000-03-20", 'D'),
#     np.datetime64("2001-02-11", 'D'),
#     np.datetime64("2001-12-30", 'D'),
#     np.datetime64("2002-02-17", 'D'),
#     np.datetime64("2003-01-18", 'D'),
#     np.datetime64("2004-01-05", 'D'),
#     np.datetime64("2006-01-21", 'D'),
#     np.datetime64("2007-02-24", 'D'),
#     np.datetime64("2008-02-22", 'D'),
#     np.datetime64("2009-01-24", 'D')
# ]

# # 1990-2000
# ssw_event_dates = [
#     np.datetime64("1998-12-15", 'D'),
#     np.datetime64("1999-02-26", 'D')
# ]

# # 1980-1990
# ssw_event_dates = [
#     np.datetime64("1980-02-29", 'D'),
#     np.datetime64("1981-03-04", 'D'),
#     np.datetime64("1981-12-04", 'D'),
#     np.datetime64("1984-02-24", 'D'),
#     np.datetime64("1985-01-01", 'D'),
#     np.datetime64("1987-01-23", 'D'),
#     np.datetime64("1987-12-08", 'D'),
#     np.datetime64("1988-03-14", 'D'),
#     np.datetime64("1989-02-21", 'D')
# ]

# # 1970-1980
# ssw_event_dates = [
#     np.datetime64("1970-01-02", 'D'),
#     np.datetime64("1971-01-18", 'D'),
#     np.datetime64("1971-03-20", 'D'),
#     np.datetime64("1973-01-31", 'D'),
#     np.datetime64("1977-01-09", 'D'),
#     np.datetime64("1979-02-22", 'D')
# ]

# # 1958-1970
# ssw_event_dates = [
#     np.datetime64("1958-02-08", 'D'),
#     np.datetime64("1960-01-17", 'D'),
#     np.datetime64("1963-01-27", 'D'),
#     np.datetime64("1965-12-16", 'D'),
#     np.datetime64("1966-02-22", 'D'),
#     np.datetime64("1968-01-07", 'D'),
#     np.datetime64("1968-11-28", 'D'),
#     np.datetime64("1969-03-13", 'D')
# ]

all_events = [{'date': date, 'ds': ds_ssw} for date in ssw_event_dates]

#############################################
# 4. Precompute thickness percentiles       #
#############################################

# Pre-subset the SSW data spatially and by pressure level for efficiency.
ds_subset = ds_ssw['z'].sel(
    pressure_level=[500.0, 1000.0],
    latitude=slice(new_top, new_bottom),
    longitude=slice(new_left, new_right)
)

total_days = 14
days_to_animate = np.arange(0, total_days + 1)
precomputed_percentiles = {}

for d in days_to_animate:
    thickness_fields = []
    for event in all_events:
        event_date = event['date'] + np.timedelta64(d, 'D')
        try:
            # Use the pre-subset data and select the nearest valid_time
            geopotential = ds_subset.sel(valid_time=event_date, method='nearest')
        except Exception:
            continue
        
        if geopotential.isnull().all():
            continue
        
        # Convert geopotential to geopotential height (m)
        geopotential_height = geopotential / 9.80665
        
        try:
            hgt_500 = geopotential_height.sel(pressure_level=500.0)
            hgt_1000 = geopotential_height.sel(pressure_level=1000.0)
        except Exception:
            continue
        
        # Calculate thickness
        thickness = hgt_500 - hgt_1000
        thickness_fields.append(thickness)
    
    if thickness_fields:
        # Concatenate along a new 'event' dimension and compute percentiles
        thickness_all = xr.concat(thickness_fields, dim='event')
        q25 = thickness_all.quantile(0.25, dim='event')
        q50 = thickness_all.quantile(0.50, dim='event')
        q75 = thickness_all.quantile(0.75, dim='event')
        precomputed_percentiles[d] = (q25, q50, q75)
    else:
        precomputed_percentiles[d] = None

#############################################
# 5. Prepare the figure and static base map  #
#############################################

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())

def draw_static_background(ax):
    ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    ax.imshow(
        perma_subset,
        transform=ccrs.PlateCarree(),
        extent=perma_extent,
        origin='upper',
        cmap='viridis',
        alpha=0.5
    )

# Draw the static base map once.
draw_static_background(ax)

# Add a static legend.
legend_handles = [
    Line2D([0], [0], color='blue', label='25th Percentile'),
    Line2D([0], [0], color='green', label='Median'),
    Line2D([0], [0], color='red', label='75th Percentile')
]
ax.legend(handles=legend_handles, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)

#############################################
# 6. Define animation update function       #
#############################################

# We'll keep track of dynamic (contour) artists so we can remove them on each frame.
dynamic_artists = []

def update(day_offset):
    global dynamic_artists
    # Remove any previously drawn dynamic contour artists.
    for artist in dynamic_artists:
        try:
            artist.remove()
        except Exception:
            pass
    dynamic_artists = []
    
    # Plot precomputed contours if data is available for this day offset.
    percentiles = precomputed_percentiles.get(day_offset)
    if percentiles is not None:
        q25, q50, q75 = percentiles
        cs25 = ax.contour(
            q25.longitude,
            q25.latitude,
            q25,
            levels=[5400],
            colors='blue',
            transform=ccrs.PlateCarree(),
            linewidths=1.5
        )
        cs50 = ax.contour(
            q50.longitude,
            q50.latitude,
            q50,
            levels=[5400],
            colors='green',
            transform=ccrs.PlateCarree(),
            linewidths=1.5
        )
        cs75 = ax.contour(
            q75.longitude,
            q75.latitude,
            q75,
            levels=[5400],
            colors='red',
            transform=ccrs.PlateCarree(),
            linewidths=1.5
        )
        # Store the contour collections so they can be removed in the next update.
        dynamic_artists.extend(cs25.collections)
        dynamic_artists.extend(cs50.collections)
        dynamic_artists.extend(cs75.collections)
    
    # Update the title to reflect the current day offset.
    ax.set_title(f'SSW-Only 5400-m Thickness IQR\nDay {day_offset} After Event Dates', fontsize=14)
    
    return dynamic_artists  # optional, if blitting or further handling is desired

#############################################
# 7. Create and save the animation           #
#############################################

ani = animation.FuncAnimation(
    fig,
    update,
    frames=days_to_animate,
    interval=1000,  # 1 second per frame
    blit=False
)

# Save the animation with FFMpegWriter
from matplotlib.animation import FFMpegWriter
import matplotlib

matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)

ani.save('thickness_IQR_animation_SSW_1958-1970.mp4', writer=writer)
plt.show()
