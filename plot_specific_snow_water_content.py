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
#     pixel_size_y = abs(perma_ds.transform.e)
#     nrows = perma_data.shape[0]
#     rows = np.arange(nrows)
#     center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
#     valid_rows = np.where(center_lats > 25)[0]
#     if valid_rows.size == 0:
#         raise ValueError("No rows found with latitude above 25°N.")
#     last_valid_row = valid_rows[-1]
#     perma_subset = perma_data[:last_valid_row + 1, :]
#     new_left = perma_ds.bounds.left
#     new_right = perma_ds.bounds.right
#     new_top = perma_ds.bounds.top
#     new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
#     perma_extent = (new_left, new_right, new_bottom, new_top)

# #############################################
# # 2. Load ERA5 pressure level dataset        #
# #############################################

# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\snow_ice_content.nc")

# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(
#             longitude=(((ds.longitude + 180) % 360) - 180)
#         )
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)

# #############################################
# # 3. Define SSW event dates, labels, and color#
# #############################################

# ssw_event_dates = [
#     np.datetime64("2010-02-10", 'D'),
#     np.datetime64("2010-03-24", 'D'),
#     np.datetime64("2013-01-07", 'D'),
#     np.datetime64("2018-02-12", 'D'),
#     np.datetime64("2019-01-02", 'D'),
#     np.datetime64("2021-01-05", 'D'),
#     np.datetime64("2023-02-16", 'D')
# ]
# ssw_event_labels = [
#     "SSW 2010-02-10",
#     "SSW 2010-03-24",
#     "SSW 2013-01-07",
#     "SSW 2018-02-12",
#     "SSW 2019-01-02",
#     "SSW 2021-01-05",
#     "SSW 2023-02-16"
# ]

# colors_ssw = ['red'] * len(ssw_event_dates)

# all_events = []
# for date, label in zip(ssw_event_dates, ssw_event_labels):
#     all_events.append({'date': date, 'label': label, 'color': 'red', 'ds': ds_ssw})

# #############################################
# # 4. Prepare the figure, axis, and colorbar  #
# #############################################

# fig = plt.figure(figsize=(12, 12))
# ax = plt.axes(projection=ccrs.NorthPolarStereo())
# ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
# ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
# gl.top_labels = False
# gl.right_labels = False

# legend_elements = [Line2D([0], [0], color='red', lw=2, label='SSW Events')]
# cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])

# #############################################
# # 5. Define animation parameters and update  #
# #############################################

# total_days = 14
# days_to_animate = np.arange(0, total_days + 1)
# pressure_level = 850.0  # Specify the pressure level (hPa) to plot

# def update(day_offset):
#     ax.clear()
#     ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
#     ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
#     ax.add_feature(cfeature.COASTLINE)
#     ax.add_feature(cfeature.BORDERS, linestyle=':')
#     gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False

#     # Plot permafrost background.
#     ax.imshow(
#         perma_subset,
#         transform=ccrs.PlateCarree(),
#         extent=perma_extent,
#         origin='upper',
#         cmap='viridis',
#         alpha=0.5
#     )
#     ax.set_title(f'Specific Snow Water Content at {pressure_level} hPa\nDay {day_offset} After SSW Events', fontsize=14)
    
#     mappables = []
#     for event in all_events:
#         event_date = event['date'] + np.timedelta64(day_offset, 'D')
#         ds = event['ds']
#         try:
#             # Select the specific snow water content.
#             field = ds['cswc'].sel(valid_time=event_date, method='nearest')
#         except KeyError:
#             print(f"No data available for event date: {event_date}")
#             continue
#         try:
#             field = field.sel(
#                 pressure_level=pressure_level,
#                 latitude=slice(new_top, new_bottom),
#                 longitude=slice(new_left, new_right)
#             )
#         except Exception as e:
#             print(f"Error selecting data for event date {event_date}: {e}")
#             continue
#         if field.isnull().all():
#             print(f"No data available (all NaN) for event date: {event_date}")
#             continue
#         # Plot filled contours of cswc.
#         cf = ax.contourf(
#             field.longitude,
#             field.latitude,
#             field,
#             levels=20,
#             cmap='viridis',
#             transform=ccrs.PlateCarree(),
#             alpha=0.8
#         )
#         mappables.append(cf)
#     ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)
#     if mappables:
#         cax.cla()
#         plt.colorbar(mappables[0], cax=cax)

# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=days_to_animate,
#     interval=1000,
#     blit=False
# )

# from matplotlib.animation import FFMpegWriter
# import matplotlib
# matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
# writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
# ani.save('cswc_permafrost_animation_SSW.mp4', writer=writer)

# plt.show()

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
#     pixel_size_y = abs(perma_ds.transform.e)
#     nrows = perma_data.shape[0]
#     rows = np.arange(nrows)
#     center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
#     valid_rows = np.where(center_lats > 25)[0]
#     if valid_rows.size == 0:
#         raise ValueError("No rows found with latitude above 25°N.")
#     last_valid_row = valid_rows[-1]
#     perma_subset = perma_data[:last_valid_row + 1, :]
#     new_left = perma_ds.bounds.left
#     new_right = perma_ds.bounds.right
#     new_top = perma_ds.bounds.top
#     new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
#     perma_extent = (new_left, new_right, new_bottom, new_top)

# #############################################
# # 2. Load ERA5 pressure level dataset        #
# #############################################
# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\snow_ice_content.nc")

# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)

# #############################################
# # 3. Define SSW event dates, labels, and color#
# #############################################
# ssw_event_dates = [
#     np.datetime64("2010-02-10", 'D'),
#     np.datetime64("2010-03-24", 'D'),
#     np.datetime64("2013-01-07", 'D'),
#     np.datetime64("2018-02-12", 'D'),
#     np.datetime64("2019-01-02", 'D'),
#     np.datetime64("2021-01-05", 'D'),
#     np.datetime64("2023-02-16", 'D')
# ]
# ssw_event_labels = [
#     "SSW 2010-02-10",
#     "SSW 2010-03-24",
#     "SSW 2013-01-07",
#     "SSW 2018-02-12",
#     "SSW 2019-01-02",
#     "SSW 2021-01-05",
#     "SSW 2023-02-16"
# ]
# colors_ssw = ['red'] * len(ssw_event_dates)
# all_events = []
# for date, label in zip(ssw_event_dates, ssw_event_labels):
#     all_events.append({'date': date, 'label': label, 'color': 'red', 'ds': ds_ssw})

# #############################################
# # 4. Prepare the figure and animation axis   #
# #############################################
# fig = plt.figure(figsize=(12, 12))
# ax = plt.axes(projection=ccrs.NorthPolarStereo())
# ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
# ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
# gl.top_labels = False
# gl.right_labels = False

# legend_elements = [Line2D([0], [0], color='red', lw=2, label='SSW Events')]
# cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])

# #############################################
# # 5. Define animation parameters and update  #
# #############################################
# total_days = 14
# days_to_animate = np.arange(0, total_days + 1)
# pressure_level = 850.0  # in hPa

# def update(day_offset):
#     ax.clear()
#     ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
#     ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
#     ax.add_feature(cfeature.COASTLINE)
#     ax.add_feature(cfeature.BORDERS, linestyle=':')
#     gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False

#     # Plot permafrost overlay
#     ax.imshow(
#         perma_subset,
#         transform=ccrs.PlateCarree(),
#         extent=perma_extent,
#         origin='upper',
#         cmap='viridis',
#         alpha=0.5
#     )
#     ax.set_title(f'Specific Snow Water Content at {pressure_level} hPa\nDay {day_offset} After SSW Events', fontsize=14)

#     # Loop through events and plot the cswc field with pcolormesh
#     for event in all_events:
#         event_date = event['date'] + np.timedelta64(day_offset, 'D')
#         ds = event['ds']
#         try:
#             field = ds['cswc'].sel(valid_time=event_date, method='nearest')
#         except KeyError:
#             print(f"No data available for event date: {event_date}")
#             continue
#         try:
#             field = field.sel(
#                 pressure_level=pressure_level,
#                 latitude=slice(new_top, new_bottom),
#                 longitude=slice(new_left, new_right)
#             )
#         except Exception as e:
#             print(f"Error selecting data for event date {event_date}: {e}")
#             continue
#         if field.isnull().all():
#             print(f"No data available (all NaN) for event date: {event_date}")
#             continue

#         pm = ax.pcolormesh(
#             field.longitude,
#             field.latitude,
#             field,
#             cmap='viridis',
#             transform=ccrs.PlateCarree(),
#             alpha=0.7
#         )
#     ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)
#     if pm is not None:
#         cax.cla()
#         plt.colorbar(pm, cax=cax)

# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=days_to_animate,
#     interval=1000,
#     blit=False
# )

# from matplotlib.animation import FFMpegWriter
# import matplotlib
# matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
# writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
# ani.save('cswc_permafrost_animation_SSW3.mp4', writer=writer)

# plt.show()


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
    # Get the pixel size (assumes north-up image)
    pixel_size_y = abs(perma_ds.transform.e)
    nrows = perma_data.shape[0]
    rows = np.arange(nrows)
    # Compute center latitudes for each row
    center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    valid_rows = np.where(center_lats > 25)[0]
    if valid_rows.size == 0:
        raise ValueError("No rows found with latitude above 25°N.")
    last_valid_row = valid_rows[-1]
    perma_subset = perma_data[:last_valid_row + 1, :]
    # Recalculate geographic extent
    new_left = perma_ds.bounds.left
    new_right = perma_ds.bounds.right
    new_top = perma_ds.bounds.top
    new_bottom = perma_ds.bounds.top - (last_valid_row + 1) * pixel_size_y
    perma_extent = (new_left, new_right, new_bottom, new_top)

#############################################
# 2. Load ERA5 pressure level dataset        #
#############################################
ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\snow_ice_content.nc")

def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
        ds = ds.sortby('longitude')
    return ds

ds_ssw = adjust_longitude(ds_ssw)

#############################################
# 3. Define SSW event dates, labels, and color #
#############################################
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
colors_ssw = ['red'] * len(ssw_event_dates)
all_events = []
for date, label in zip(ssw_event_dates, ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': 'red', 'ds': ds_ssw})

#############################################
# 4. Prepare the figure and animation axis   #
#############################################
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

# Legend for SSW events
legend_elements = [Line2D([0], [0], color='red', lw=2, label='SSW Events')]

# (Optional) Create an axis for the colorbar if desired
cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])

#############################################
# 5. Define animation parameters and update  #
#############################################
total_days = 14
days_to_animate = np.arange(0, total_days + 1)  # Days 0 to 14
pressure_level = 850.0  # in hPa

def update(day_offset):
    ax.clear()
    ax.set_extent(perma_extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

    # Plot the permafrost background (using imshow as in your provided code)
    ax.imshow(
        perma_subset,
        transform=ccrs.PlateCarree(),
        extent=perma_extent,
        origin='upper',
        cmap='viridis',
        alpha=0.5
    )
    ax.set_title(f'Specific Cloud Ice Water Content at {pressure_level} hPa\nDay {day_offset} After SSW Events', fontsize=14)

    # For each SSW event, select the parameter field and plot it using pcolormesh
    for event in all_events:
        event_date = event['date'] + np.timedelta64(day_offset, 'D')
        ds = event['ds']
        try:
            field = ds['ciwc'].sel(valid_time=event_date, method='nearest')
        except KeyError:
            print(f"No data available for event date: {event_date}")
            continue
        try:
            field = field.sel(
                pressure_level=pressure_level,
                latitude=slice(new_top, new_bottom),  # high (top) to low (bottom)
                longitude=slice(new_left, new_right)
            )
        except Exception as e:
            print(f"Error selecting data for event date {event_date}: {e}")
            continue
        if field.isnull().all():
            print(f"No data available (all NaN) for event date: {event_date}")
            continue

        # Plot the parameter field using pcolormesh.
        # Adjust alpha as needed so that the permafrost overlay remains visible.
        pm = ax.pcolormesh(
            field.longitude,
            field.latitude,
            field,
            cmap='viridis',
            transform=ccrs.PlateCarree(),
            alpha=0.7
        )
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)
    # Update the colorbar (if using a separate axis)
    if pm is not None:
        cax.cla()
        plt.colorbar(pm, cax=cax)

ani = animation.FuncAnimation(
    fig,
    update,
    frames=days_to_animate,
    interval=1000,
    blit=False
)

# (Optional) Save the animation to a video file.
from matplotlib.animation import FFMpegWriter
import matplotlib
matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
ani.save('ciwc_permafrost_animation5_SSW.mp4', writer=writer)

plt.show()
