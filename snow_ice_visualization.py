import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D

#############################################
# 1. Load ERA5 dataset (ciwc: cloud ice water) #
#############################################
# Adjust the file path as needed.
ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\snow_ice_content.nc")

def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
        ds = ds.sortby('longitude')
    return ds

ds_ssw = adjust_longitude(ds_ssw)

#############################################
# 2. Define domain based on ERA5 data         #
#############################################
# Use the ERA5 dataset's native domain
new_left  = float(ds_ssw.longitude.min())
new_right = float(ds_ssw.longitude.max())
new_top   = float(ds_ssw.latitude.max())
new_bottom= float(ds_ssw.latitude.min())
domain_extent = (new_left, new_right, new_bottom, new_top)

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
legend_elements = [Line2D([0], [0], color='red', lw=2, label='SSW Events')]
all_events = []
for date, label in zip(ssw_event_dates, ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': 'red', 'ds': ds_ssw})

#############################################
# 4. Precompute Climatology for Anomalies      #
#############################################
# Choose a pressure level (hPa) for your analysis.
pressure_level = 850.0

# Select the field over the full ERA5 domain.
field_all = ds_ssw['ciwc'].sel(
    pressure_level=pressure_level,
    latitude=slice(new_top, new_bottom),
    longitude=slice(new_left, new_right)
)
# Compute the climatology as the mean over the full valid_time period.
# (You could further restrict this to winter months if desired.)
climatology = field_all.mean(dim='valid_time')

#############################################
# 5. Prepare the figure, map, and colorbar     #
#############################################
fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent(domain_extent, crs=ccrs.PlateCarree())

def draw_base_map(ax):
    ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

draw_base_map(ax)

# Create a colorbar axis.
cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])

#############################################
# 6. Define animation update function        #
#############################################
total_days = 14
days_to_animate = np.arange(0, total_days + 1)  # Days 0 to 14

def update(day_offset):
    # Clear and redraw the base map.
    ax.clear()
    ax.set_extent(domain_extent, crs=ccrs.PlateCarree())
    draw_base_map(ax)
    
    # List to store each event's field for the current day offset.
    fields = []
    
    for event in all_events:
        event_date = event['date'] + np.timedelta64(day_offset, 'D')
        ds = event['ds']
        try:
            # Select the nearest available time.
            field = ds['ciwc'].sel(valid_time=event_date, method='nearest')
        except Exception as e:
            print(f"No data for event date {event_date}: {e}")
            continue
        try:
            field = field.sel(
                pressure_level=pressure_level,
                latitude=slice(new_top, new_bottom),
                longitude=slice(new_left, new_right)
            )
        except Exception as e:
            print(f"Error selecting data for event date {event_date}: {e}")
            continue
        if field.isnull().all():
            print(f"All values are NaN for event date: {event_date}")
            continue
        fields.append(field)
    
    if fields:
        # Compute the composite mean and standard deviation (variability).
        composite = xr.concat(fields, dim='event').mean(dim='event')
        composite_std = xr.concat(fields, dim='event').std(dim='event')
        
        # Compute the anomaly relative to the climatology.
        anomaly = composite - climatology
        
        # Set symmetric color limits based on the anomaly range.
        vmax = np.nanmax(np.abs(anomaly.values))
        vmin = -vmax
        
        # Plot the anomaly using a diverging colormap centered at zero.
        pm = ax.pcolormesh(
            composite.longitude,
            composite.latitude,
            anomaly,
            cmap='RdBu_r',
            vmin=vmin,
            vmax=vmax,
            transform=ccrs.PlateCarree(),
            alpha=0.8
        )
        
        # Optionally, overlay a contour for the composite standard deviation.
        cs = ax.contour(
            composite.longitude,
            composite.latitude,
            composite_std,
            levels=[np.nanmean(composite_std.values)],
            colors='k',
            transform=ccrs.PlateCarree(),
            linewidths=1.5
        )
        
        # Update the colorbar.
        cax.cla()
        plt.colorbar(pm, cax=cax)
    
    # Update the title.
    ax.set_title(f'Anomaly of Specific Cloud Ice Water Content at {pressure_level} hPa\n'
                 f'Day {day_offset} After SSW Events', fontsize=14)
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)

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

from matplotlib.animation import FFMpegWriter
import matplotlib
matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
ani.save('ciwc_anomaly_variability_SSW_no_permafrost.mp4', writer=writer)

plt.show()
