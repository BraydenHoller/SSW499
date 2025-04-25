# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from matplotlib.lines import Line2D

# # Load the consolidated dataset
# ds = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Since_2010_Temps_1000_500.nc")

# # Adjust longitudes from 0-360 to -180 to 180 if necessary
# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(
#             longitude=(((ds.longitude + 180) % 360) - 180)
#         )
#         ds = ds.sortby('longitude')
#     return ds

# ds = adjust_longitude(ds)

# # Define event dates and labels for both SSW and Non-SSW events
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
# ssw_event_labels = [
#     "SSW 2010-02-10",
#     "SSW 2010-03-24",
#     "SSW 2013-01-07",
#     "SSW 2018-02-12",
#     "SSW 2019-01-02",
#     "SSW 2021-01-05",
#     "SSW 2023-02-16"
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
# non_ssw_event_labels = [
#     "Non-SSW 2012-01-05",
#     "Non-SSW 2014-01-05",
#     "Non-SSW 2015-01-05",
#     "Non-SSW 2016-01-05",
#     "Non-SSW 2017-01-05",
#     "Non-SSW 2020-01-05",
#     "Non-SSW 2022-01-05"
# ]

# # Define the latitude and longitude boundaries
# lat_min, lat_max = 20, 90  # From 20 degrees latitude up to the North Pole
# lon_min, lon_max = -180, 180  # Global longitude coverage

# # Colors for SSW and Non-SSW events
# color_ssw = 'red'
# color_non_ssw = 'navy'

# # Prepare the figure and axes
# fig = plt.figure(figsize=(14, 14))
# ax = plt.axes(projection=ccrs.NorthPolarStereo())
# ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# # Add map features
# ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
# gl.top_labels = False
# gl.right_labels = False

# # Collect legend elements
# legend_elements = [
#     Line2D([0], [0], color=color_ssw, lw=2, label='SSW Events'),
#     Line2D([0], [0], color=color_non_ssw, lw=2, label='Non-SSW Events')
# ]

# # Combine SSW and Non-SSW events into one list with their corresponding colors
# all_events = []
# for date, label in zip(ssw_event_dates, ssw_event_labels):
#     all_events.append({'date': date, 'label': label, 'color': color_ssw})

# for date, label in zip(non_ssw_event_dates, non_ssw_event_labels):
#     all_events.append({'date': date, 'label': label, 'color': color_non_ssw})

# # Define the total number of days to animate
# total_days = 14  # Animate for 14 days after event dates

# # Create a list of days to animate
# days_to_animate = np.arange(0, total_days + 1)  # From day 0 to day 14

# # Define the cities to label with their coordinates
# cities = [
#     {'name': 'Washington D.C.', 'lat': 38.9072, 'lon': -77.0369},
#     {'name': 'Beijing', 'lat': 39.9042, 'lon': 116.4074},
#     {'name': 'Moscow', 'lat': 55.7558, 'lon': 37.6173}
# ]

# def update(day_offset):
#     ax.clear()
#     ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

#     # Add map features
#     ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
#     ax.add_feature(cfeature.COASTLINE)
#     ax.add_feature(cfeature.BORDERS, linestyle=':')
#     gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False

#     # Title for the current frame
#     current_day = day_offset
#     ax.set_title(f'0°C Isotherm at 1000 mb - Day {current_day} After Event Dates')

#     # Plot contours for all events at the current day offset
#     for event in all_events:
#         event_date = event['date'] + np.timedelta64(day_offset, 'D')
#         color = event['color']

#         # Select the valid_time using method='nearest'
#         try:
#             temperature = ds['t'].sel(valid_time=event_date, method='nearest')
#         except KeyError:
#             print(f"No data available for event date: {event_date}")
#             continue

#         # Select pressure level and slice latitude and longitude
#         temperature = temperature.sel(
#             pressure_level=1000.0,
#             latitude=slice(lat_max, lat_min),
#             longitude=slice(lon_min, lon_max)
#         )

#         # Check if data is available
#         if temperature.isnull().all():
#             print(f"No data available for event date: {event_date}")
#             continue

#         # Convert temperature from Kelvin to Celsius
#         temperature_celsius = temperature - 273.15

#         # Plot the 0°C isotherm
#         cs = ax.contour(
#             temperature_celsius.longitude,
#             temperature_celsius.latitude,
#             temperature_celsius,
#             levels=[0],
#             colors=color,
#             transform=ccrs.PlateCarree(),
#             linewidths=1
#         )

#     # Add labels for Washington D.C., Beijing, and Moscow
#     for city in cities:
#         ax.text(
#             city['lon'],
#             city['lat'],
#             city['name'],
#             transform=ccrs.PlateCarree(),
#             fontsize=12,
#             fontweight='bold',
#             ha='center',
#             va='center',
#             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2')
#         )

#     # Add the legend
#     ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.2)

# # Create the animation
# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=days_to_animate,
#     interval=1000,
#     blit=False
# )

# # To save the animation as a video file (optional)
# from matplotlib.animation import FFMpegWriter
# import matplotlib
# matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
# writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
# ani.save('temperature_animation_with_cities1000.mp4', writer=writer)

# # Display the animation
# plt.show()
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs  # For map projections
import cartopy.feature as cfeature  # For map features like coastlines
from matplotlib.lines import Line2D  # For creating custom legend entries

# Load the consolidated dataset containing temperature data
ds = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Since_2010_Temps_1000_500.nc")

# Function to adjust longitudes from 0-360 to -180 to 180 if necessary
def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        # Adjust longitude coordinates to range from -180 to 180
        ds = ds.assign_coords(
            longitude=(((ds.longitude + 180) % 360) - 180)
        )
        # Sort dataset by longitude for consistent plotting
        ds = ds.sortby('longitude')
    return ds

# Apply longitude adjustment to the dataset
ds = adjust_longitude(ds)

# Define event dates and labels for Sudden Stratospheric Warming (SSW) and Non-SSW events

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

# Define the latitude and longitude boundaries for the map
lat_min, lat_max = 20, 90  # From 20°N latitude up to the North Pole
lon_min, lon_max = -180, 180  # Full global longitude coverage

# Colors for plotting SSW and Non-SSW events
color_ssw = 'red'      # SSW events will be plotted in red
color_non_ssw = 'navy' # Non-SSW events will be plotted in navy

# Prepare the figure and axes with a North Polar Stereographic projection
fig = plt.figure(figsize=(14, 14))
ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())  # Set map boundaries

# Add map features to the axes
ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')  # Land areas
ax.add_feature(cfeature.COASTLINE)  # Coastlines
ax.add_feature(cfeature.BORDERS, linestyle=':')  # Country borders with dotted lines

# Add gridlines without labels
gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

# Create legend elements for SSW and Non-SSW events
legend_elements = [
    Line2D([0], [0], color=color_ssw, lw=2, label='SSW Events'),
    Line2D([0], [0], color=color_non_ssw, lw=2, label='Non-SSW Events')
]

# Combine SSW and Non-SSW events into one list with their dates, labels, and colors
all_events = []
for date, label in zip(ssw_event_dates, ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': color_ssw})

for date, label in zip(non_ssw_event_dates, non_ssw_event_labels):
    all_events.append({'date': date, 'label': label, 'color': color_non_ssw})

# Define the total number of days to animate after each event date
total_days = 14  # Animate for 14 days

# Create a list of days to animate (from day 0 to day 14)
days_to_animate = np.arange(0, total_days + 1)

# Define the cities to label on the map with their coordinates
cities = [
    {'name': 'Washington D.C.', 'lat': 38.9072, 'lon': -77.0369},
    {'name': 'Beijing', 'lat': 39.9042, 'lon': 116.4074},
    {'name': 'Moscow', 'lat': 55.7558, 'lon': 37.6173}
]

# Define the update function that will be called for each frame of the animation
def update(day_offset):
    ax.clear()  # Clear the axes for the new frame
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())  # Reset map boundaries

    # Re-add map features after clearing
    ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

    # Set the title for the current frame
    current_day = day_offset
    ax.set_title(f'0°C Isotherm at 1000 mb - Day {current_day} After Event Dates')

    # Loop through all events and plot the 0°C isotherm for each
    for event in all_events:
        event_date = event['date'] + np.timedelta64(day_offset, 'D')  # Calculate the date for this frame
        color = event['color']  # Get the color associated with the event

        # Try to select temperature data for the event date using the nearest available time
        try:
            temperature = ds['t'].sel(valid_time=event_date, method='nearest')
        except KeyError:
            print(f"No data available for event date: {event_date}")
            continue  # Skip to the next event if data is missing

        # Select the temperature data at 1000 mb pressure level within the latitude and longitude bounds
        temperature = temperature.sel(
            pressure_level=1000.0,
            latitude=slice(lat_max, lat_min),
            longitude=slice(lon_min, lon_max)
        )

        # Check if the selected data is valid (not all NaNs)
        if temperature.isnull().all():
            print(f"No data available for event date: {event_date}")
            continue  # Skip if data is invalid

        # Convert temperature from Kelvin to Celsius
        temperature_celsius = temperature - 273.15

        # Plot the 0°C isotherm (contour where temperature equals 0°C)
        cs = ax.contour(
            temperature_celsius.longitude,
            temperature_celsius.latitude,
            temperature_celsius,
            levels=[0],  # Contour level at 0°C
            colors=color,  # Use the event's color
            transform=ccrs.PlateCarree(),  # Coordinate reference system
            linewidths=1  # Line width of the contour
        )

    # Add labels for each city on the map
    for city in cities:
        ax.text(
            city['lon'], city['lat'], city['name'],
            transform=ccrs.PlateCarree(),
            fontsize=12,
            fontweight='bold',
            ha='center',  # Horizontal alignment
            va='center',  # Vertical alignment
            bbox=dict(
                facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'
            )  # Background box for the text
        )

    # Add the legend to the plot, positioned outside the main axes
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.2)

# Create the animation using the update function and specified frames
ani = animation.FuncAnimation(
    fig,
    update,            # Function to update each frame
    frames=days_to_animate,  # Frames to animate (days 0 to 14)
    interval=1000,     # Time between frames in milliseconds
    blit=False         # Disable blitting for simplicity
)

# Import FFMpegWriter to save the animation as a video file
from matplotlib.animation import FFMpegWriter
import matplotlib
# Specify the path to the ffmpeg executable (required for saving animations)
matplotlib.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\ffmpeg\ffmpeg-7.1-full_build\bin\ffmpeg.exe"
# Create an FFMpegWriter object with specified frames per second and bitrate
writer = FFMpegWriter(fps=1, metadata=dict(artist='Brayden Holler'), bitrate=1800)
# Save the animation to a video file using the writer
ani.save('temperature_animation_with_cities1000.mp4', writer=writer)

# Display the animation in the matplotlib viewer
plt.show()
