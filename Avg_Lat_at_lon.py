# Working Code just for SSW Events

# import rasterio 
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt

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

# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All")

# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)

# # Pre-subset the SSW data spatially and by pressure level to match the permafrost domain.
# ds_subset = ds_ssw['z'].sel(
#     pressure_level=[500.0, 1000.0],
#     latitude=slice(new_top, new_bottom),
#     longitude=slice(new_left, new_right)
# )

# #############################################
# # 3. Define SSW event dates (2010-present)   #
# #############################################

# ssw_event_dates = [
#     np.datetime64("2010-02-10", 'D'),
#     np.datetime64("2010-03-24", 'D'),
#     np.datetime64("2013-01-06", 'D'),
#     np.datetime64("2018-02-12", 'D'),
#     np.datetime64("2019-01-01", 'D'),
#     np.datetime64("2021-01-05", 'D'),
#     np.datetime64("2023-02-16", 'D'),
#     np.datetime64("2000-03-20", 'D'),
#     np.datetime64("2001-02-11", 'D'),
#     np.datetime64("2001-12-30", 'D'),
#     np.datetime64("2002-02-17", 'D'),
#     np.datetime64("2003-01-18", 'D'),
#     np.datetime64("2004-01-05", 'D'),
#     np.datetime64("2006-01-21", 'D'),
#     np.datetime64("2007-02-24", 'D'),
#     np.datetime64("2008-02-22", 'D'),
#     np.datetime64("2009-01-24", 'D'),
#     np.datetime64("1998-12-15", 'D'),
#     np.datetime64("1999-02-26", 'D'),
#     np.datetime64("1980-02-29", 'D'),
#     np.datetime64("1981-03-04", 'D'),
#     np.datetime64("1981-12-04", 'D'),
#     np.datetime64("1984-02-24", 'D'),
#     np.datetime64("1985-01-01", 'D'),
#     np.datetime64("1987-01-23", 'D'),
#     np.datetime64("1987-12-08", 'D'),
#     np.datetime64("1988-03-14", 'D'),
#     np.datetime64("1989-02-21", 'D'),
#     np.datetime64("1970-01-02", 'D'),
#     np.datetime64("1971-01-18", 'D'),
#     np.datetime64("1971-03-20", 'D'),
#     np.datetime64("1973-01-31", 'D'),
#     np.datetime64("1977-01-09", 'D'),
#     np.datetime64("1979-02-22", 'D'),
#     np.datetime64("1958-02-08", 'D'),
#     np.datetime64("1960-01-17", 'D'),
#     np.datetime64("1963-01-27", 'D'),
#     np.datetime64("1965-12-16", 'D'),
#     np.datetime64("1966-02-22", 'D'),
#     np.datetime64("1968-01-07", 'D'),
#     np.datetime64("1968-11-28", 'D'),
#     np.datetime64("1969-03-13", 'D')
# ]

# #############################################
# # 4. Set parameters for the analysis         #
# #############################################
# lon_of_interest = -105.0   # Specific longitude at which to extract the contour
# total_days = 14            # Post-event duration to consider
# num_days = total_days + 1  # Include day 0

# #############################################
# # 5. Compute average contour latitude         #
# #############################################
# # For each event date, we now select all 15 valid_time slices in one .sel call.
# results = []  # List to store results per event

# for event_date in ssw_event_dates:
#     # Build an array of valid times for the event (vectorized)
#     times = event_date + np.arange(num_days)
#     try:
#         # Select all times at once for this event.
#         data = ds_subset.sel(valid_time=times, method='nearest')
#     except Exception:
#         continue

#     # Convert geopotential to geopotential height (m)
#     geopot_height = data / 9.80665
#     hgt_500 = geopot_height.sel(pressure_level=500.0)
#     hgt_1000 = geopot_height.sel(pressure_level=1000.0)
    
#     # Compute thickness (500 hPa minus 1000 hPa)
#     thickness = hgt_500 - hgt_1000

#     # Interpolate thickness along the specified longitude for all valid times at once.
#     thickness_line = thickness.interp(longitude=lon_of_interest)

#     daily_crossings = []  # To hold the crossing latitude for each day

#     # Loop over each time slice (there are only 15, so this loop is inexpensive)
#     for valid_time in thickness_line.valid_time.values:
#         profile = thickness_line.sel(valid_time=valid_time).values  # 1D array along latitude
#         lats = thickness_line.latitude.values
#         diff = profile - 5400.0

#         # Identify indices where the sign changes (zero-crossings)
#         indices = np.where(diff[:-1] * diff[1:] < 0)[0]
#         crossing_lats = []

#         # For each zero crossing, do linear interpolation
#         for i in indices:
#             # Avoid division by zero and perform linear interpolation
#             d1, d2 = diff[i], diff[i+1]
#             if d2 - d1 == 0:
#                 continue
#             frac = -d1 / (d2 - d1)
#             lat1, lat2 = lats[i], lats[i+1]
#             crossing_lats.append(lat1 + frac * (lat2 - lat1))
        
#         if crossing_lats:
#             daily_crossings.append(np.mean(crossing_lats))
    
#     if daily_crossings:
#         overall_avg = np.mean(daily_crossings)
#         results.append({'date': event_date, 'avg_lat': overall_avg})

# #############################################
# # 6. Plot the average crossing latitudes     #
# #############################################
# fig, ax = plt.subplots(figsize=(6, 8))

# # Draw a vertical reference line at the specified longitude.
# ax.axvline(x=lon_of_interest, color='gray', linestyle='--')

# # Plot each event's average crossing latitude and annotate with the event date.
# for res in results:
#     ax.plot(lon_of_interest, res['avg_lat'], 'o', markersize=8)
#     ax.text(lon_of_interest + 0.1, res['avg_lat'], str(res['date'].astype('M8[D]')),
#             fontsize=9, verticalalignment='center')

# # Limit the x-axis to a narrow window around the longitude of interest.
# ax.set_xlim(lon_of_interest - 2, lon_of_interest + 2)

# # Adjust y-axis limits based on the computed latitudes.
# if results:
#     all_lats = [res['avg_lat'] for res in results]
#     margin = 2
#     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)

# ax.set_xlabel("Longitude (°)")
# ax.set_ylabel("Latitude (°)")
# ax.set_title(f"Average 5400 m Thickness Contour Latitude\nat lon {lon_of_interest} (averaged over {total_days} days post-SSW)")
# plt.tight_layout()
# plt.show()

# Very slow, did not finish running after 10 minutes
# The code is inefficient because it loops over each event date and then loops over each time slice.
# import rasterio 
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt

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

# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All")

# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)

# # Pre-subset the SSW data spatially and by pressure level to match the permafrost domain.
# ds_subset = ds_ssw['z'].sel(
#     pressure_level=[500.0, 1000.0],
#     latitude=slice(new_top, new_bottom),
#     longitude=slice(new_left, new_right)
# )

# #############################################
# # 3. Define SSW event and non-event dates   #
# #############################################

# # SSW event dates (from various periods)
# ssw_event_dates = [
#     np.datetime64("2010-02-10", 'D'),
#     np.datetime64("2010-03-24", 'D'),
#     np.datetime64("2013-01-06", 'D'),
#     np.datetime64("2018-02-12", 'D'),
#     np.datetime64("2019-01-01", 'D'),
#     np.datetime64("2021-01-05", 'D'),
#     np.datetime64("2023-02-16", 'D'),
#     np.datetime64("2000-03-20", 'D'),
#     np.datetime64("2001-02-11", 'D'),
#     np.datetime64("2001-12-30", 'D'),
#     np.datetime64("2002-02-17", 'D'),
#     np.datetime64("2003-01-18", 'D'),
#     np.datetime64("2004-01-05", 'D'),
#     np.datetime64("2006-01-21", 'D'),
#     np.datetime64("2007-02-24", 'D'),
#     np.datetime64("2008-02-22", 'D'),
#     np.datetime64("2009-01-24", 'D'),
#     np.datetime64("1998-12-15", 'D'),
#     np.datetime64("1999-02-26", 'D'),
#     np.datetime64("1980-02-29", 'D'),
#     np.datetime64("1981-03-04", 'D'),
#     np.datetime64("1981-12-04", 'D'),
#     np.datetime64("1984-02-24", 'D'),
#     np.datetime64("1985-01-01", 'D'),
#     np.datetime64("1987-01-23", 'D'),
#     np.datetime64("1987-12-08", 'D'),
#     np.datetime64("1988-03-14", 'D'),
#     np.datetime64("1989-02-21", 'D'),
#     np.datetime64("1970-01-02", 'D'),
#     np.datetime64("1971-01-18", 'D'),
#     np.datetime64("1971-03-20", 'D'),
#     np.datetime64("1973-01-31", 'D'),
#     np.datetime64("1977-01-09", 'D'),
#     np.datetime64("1979-02-22", 'D'),
#     np.datetime64("1958-02-08", 'D'),
#     np.datetime64("1960-01-17", 'D'),
#     np.datetime64("1963-01-27", 'D'),
#     np.datetime64("1965-12-16", 'D'),
#     np.datetime64("1966-02-22", 'D'),
#     np.datetime64("1968-01-07", 'D'),
#     np.datetime64("1968-11-28", 'D'),
#     np.datetime64("1969-03-13", 'D')
# ]

# # Non-event dates: for selected years, we choose January 30 as the reference date.
# non_event_dates = [
#     np.datetime64("1959-01-30", 'D'),
#     np.datetime64("1961-01-30", 'D'),
#     np.datetime64("1962-01-30", 'D'),
#     np.datetime64("1964-01-30", 'D'),
#     np.datetime64("1967-01-30", 'D'),
#     np.datetime64("1972-01-30", 'D'),
#     np.datetime64("1974-01-30", 'D'),
#     np.datetime64("1975-01-30", 'D'),
#     np.datetime64("1976-01-30", 'D'),
#     np.datetime64("1978-01-30", 'D'),
#     np.datetime64("1982-01-30", 'D'),
#     np.datetime64("1983-01-30", 'D'),
#     np.datetime64("1986-01-30", 'D'),
#     np.datetime64("1990-01-30", 'D'),
#     np.datetime64("1991-01-30", 'D'),
#     np.datetime64("1992-01-30", 'D'),
#     np.datetime64("1993-01-30", 'D'),
#     np.datetime64("1994-01-30", 'D'),
#     np.datetime64("1995-01-30", 'D'),
#     np.datetime64("1996-01-30", 'D'),
#     np.datetime64("1997-01-30", 'D'),
#     np.datetime64("2005-01-30", 'D'),
#     np.datetime64("2011-01-30", 'D'),
#     np.datetime64("2012-01-30", 'D'),
#     np.datetime64("2014-01-30", 'D'),
#     np.datetime64("2015-01-30", 'D'),
#     np.datetime64("2016-01-30", 'D'),
#     np.datetime64("2017-01-30", 'D'),
#     np.datetime64("2020-01-30", 'D'),
#     np.datetime64("2022-01-30", 'D')
# ]

# #############################################
# # 4. Set parameters for the analysis         #
# #############################################
# lon_of_interest = -105.0   # Specific longitude at which to extract the contour
# total_days = 14            # Post-date duration to consider
# num_days = total_days + 1  # Include day 0

# #############################################
# # 5. Define a helper function to compute     #
# #    average contour latitude over num_days  #
# #############################################

# def compute_avg_contour(start_date):
#     """
#     Given a starting date, compute the average latitude where the 5400 m thickness
#     contour occurs over a num_days period.
#     """
#     # Build an array of valid times for the period
#     times = start_date + np.arange(num_days)
#     try:
#         data = ds_subset.sel(valid_time=times, method='nearest')
#     except Exception:
#         return None
#     # Convert geopotential (m^2/s^2) to geopotential height (m)
#     geopot_height = data / 9.80665
#     hgt_500 = geopot_height.sel(pressure_level=500.0)
#     hgt_1000 = geopot_height.sel(pressure_level=1000.0)
#     thickness = hgt_500 - hgt_1000

#     # Interpolate thickness along the specified longitude for all times
#     thickness_line = thickness.interp(longitude=lon_of_interest)
#     daily_crossings = []

#     # Loop over each time slice (only 15 slices)
#     for valid_time in thickness_line.valid_time.values:
#         profile = thickness_line.sel(valid_time=valid_time).values  # 1D array along latitude
#         lats = thickness_line.latitude.values
#         diff = profile - 5400.0

#         # Find indices where diff changes sign
#         indices = np.where(diff[:-1] * diff[1:] < 0)[0]
#         crossing_lats = []
#         for i in indices:
#             d1, d2 = diff[i], diff[i+1]
#             if d2 - d1 == 0:
#                 continue
#             frac = -d1 / (d2 - d1)
#             lat1, lat2 = lats[i], lats[i+1]
#             crossing_lats.append(lat1 + frac * (lat2 - lat1))
#         if crossing_lats:
#             daily_crossings.append(np.mean(crossing_lats))
#     if daily_crossings:
#         return np.mean(daily_crossings)
#     else:
#         return None

# #############################################
# # 6. Compute average contour latitudes for   #
# #    both event and non-event dates          #
# #############################################
# results = []  # List to store results as dicts with date, avg_lat, and type

# # Process event dates (mark type as 'event')
# for event_date in ssw_event_dates:
#     avg_lat = compute_avg_contour(event_date)
#     if avg_lat is not None:
#         results.append({'date': event_date, 'avg_lat': avg_lat, 'type': 'event'})

# # Process non-event dates (mark type as 'non-event')
# for non_event_date in non_event_dates:
#     avg_lat = compute_avg_contour(non_event_date)
#     if avg_lat is not None:
#         results.append({'date': non_event_date, 'avg_lat': avg_lat, 'type': 'non-event'})

# #############################################
# # 7. Plot the average crossing latitudes     #
# #############################################
# fig, ax = plt.subplots(figsize=(6, 8))

# # Draw a vertical reference line at the specified longitude.
# ax.axvline(x=lon_of_interest, color='gray', linestyle='--')

# # Plot each result using a different color based on type.
# for res in results:
#     if res['type'] == 'event':
#         marker_color = 'blue'
#     else:
#         marker_color = 'red'
#     ax.plot(lon_of_interest, res['avg_lat'], 'o', markersize=8, color=marker_color)
#     ax.text(lon_of_interest + 0.1, res['avg_lat'], str(res['date'].astype('M8[D]')),
#             fontsize=9, verticalalignment='center', color=marker_color)

# # Limit the x-axis to a narrow window around the longitude of interest.
# ax.set_xlim(lon_of_interest - 2, lon_of_interest + 2)

# # Adjust y-axis limits based on the computed latitudes.
# if results:
#     all_lats = [res['avg_lat'] for res in results]
#     margin = 2
#     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)

# ax.set_xlabel("Longitude (°)")
# ax.set_ylabel("Latitude (°)")
# ax.set_title(f"Average 5400 m Thickness Contour Latitude\nat lon {lon_of_interest} (averaged over {total_days} days)")
# plt.tight_layout()

# # Add legend manually
# blue_marker = plt.Line2D([0], [0], marker='o', color='w', label='Event Years',
#                            markerfacecolor='blue', markersize=8)
# red_marker = plt.Line2D([0], [0], marker='o', color='w', label='Non-Event Years',
#                            markerfacecolor='red', markersize=8)
# ax.legend(handles=[blue_marker, red_marker], loc='best')

# plt.show()

# #!/usr/bin/env python3
# import rasterio 
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt

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
# # 2. Load geopotential height dataset (SSW) #
# #############################################
# ds_ssw = xr.open_dataset(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All")

# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
#         ds = ds.sortby('longitude')
#     return ds

# ds_ssw = adjust_longitude(ds_ssw)

# # Pre-subset spatially and by pressure level.
# ds_subset = ds_ssw['z'].sel(
#     pressure_level=[500.0, 1000.0],
#     latitude=slice(new_top, new_bottom),
#     longitude=slice(new_left, new_right)
# )

# #############################################
# # 3. Define SSW event and non-event dates   #
# #############################################
# # SSW event dates (from various periods)
# ssw_event_dates = [
#     np.datetime64("2010-02-10", 'D'),
#     np.datetime64("2010-03-24", 'D'),
#     np.datetime64("2013-01-06", 'D'),
#     np.datetime64("2018-02-12", 'D'),
#     np.datetime64("2019-01-01", 'D'),
#     np.datetime64("2021-01-05", 'D'),
#     np.datetime64("2023-02-16", 'D'),
#     np.datetime64("2000-03-20", 'D'),
#     np.datetime64("2001-02-11", 'D'),
#     np.datetime64("2001-12-30", 'D'),
#     np.datetime64("2002-02-17", 'D'),
#     np.datetime64("2003-01-18", 'D'),
#     np.datetime64("2004-01-05", 'D'),
#     np.datetime64("2006-01-21", 'D'),
#     np.datetime64("2007-02-24", 'D'),
#     np.datetime64("2008-02-22", 'D'),
#     np.datetime64("2009-01-24", 'D'),
#     np.datetime64("1998-12-15", 'D'),
#     np.datetime64("1999-02-26", 'D'),
#     np.datetime64("1980-02-29", 'D'),
#     np.datetime64("1981-03-04", 'D'),
#     np.datetime64("1981-12-04", 'D'),
#     np.datetime64("1984-02-24", 'D'),
#     np.datetime64("1985-01-01", 'D'),
#     np.datetime64("1987-01-23", 'D'),
#     np.datetime64("1987-12-08", 'D'),
#     np.datetime64("1988-03-14", 'D'),
#     np.datetime64("1989-02-21", 'D'),
#     np.datetime64("1970-01-02", 'D'),
#     np.datetime64("1971-01-18", 'D'),
#     np.datetime64("1971-03-20", 'D'),
#     np.datetime64("1973-01-31", 'D'),
#     np.datetime64("1977-01-09", 'D'),
#     np.datetime64("1979-02-22", 'D'),
#     np.datetime64("1958-02-08", 'D'),
#     np.datetime64("1960-01-17", 'D'),
#     np.datetime64("1963-01-27", 'D'),
#     np.datetime64("1965-12-16", 'D'),
#     np.datetime64("1966-02-22", 'D'),
#     np.datetime64("1968-01-07", 'D'),
#     np.datetime64("1968-11-28", 'D'),
#     np.datetime64("1969-03-13", 'D')
# ]

# # Non-event dates: January 30 for selected years.
# non_event_dates = [
#     np.datetime64("1959-01-30", 'D'),
#     np.datetime64("1961-01-30", 'D'),
#     np.datetime64("1962-01-30", 'D'),
#     np.datetime64("1964-01-30", 'D'),
#     np.datetime64("1967-01-30", 'D'),
#     np.datetime64("1972-01-30", 'D'),
#     np.datetime64("1974-01-30", 'D'),
#     np.datetime64("1975-01-30", 'D'),
#     np.datetime64("1976-01-30", 'D'),
#     np.datetime64("1978-01-30", 'D'),
#     np.datetime64("1982-01-30", 'D'),
#     np.datetime64("1983-01-30", 'D'),
#     np.datetime64("1986-01-30", 'D'),
#     np.datetime64("1990-01-30", 'D'),
#     np.datetime64("1991-01-30", 'D'),
#     np.datetime64("1992-01-30", 'D'),
#     np.datetime64("1993-01-30", 'D'),
#     np.datetime64("1994-01-30", 'D'),
#     np.datetime64("1995-01-30", 'D'),
#     np.datetime64("1996-01-30", 'D'),
#     np.datetime64("1997-01-30", 'D'),
#     np.datetime64("2005-01-30", 'D'),
#     np.datetime64("2011-01-30", 'D'),
#     np.datetime64("2012-01-30", 'D'),
#     np.datetime64("2014-01-30", 'D'),
#     np.datetime64("2015-01-30", 'D'),
#     np.datetime64("2016-01-30", 'D'),
#     np.datetime64("2017-01-30", 'D'),
#     np.datetime64("2020-01-30", 'D'),
#     np.datetime64("2022-01-30", 'D')
# ]

# #############################################
# # 4. Set parameters for the analysis         #
# #############################################
# lon_of_interest = 35.0   # Specific longitude at which to extract the contour
# total_days = 14            # Post-date duration to consider
# num_days = total_days + 1  # Include day 0

# #############################################
# # 5. Define a helper function to compute     #
# #    average contour latitude over num_days  #
# #############################################
# def compute_avg_contour(start_date):
#     """
#     For a given start_date, compute the average latitude where the 5400 m thickness
#     contour occurs over a period of num_days.
#     """
#     times = start_date + np.arange(num_days)
#     try:
#         # Select only the 15 required time slices (this loads a small subset)
#         data = ds_subset.sel(valid_time=times, method='nearest')
#     except Exception:
#         return None
    
#     # Compute geopotential height and thickness on the fly.
#     geopot_height = data / 9.80665
#     hgt_500 = geopot_height.sel(pressure_level=500.0)
#     hgt_1000 = geopot_height.sel(pressure_level=1000.0)
#     thickness = hgt_500 - hgt_1000

#     # Interpolate thickness along the desired longitude.
#     thickness_line = thickness.interp(longitude=lon_of_interest)
#     daily_crossings = []
    
#     # Loop over the 15 time slices (this loop is inexpensive)
#     for valid_time in thickness_line.valid_time.values:
#         profile = thickness_line.sel(valid_time=valid_time).values  # 1D array along latitude
#         lats = thickness_line.latitude.values
#         diff = profile - 5400.0
#         indices = np.where(diff[:-1] * diff[1:] < 0)[0]
#         crossing_lats = []
#         for i in indices:
#             d1, d2 = diff[i], diff[i+1]
#             if d2 - d1 == 0:
#                 continue
#             frac = -d1 / (d2 - d1)
#             lat1, lat2 = lats[i], lats[i+1]
#             crossing_lats.append(lat1 + frac * (lat2 - lat1))
#         if crossing_lats:
#             daily_crossings.append(np.mean(crossing_lats))
#     return np.mean(daily_crossings) if daily_crossings else None

# #############################################
# # 6. Compute average contour latitudes for   #
# #    both event and non-event dates          #
# #############################################
# results = []  # Will hold dicts with date, avg_lat, and type

# # Process event dates (type 'event')
# for event_date in ssw_event_dates:
#     avg_lat = compute_avg_contour(event_date)
#     if avg_lat is not None:
#         results.append({'date': event_date, 'avg_lat': avg_lat, 'type': 'event'})

# # Process non-event dates (type 'non-event')
# for non_event_date in non_event_dates:
#     avg_lat = compute_avg_contour(non_event_date)
#     if avg_lat is not None:
#         results.append({'date': non_event_date, 'avg_lat': avg_lat, 'type': 'non-event'})

# # #############################################
# # # 7. Plot the average crossing latitudes     #
# # #############################################
# # fig, ax = plt.subplots(figsize=(6, 8))
# # ax.axvline(x=lon_of_interest, color='gray', linestyle='--')

# # for res in results:
# #     marker_color = 'blue' if res['type'] == 'event' else 'red'
# #     ax.plot(lon_of_interest, res['avg_lat'], 'o', markersize=8, color=marker_color)
# #     ax.text(lon_of_interest + 0.1, res['avg_lat'], str(res['date'].astype('M8[D]')),
# #             fontsize=9, verticalalignment='center', color=marker_color)

# # ax.set_xlim(lon_of_interest - 2, lon_of_interest + 2)
# # if results:
# #     all_lats = [res['avg_lat'] for res in results]
# #     margin = 2
# #     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
# # ax.set_xlabel("Longitude (°)")
# # ax.set_ylabel("Latitude (°)")
# # ax.set_title(f"Average 5400 m Thickness Contour Latitude\nat lon {lon_of_interest} (averaged over {total_days} days)")
# # plt.tight_layout()

# # # Add legend manually
# # blue_marker = plt.Line2D([0], [0], marker='o', color='w', label='Event Years',
# #                            markerfacecolor='blue', markersize=8)
# # red_marker = plt.Line2D([0], [0], marker='o', color='w', label='Non-Event Years',
# #                            markerfacecolor='red', markersize=8)
# # ax.legend(handles=[blue_marker, red_marker], loc='best')
# # plt.show()

# #############################################
# # 7. Plot the maximum contour latitudes as  #
# #    a strip plot to better show density    #
# #############################################
# fig, ax = plt.subplots(figsize=(6, 8))
# jitter_strength = 0.1  # controls horizontal spread

# # Instead of a fixed x value, assign group positions:
# for res in results:
#     if res['type'] == 'event':
#         x = 1 + np.random.uniform(-jitter_strength, jitter_strength)
#         color = 'blue'
#     else:
#         x = 2 + np.random.uniform(-jitter_strength, jitter_strength)
#         color = 'red'
#     ax.plot(x, res['avg_lat'], 'o', markersize=8, color=color)
#     ax.text(x + 0.05, res['avg_lat'], str(res['date'].astype('M8[D]')),
#             fontsize=9, verticalalignment='center', color=color)

# # Set x-axis to show group labels
# ax.set_xlim(0.5, 2.5)
# ax.set_xticks([1, 2])
# ax.set_xticklabels(['Event Years', 'Non-Event Years'])

# # Set y-axis limits based on latitude range
# if results:
#     all_lats = [res['avg_lat'] for res in results]
#     margin = 2
#     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
# ax.set_ylabel("Latitude (°)")
# ax.set_title(f"Density of Average 5400 m Thickness Contour Latitude\n(over {total_days} days)")
# plt.tight_layout()
# plt.show()

#!/usr/bin/env python3
import json
import numpy as np
import matplotlib.pyplot as plt

#############################################
# 1. Load precomputed raw contour JSON data  #
#############################################
# (This JSON file was produced by your uploaded script.)
json_file = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\all_5400m_contours_15day_raw_all.json"
with open(json_file, "r") as f:
    raw_data = json.load(f)

#############################################
# 2. Define SSW event and non-event dates   #
#############################################
ssw_event_dates = [
    np.datetime64("2010-02-10", 'D'),
    np.datetime64("2010-03-24", 'D'),
    np.datetime64("2013-01-06", 'D'),
    np.datetime64("2018-02-12", 'D'),
    np.datetime64("2019-01-01", 'D'),
    np.datetime64("2021-01-05", 'D'),
    np.datetime64("2023-02-16", 'D'),
    np.datetime64("2000-03-20", 'D'),
    np.datetime64("2001-02-11", 'D'),
    np.datetime64("2001-12-30", 'D'),
    np.datetime64("2002-02-17", 'D'),
    np.datetime64("2003-01-18", 'D'),
    np.datetime64("2004-01-05", 'D'),
    np.datetime64("2006-01-21", 'D'),
    np.datetime64("2007-02-24", 'D'),
    np.datetime64("2008-02-22", 'D'),
    np.datetime64("2009-01-24", 'D'),
    np.datetime64("1998-12-15", 'D'),
    np.datetime64("1999-02-26", 'D'),
    np.datetime64("1980-02-29", 'D'),
    np.datetime64("1981-03-04", 'D'),
    np.datetime64("1981-12-04", 'D'),
    np.datetime64("1984-02-24", 'D'),
    np.datetime64("1985-01-01", 'D'),
    np.datetime64("1987-01-23", 'D'),
    np.datetime64("1987-12-08", 'D'),
    np.datetime64("1988-03-14", 'D'),
    np.datetime64("1989-02-21", 'D'),
    np.datetime64("1970-01-02", 'D'),
    np.datetime64("1971-01-18", 'D'),
    np.datetime64("1971-03-20", 'D'),
    np.datetime64("1973-01-31", 'D'),
    np.datetime64("1977-01-09", 'D'),
    np.datetime64("1979-02-22", 'D'),
    np.datetime64("1958-02-08", 'D'),
    np.datetime64("1960-01-17", 'D'),
    np.datetime64("1963-01-27", 'D'),
    np.datetime64("1965-12-16", 'D'),
    np.datetime64("1966-02-22", 'D'),
    np.datetime64("1968-01-07", 'D'),
    np.datetime64("1968-11-28", 'D'),
    np.datetime64("1969-03-13", 'D')
]

non_event_dates = [
    np.datetime64("1959-01-30", 'D'),
    np.datetime64("1961-01-30", 'D'),
    np.datetime64("1962-01-30", 'D'),
    np.datetime64("1964-01-30", 'D'),
    np.datetime64("1967-01-30", 'D'),
    np.datetime64("1972-01-30", 'D'),
    np.datetime64("1974-01-30", 'D'),
    np.datetime64("1975-01-30", 'D'),
    np.datetime64("1976-01-30", 'D'),
    np.datetime64("1978-01-30", 'D'),
    np.datetime64("1982-01-30", 'D'),
    np.datetime64("1983-01-30", 'D'),
    np.datetime64("1986-01-30", 'D'),
    np.datetime64("1990-01-30", 'D'),
    np.datetime64("1991-01-30", 'D'),
    np.datetime64("1992-01-30", 'D'),
    np.datetime64("1993-01-30", 'D'),
    np.datetime64("1994-01-30", 'D'),
    np.datetime64("1995-01-30", 'D'),
    np.datetime64("1996-01-30", 'D'),
    np.datetime64("1997-01-30", 'D'),
    np.datetime64("2005-01-30", 'D'),
    np.datetime64("2011-01-30", 'D'),
    np.datetime64("2012-01-30", 'D'),
    np.datetime64("2014-01-30", 'D'),
    np.datetime64("2015-01-30", 'D'),
    np.datetime64("2016-01-30", 'D'),
    np.datetime64("2017-01-30", 'D'),
    np.datetime64("2020-01-30", 'D'),
    np.datetime64("2022-01-30", 'D')
]

#############################################
# 3. Set parameters for the analysis         #
#############################################
lon_of_interest = -180.0   # Specific longitude at which to extract the contour
total_days = 14            # 15-day period (day 0 through day 14)
num_days = total_days + 1

#############################################
# 4. Define helper function to compute the   #
#    average contour crossing latitude from  #
#    JSON data (using precomputed 'crossings') #
#############################################
def compute_avg_contour_from_json(start_date, raw_data, lon_of_interest):
    """
    For a given start date (as a numpy.datetime64), find the matching JSON entry,
    then for each day in the period, extract the list of crossing latitudes (from
    the precomputed 'crossings' dictionary) at the target longitude and compute
    the daily average. Finally, return the overall average across the period.
    """
    start_date_str = str(start_date)
    matching = [entry for entry in raw_data if entry["start_date"] == start_date_str]
    if not matching:
        print(f"No data found for start date {start_date_str}")
        return None

    period = matching[0]["period"]
    daily_avgs = []
    for day in period:
        # 'crossings' is a dictionary with keys as target longitudes (strings)
        crossings = day["crossings"].get(str(lon_of_interest), [])
        if crossings:
            daily_avgs.append(np.mean(crossings))
    return np.mean(daily_avgs) if daily_avgs else None

#############################################
# 5. Compute average contour latitudes for     #
#    both event and non-event dates            #
#############################################
results = []  # List to store results as dicts with date, avg_lat, and type

# Process event dates
for event_date in ssw_event_dates:
    avg_lat = compute_avg_contour_from_json(event_date, raw_data, lon_of_interest)
    if avg_lat is not None:
        results.append({'date': event_date, 'avg_lat': avg_lat, 'type': 'event'})

# Process non-event dates
for non_event_date in non_event_dates:
    avg_lat = compute_avg_contour_from_json(non_event_date, raw_data, lon_of_interest)
    if avg_lat is not None:
        results.append({'date': non_event_date, 'avg_lat': avg_lat, 'type': 'non-event'})

#############################################
# 6. Plot the average crossing latitudes       #
#############################################
fig, ax = plt.subplots(figsize=(6, 8))

# Draw a vertical reference line at the target longitude.
ax.axvline(x=lon_of_interest, color='gray', linestyle='--')

# Plot each result with a marker color based on event type.
for res in results:
    marker_color = 'blue' if res['type'] == 'event' else 'red'
    ax.plot(lon_of_interest, res['avg_lat'], 'o', markersize=8, color=marker_color)
    ax.text(lon_of_interest + 0.1, res['avg_lat'], str(res['date'].astype('M8[D]')),
            fontsize=9, verticalalignment='center', color=marker_color)

# Set a narrow x-axis window around the target longitude.
ax.set_xlim(lon_of_interest - 2, lon_of_interest + 2)

# Adjust y-axis limits based on computed latitudes.
if results:
    all_lats = [res['avg_lat'] for res in results]
    margin = 2
    ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)

ax.set_xlabel("Longitude (°)")
ax.set_ylabel("Latitude (°)")
ax.set_title(f"Average 5400 m Thickness Contour Latitude\nat lon {lon_of_interest} (averaged over {total_days} days)")
plt.tight_layout()

# Add a manual legend.
blue_marker = plt.Line2D([0], [0], marker='o', color='w', label='Event Years',
                           markerfacecolor='blue', markersize=8)
red_marker = plt.Line2D([0], [0], marker='o', color='w', label='Non-Event Years',
                           markerfacecolor='red', markersize=8)
ax.legend(handles=[blue_marker, red_marker], loc='best')

plt.show()
