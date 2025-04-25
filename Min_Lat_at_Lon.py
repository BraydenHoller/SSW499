#!/usr/bin/env python3
import rasterio 
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#############################################
# 1. Load and subset the permafrost raster  #
#############################################
byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\PermaFrost\llipa.byte'

with rasterio.open(byte_file_path) as perma_ds:
    perma_data = perma_ds.read(1)
    pixel_size_y = abs(perma_ds.transform.e)
    nrows = perma_data.shape[0]
    rows = np.arange(nrows)
    center_lats = perma_ds.bounds.top - (rows + 0.5) * pixel_size_y
    valid_rows = np.where(center_lats > 25)[0]
    if valid_rows.size == 0:
        raise ValueError("No rows found with latitude above 25°N.")
    last_valid_row = valid_rows[-1]
    perma_subset = perma_data[:last_valid_row + 1, :]
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

ds_subset = ds_ssw['z'].sel(
    pressure_level=[500.0, 1000.0],
    latitude=slice(new_top, new_bottom),
    longitude=slice(new_left, new_right)
)

#############################################
# 3. Define SSW event and non-event dates   #
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
# 4. Set parameters for the analysis         #
#############################################
lon_of_interest = -105.0   # Specific longitude at which to extract the contour
total_days = 14            # Post-date duration to consider
num_days = total_days + 1  # Include day 0

#############################################
# 5. Helper function to compute minimum      #
#    contour crossing latitude over num_days #
#############################################
def compute_min_contour(start_date):
    """
    For a given start_date, compute the minimum latitude where the 5400 m thickness
    contour occurs over a period of num_days.
    """
    times = start_date + np.arange(num_days)
    try:
        data = ds_subset.sel(valid_time=times, method='nearest')
    except Exception:
        return None
    
    # Compute geopotential height and thickness on the fly.
    geopot_height = data / 9.80665
    hgt_500 = geopot_height.sel(pressure_level=500.0)
    hgt_1000 = geopot_height.sel(pressure_level=1000.0)
    thickness = hgt_500 - hgt_1000

    # Interpolate thickness along the desired longitude.
    thickness_line = thickness.interp(longitude=lon_of_interest)
    daily_crossings = []
    
    # Loop over the 15 time slices.
    for valid_time in thickness_line.valid_time.values:
        profile = thickness_line.sel(valid_time=valid_time).values  # 1D array along latitude
        lats = thickness_line.latitude.values
        diff = profile - 5400.0
        indices = np.where(diff[:-1] * diff[1:] < 0)[0]
        crossing_lats = []
        for i in indices:
            d1, d2 = diff[i], diff[i+1]
            if d2 - d1 == 0:
                continue
            frac = -d1 / (d2 - d1)
            lat1, lat2 = lats[i], lats[i+1]
            crossing_lats.append(lat1 + frac * (lat2 - lat1))
        if crossing_lats:
            # Use the daily average if multiple crossings exist.
            daily_crossings.append(np.mean(crossing_lats))
    
    # Return the minimum daily crossing latitude over the period.
    return np.min(daily_crossings) if daily_crossings else None

#############################################
# 6. Compute minimum contour latitudes for    #
#    both event and non-event dates           #
#############################################
results = []  # Will hold dicts with date, min_lat, and type

# Process event dates (type 'event')
for event_date in ssw_event_dates:
    min_lat = compute_min_contour(event_date)
    if min_lat is not None:
        results.append({'date': event_date, 'min_lat': min_lat, 'type': 'event'})

# Process non-event dates (type 'non-event')
for non_event_date in non_event_dates:
    min_lat = compute_min_contour(non_event_date)
    if min_lat is not None:
        results.append({'date': non_event_date, 'min_lat': min_lat, 'type': 'non-event'})

# #############################################
# # 7. Plot the minimum contour latitudes       #
# #############################################
# fig, ax = plt.subplots(figsize=(6, 8))
# ax.axvline(x=lon_of_interest, color='gray', linestyle='--')

# for res in results:
#     marker_color = 'blue' if res['type'] == 'event' else 'red'
#     ax.plot(lon_of_interest, res['min_lat'], 'o', markersize=8, color=marker_color)
#     ax.text(lon_of_interest + 0.1, res['min_lat'], str(res['date'].astype('M8[D]')),
#             fontsize=9, verticalalignment='center', color=marker_color)

# ax.set_xlim(lon_of_interest - 2, lon_of_interest + 2)
# if results:
#     all_lats = [res['min_lat'] for res in results]
#     margin = 2
#     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
# ax.set_xlabel("Longitude (°)")
# ax.set_ylabel("Latitude (°)")
# ax.set_title(f"Minimum 5400 m Thickness Contour Latitude\nat lon {lon_of_interest} (over {total_days} days)")
# plt.tight_layout()

# blue_marker = plt.Line2D([0], [0], marker='o', color='w', label='Event Years',
#                            markerfacecolor='blue', markersize=8)
# red_marker = plt.Line2D([0], [0], marker='o', color='w', label='Non-Event Years',
#                            markerfacecolor='red', markersize=8)
# ax.legend(handles=[blue_marker, red_marker], loc='best')
# plt.show()

#############################################
# 7. Plot the maximum contour latitudes as  #
#    a strip plot to better show density    #
#############################################
fig, ax = plt.subplots(figsize=(6, 8))
jitter_strength = 0.1  # controls horizontal spread

# Instead of a fixed x value, assign group positions:
for res in results:
    if res['type'] == 'event':
        x = 1 + np.random.uniform(-jitter_strength, jitter_strength)
        color = 'blue'
    else:
        x = 2 + np.random.uniform(-jitter_strength, jitter_strength)
        color = 'red'
    ax.plot(x, res['min_lat'], 'o', markersize=8, color=color)
    ax.text(x + 0.05, res['min_lat'], str(res['date'].astype('M8[D]')),
            fontsize=9, verticalalignment='center', color=color)

# Set x-axis to show group labels
ax.set_xlim(0.5, 2.5)
ax.set_xticks([1, 2])
ax.set_xticklabels(['Event Years', 'Non-Event Years'])

# Set y-axis limits based on latitude range
if results:
    all_lats = [res['min_lat'] for res in results]
    margin = 2
    ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
ax.set_ylabel("Latitude (°)")
ax.set_title(f"Density of Minimum 5400 m Thickness Contour Latitude\n(over {total_days} days)")
plt.tight_layout()
plt.show()
