# !/usr/bin/env python3
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
# 5. Helper function to compute maximum      #
#    contour crossing latitude over num_days #
#############################################
def compute_max_contour(start_date):
    """
    For a given start_date, compute the maximum latitude where the 5400 m thickness
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
            # For the day, we take the average if multiple crossings exist.
            daily_crossings.append(np.mean(crossing_lats))
    
    # Return the maximum daily crossing latitude over the period.
    return np.max(daily_crossings) if daily_crossings else None

#############################################
# 6. Compute maximum contour latitudes for    #
#    both event and non-event dates           #
#############################################
results = []  # Will hold dicts with date, max_lat, and type

# Process event dates (type 'event')
for event_date in ssw_event_dates:
    max_lat = compute_max_contour(event_date)
    if max_lat is not None:
        results.append({'date': event_date, 'max_lat': max_lat, 'type': 'event'})

# Process non-event dates (type 'non-event')
for non_event_date in non_event_dates:
    max_lat = compute_max_contour(non_event_date)
    if max_lat is not None:
        results.append({'date': non_event_date, 'max_lat': max_lat, 'type': 'non-event'})

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
    ax.plot(x, res['max_lat'], 'o', markersize=8, color=color)
    ax.text(x + 0.05, res['max_lat'], str(res['date'].astype('M8[D]')),
            fontsize=9, verticalalignment='center', color=color)

# Set x-axis to show group labels
ax.set_xlim(0.5, 2.5)
ax.set_xticks([1, 2])
ax.set_xticklabels(['Event Years', 'Non-Event Years'])

# Set y-axis limits based on latitude range
if results:
    all_lats = [res['max_lat'] for res in results]
    margin = 2
    ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
ax.set_ylabel("Latitude (°)")
ax.set_title(f"Density of Maximum 5400 m Thickness Contour Latitude\n(over {total_days} days)")
plt.tight_layout()
plt.show()

# import json
# import numpy as np
# import matplotlib.pyplot as plt

# #############################################
# # 1. Load pre-processed JSON data           #
# #############################################
# # Replace with the path to your JSON file.
# json_file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\all_5400m_contours_15day_raw_all.json"
# with open(json_file_path, "r") as f:
#     json_data = json.load(f)

# #############################################
# # 2. Define dates for analysis              #
# #############################################
# # SSW event dates and non-event dates (as np.datetime64)
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
# # 3. Set analysis parameters                #
# #############################################
# # Specify the target longitude (ensure it matches one of the keys in the JSON data).
# lon_of_interest = -105.0  
# total_days = 14           # Analysis period: 15 days (day 0 through day 14)

# #############################################
# # 4. Function to compute maximum crossing    #
# #    latitude from JSON data                #
# #############################################
# def compute_max_contour_from_json(start_date_str, json_data, lon_of_interest):
#     """
#     For a given start_date (string), extract the precomputed contour crossing
#     data from the JSON and compute the maximum daily average crossing latitude
#     at the specified longitude over the period.
#     """
#     # Locate the entry matching the start date.
#     entry = next((item for item in json_data if item["start_date"] == start_date_str), None)
#     if entry is None:
#         print(f"No data found for start date {start_date_str}")
#         return None
    
#     period = entry["period"]
#     daily_avg = []
#     key = str(lon_of_interest)
#     for day in period:
#         crossings_dict = day.get("crossings", {})
#         crossings = crossings_dict.get(key, [])
#         if crossings:
#             daily_avg.append(np.mean(crossings))
#     return np.max(daily_avg) if daily_avg else None

# #############################################
# # 5. Compute maximum crossing latitudes      #
# #############################################
# results = []  # To hold dictionaries with date, max_lat, and type

# # Process SSW event dates.
# for event_date in ssw_event_dates:
#     date_str = str(event_date)
#     max_lat = compute_max_contour_from_json(date_str, json_data, lon_of_interest)
#     if max_lat is not None:
#         results.append({'date': event_date, 'max_lat': max_lat, 'type': 'event'})

# # Process non-event dates.
# for non_event_date in non_event_dates:
#     date_str = str(non_event_date)
#     max_lat = compute_max_contour_from_json(date_str, json_data, lon_of_interest)
#     if max_lat is not None:
#         results.append({'date': non_event_date, 'max_lat': max_lat, 'type': 'non-event'})

# #############################################
# # 6. Plot the maximum contour latitudes      #
# #    as a strip plot                        #
# #############################################
# fig, ax = plt.subplots(figsize=(6, 8))
# jitter_strength = 0.1  # Controls horizontal jitter for visualization

# for res in results:
#     if res['type'] == 'event':
#         x = 1 + np.random.uniform(-jitter_strength, jitter_strength)
#         color = 'blue'
#     else:
#         x = 2 + np.random.uniform(-jitter_strength, jitter_strength)
#         color = 'red'
#     ax.plot(x, res['max_lat'], 'o', markersize=8, color=color)
#     ax.text(x + 0.05, res['max_lat'], str(res['date'].astype('M8[D]')),
#             fontsize=9, verticalalignment='center', color=color)

# ax.set_xlim(0.5, 2.5)
# ax.set_xticks([1, 2])
# ax.set_xticklabels(['Event Years', 'Non-Event Years'])
# if results:
#     all_lats = [res['max_lat'] for res in results]
#     margin = 2
#     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
# ax.set_ylabel("Latitude (°)")
# ax.set_title(f"Density of Maximum 5400 m Contour Latitude\nat lon {lon_of_interest} (over {total_days} days)")
# plt.tight_layout()
# plt.show()

# import json
# import numpy as np
# import matplotlib.pyplot as plt

# #############################################
# # 1. Load pre-processed JSON data           #
# #############################################
# json_file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\all_5400m_contours_15day_raw_all.json"
# with open(json_file_path, "r") as f:
#     json_data = json.load(f)

# #############################################
# # 2. Define known event dates for classification
# #############################################
# # (These lists are used only for labeling the points.
# #  If your JSON file was generated from a different set of dates,
# #  then only those entries will appear.)
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
# # Create a set of strings for faster lookup.
# event_dates_set = {str(date) for date in ssw_event_dates}

# #############################################
# # 3. Analysis parameters                    #
# #############################################
# lon_of_interest = -105.0  # the x-coordinate where points will be plotted
# # (Your JSON already contains a 15-day period per start_date.)

# #############################################
# # 4. Process JSON data to extract maximum crossing latitudes
# #############################################
# results = []
# for entry in json_data:
#     start_date = entry.get("start_date")
#     period = entry.get("period", [])
#     daily_avgs = []
#     key = str(lon_of_interest)
#     for day in period:
#         # Each day should have a "crossings" dict.
#         crossings = day.get("crossings", {}).get(key, [])
#         if crossings:
#             daily_avgs.append(np.mean(crossings))
#     if daily_avgs:
#         max_lat = np.max(daily_avgs)
#         # Label as "event" if the start_date is in our known event dates, else "non-event".
#         event_type = "event" if start_date in event_dates_set else "non-event"
#         results.append({'date': start_date, 'max_lat': max_lat, 'type': event_type})

# #############################################
# # 5. Plot on a single vertical line         #
# #############################################
# fig, ax = plt.subplots(figsize=(6, 8))
# jitter_strength = 0.1  # horizontal jitter to prevent complete overlap

# # Draw a vertical reference line at lon_of_interest.
# ax.axvline(x=lon_of_interest, color='gray', linestyle='--')

# for res in results:
#     # Use slight random jitter so points don't plot exactly on top of each other.
#     x_jittered = lon_of_interest + np.random.uniform(-jitter_strength, jitter_strength)
#     color = 'blue' if res['type'] == 'event' else 'red'
#     ax.plot(x_jittered, res['max_lat'], 'o', markersize=8, color=color)
#     ax.text(x_jittered + 0.1, res['max_lat'], res['date'],
#             fontsize=9, verticalalignment='center', color=color)

# ax.set_xlim(lon_of_interest - 1, lon_of_interest + 1)
# if results:
#     all_lats = [res['max_lat'] for res in results]
#     margin = 2
#     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
# ax.set_xlabel("Longitude (°)")
# ax.set_ylabel("Latitude (°)")
# ax.set_title(f"Maximum 5400 m Contour Latitude\nat lon {lon_of_interest} (over 15 days)")
# plt.tight_layout()
# plt.show()

#!/usr/bin/env python3
import json
import numpy as np
import matplotlib.pyplot as plt

#############################################
# 1. Load pre-processed JSON data           #
#############################################
json_file_path = r"all_5400m_contours_15day_raw_all.json"
with open(json_file_path, "r") as f:
    json_data = json.load(f)

#############################################
# 2. Define known event dates for classification
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
# Create a set of strings for faster lookup.
event_dates_set = {str(date) for date in ssw_event_dates}

#############################################
# 3. Analysis parameters                    #
#############################################
lon_of_interest = -105.0  # the x-coordinate where points will be plotted

#############################################
# 4. Process JSON data to extract maximum crossing latitudes
#############################################
results = []
for entry in json_data:
    start_date = entry.get("start_date")
    period = entry.get("period", [])
    daily_avgs = []
    key = str(lon_of_interest)
    for day in period:
        crossings = day.get("crossings", {}).get(key, [])
        if crossings:
            daily_avgs.append(np.mean(crossings))
    if daily_avgs:
        max_lat = np.max(daily_avgs)
        event_type = "event" if start_date in event_dates_set else "non-event"
        results.append({'date': start_date, 'max_lat': max_lat, 'type': event_type})

#############################################
# 5. Plot on a single vertical line         #
#############################################
fig, ax = plt.subplots(figsize=(6, 8))
jitter_strength = 0.1  # horizontal jitter

# Draw a vertical reference line at lon_of_interest.
ax.axvline(x=lon_of_interest, color='gray', linestyle='--')

for res in results:
    x_jittered = lon_of_interest + np.random.uniform(-jitter_strength, jitter_strength)
    color = 'blue' if res['type'] == 'event' else 'red'
    ax.plot(x_jittered, res['max_lat'], 'o', markersize=8, color=color)
    ax.text(x_jittered + 0.1, res['max_lat'], res['date'],
            fontsize=9, verticalalignment='center', color=color)

ax.set_xlim(lon_of_interest - 1, lon_of_interest + 1)
if results:
    all_lats = [res['max_lat'] for res in results]
    margin = 10  # increased margin for a larger latitude range
    ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
ax.set_xlabel("Longitude (°)")
ax.set_ylabel("Latitude (°)")
ax.set_title(f"Maximum 5400 m Contour Latitude\nat lon {lon_of_interest} (over 15 days)")
plt.tight_layout()
plt.show()
