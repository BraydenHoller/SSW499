# #!/usr/bin/env python3
# import json
# import numpy as np
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs

# #############################################
# # 1. Load precomputed raw contour JSON data  #
# #############################################
# json_filename = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\all_5400m_contours_15day_raw_all.json"
# with open(json_filename, "r") as f:
#     raw_data = json.load(f)

# #############################################
# # 2. Define the combined date list           #
# #############################################
# # (Event + Non-event), each as np.datetime64
# all_dates = [
#     # ... your full combined date list ...
# ]

# #############################################
# # 3. Set the 15-day period parameters        #
# #############################################
# total_days = 14  # day 0 through day 14
# num_days = total_days + 1

# #############################################
# # 4. Target longitudes (example: every 30°)  #
# #############################################
# target_lons = [str(lon) for lon in np.arange(-180, 181, 30)]

# #############################################
# # 5. Function to compute average crossing latitudes
# #############################################
# def compute_avg_crossings_for_date(start_date_str, raw_data):
#     matching = [entry for entry in raw_data if entry["start_date"].startswith(start_date_str)]
#     if not matching:
#         print(f"No raw data found for start date {start_date_str}")
#         return None
#     period = matching[0]["period"]
    
#     # For each target longitude, collect daily averages
#     avg_crossings = {lon: [] for lon in target_lons}
#     for day in period:
#         crossings = day.get("crossings", {})
#         for lon in target_lons:
#             # If the key exists, compute the day's average for that lon
#             if lon in crossings and crossings[lon]:
#                 avg_crossings[lon].append(np.mean(crossings[lon]))
#             else:
#                 avg_crossings[lon].append(np.nan)
#     # Now average over days (ignoring NaNs)
#     for lon in avg_crossings:
#         vals = np.array(avg_crossings[lon])
#         valid = vals[~np.isnan(vals)]
#         if valid.size > 0:
#             avg_crossings[lon] = np.mean(valid)
#         else:
#             avg_crossings[lon] = np.nan
#     return avg_crossings

# #############################################
# # 6. (Optional) Distinguish event vs. non-event
# #############################################
# event_dates_set = {
#     # Convert your event dates to string
#     "2010-02-10", "2010-03-24", "2013-01-06", "2018-02-12",
#     # ... etc. ...
# }
# # If a date is in event_dates_set => "event", else => "non-event"

# #############################################
# # 7. Collect data for plotting
# #############################################
# plot_data = []  # (lons, lats, date_str, dtype)

# for date_obj in all_dates:
#     date_str = str(date_obj).split("T")[0]
#     avg_cross = compute_avg_crossings_for_date(date_str, raw_data)
#     if avg_cross is None:
#         continue
    
#     # Convert dictionary to sorted arrays
#     # Make sure the format for keys matches how they were stored in the JSON
#     sorted_keys = sorted(avg_cross.keys(), key=lambda x: float(x))
#     lons = np.array([float(k) for k in sorted_keys])
#     lats = np.array([avg_cross[k] for k in sorted_keys])
    
#     # Determine event or non-event
#     dtype = "event" if date_str in event_dates_set else "non-event"
    
#     # Debug: print a small sample
#     print(f"Date: {date_str}, type: {dtype}, #points: {len(lons)}")
#     print("Sample lons:", lons[:5], "Sample lats:", lats[:5])
    
#     plot_data.append((lons, lats, date_str, dtype))

# #############################################
# # 8. Plot on a polar stereographic projection
# #############################################
# proj = ccrs.NorthPolarStereo()
# fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': proj})

# # Expand latitudinal extent down to 20°N
# ax.set_extent([-180, 180, 20, 90], crs=ccrs.PlateCarree())

# ax.coastlines()
# ax.gridlines()

# # Plot each date's average contour crossing
# for lons, lats, date_str, dtype in plot_data:
#     color = 'blue' if dtype == "event" else 'red'
    
#     # Filter out NaN or impossible lat values
#     mask = ~np.isnan(lats)
#     if not np.any(mask):
#         continue  # skip if all NaN
#     ax.plot(lons[mask], lats[mask], marker='o', linestyle='-',
#             color=color, transform=ccrs.PlateCarree(), label=dtype)

# # De-duplicate legend entries
# handles, labels = ax.get_legend_handles_labels()
# unique = dict(zip(labels, handles))
# ax.legend(unique.values(), unique.keys(), loc='lower left')

# ax.set_title("Average 5400 m Contour Latitude at Every 30°\n(15-day period starting on each date)")
# plt.show()

#!/usr/bin/env python3
import json
import numpy as np
import matplotlib.pyplot as plt

# Load the JSON file with raw contour crossing data.
with open(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\all_5400m_contours_15day_raw_all.json", "r") as f:
    raw_data = json.load(f)

# Define SSW event and non-event dates.
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

# Define target longitudes (every 30° across the globe).
target_lons = np.arange(-180, 181, 30)

def compute_avg_crossing_for_date_and_lon(start_date_str, raw_data, target_lon):
    """
    For a given start date (string) and target longitude,
    compute the average latitude at which the 5400 m contour crosses that longitude
    over the 15-day period. Returns np.nan if no valid data is found.
    """
    matching = [entry for entry in raw_data if entry["start_date"] == start_date_str]
    if not matching:
        print(f"No raw data found for start date {start_date_str}")
        return np.nan

    period = matching[0]["period"]
    daily_averages = []
    for day in period:
        # Look up the crossing latitudes for this target longitude.
        crossings = day["crossings"].get(str(target_lon), [])
        if crossings:
            daily_averages.append(np.mean(crossings))
    return np.mean(daily_averages) if daily_averages else np.nan

# Collect average crossing values for each target longitude for each group.
event_results = {lon: [] for lon in target_lons}
non_event_results = {lon: [] for lon in target_lons}

# Process event dates.
for date in ssw_event_dates:
    date_str = str(date)
    for lon in target_lons:
        avg_lat = compute_avg_crossing_for_date_and_lon(date_str, raw_data, lon)
        if not np.isnan(avg_lat):
            event_results[lon].append(avg_lat)

# Process non-event dates.
for date in non_event_dates:
    date_str = str(date)
    for lon in target_lons:
        avg_lat = compute_avg_crossing_for_date_and_lon(date_str, raw_data, lon)
        if not np.isnan(avg_lat):
            non_event_results[lon].append(avg_lat)

# Compute the overall average for each target longitude for each group.
event_profile = {lon: np.mean(event_results[lon]) if event_results[lon] else np.nan for lon in target_lons}
non_event_profile = {lon: np.mean(non_event_results[lon]) if non_event_results[lon] else np.nan for lon in target_lons}

# Plot the profiles.
plt.figure(figsize=(10, 6))
plt.plot(list(event_profile.keys()), list(event_profile.values()), marker='o', label='Event Dates')
plt.plot(list(non_event_profile.keys()), list(non_event_profile.values()), marker='o', label='Non-Event Dates')
plt.xlabel("Longitude (°)")
plt.ylabel("Average 5400 m Contour Crossing Latitude (°)")
plt.title("Global Profile of Average 5400 m Contour Crossing Latitude\n(15-day averages at every 30° longitude)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()