#!/usr/bin/env python3
# import json
# import numpy as np
# import matplotlib.pyplot as plt

# # Load the JSON file produced by your script.
# # (Make sure the path below points to "all_5400m_contours_15day_raw_all.json" or the file you intend to use.)
# with open(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\all_5400m_contours_15day_raw_all.json", "r") as f:
#     raw_data = json.load(f)

# # Example event and non-event date lists
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

# # Set the target longitude at which to compute the average crossing latitude.
# # (Ensure that this longitude matches one of the target values in your JSON file.)
# lon_of_interest = 35.0

# def compute_avg_crossing_from_raw(start_date_str, raw_data, lon_of_interest):
#     """
#     For a given start date (as a string), compute the average latitude at which
#     the 5400 m contour crosses the target longitude over the 15-day period,
#     using the 'crossings' data from the JSON file.
#     """
#     matching = [entry for entry in raw_data if entry["start_date"] == start_date_str]
#     if not matching:
#         print(f"No raw data found for start date {start_date_str}")
#         return None

#     period = matching[0]["period"]
#     daily_averages = []
#     for day in period:
#         # In the new JSON structure, each day has a 'crossings' dict
#         # with keys corresponding to the target longitudes (as strings).
#         crossings = day["crossings"].get(str(lon_of_interest), [])
#         if crossings:
#             daily_averages.append(np.mean(crossings))
#     return np.mean(daily_averages) if daily_averages else None

# # Process event dates.
# results = []
# for event_date in ssw_event_dates:
#     event_date_str = str(event_date)
#     avg_lat = compute_avg_crossing_from_raw(event_date_str, raw_data, lon_of_interest)
#     if avg_lat is not None:
#         results.append({'date': event_date, 'avg_lat': avg_lat, 'type': 'event'})

# # Process non-event dates.
# for non_event_date in non_event_dates:
#     non_event_date_str = str(non_event_date)
#     avg_lat = compute_avg_crossing_from_raw(non_event_date_str, raw_data, lon_of_interest)
#     if avg_lat is not None:
#         results.append({'date': non_event_date, 'avg_lat': avg_lat, 'type': 'non-event'})

# # Plot the average crossing latitudes as a strip plot.
# fig, ax = plt.subplots(figsize=(6, 8))
# jitter_strength = 0.1

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

# ax.set_xlim(0.5, 2.5)
# ax.set_xticks([1, 2])
# ax.set_xticklabels(['Event Years', 'Non-Event Years'])
# if results:
#     all_lats = [res['avg_lat'] for res in results]
#     margin = 2
#     ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
# ax.set_ylabel("Latitude (°)")
# ax.set_title(f"Density of Average 5400 m Crossing Latitude\n(over 15 days)")
# plt.tight_layout()
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

# Set the target longitude at which to compute the average crossing latitude.
# (Make sure that this longitude is among those used when generating the JSON file.)
lon_of_interest = 35.0

def compute_avg_crossing_from_raw(start_date_str, raw_data, lon_of_interest):
    """
    For a given start date (as a string), compute the average latitude at which
    the 5400 m contour crosses the target longitude over the 15-day period,
    using the 'crossings' data from the JSON file.
    """
    matching = [entry for entry in raw_data if entry["start_date"] == start_date_str]
    if not matching:
        print(f"No raw data found for start date {start_date_str}")
        return None

    period = matching[0]["period"]
    daily_averages = []
    for day in period:
        # Retrieve the crossing latitudes for the specified longitude.
        crossings = day["crossings"].get(str(lon_of_interest), [])
        if crossings:
            daily_averages.append(np.mean(crossings))
    return np.mean(daily_averages) if daily_averages else None

# Compute average crossing latitudes for event and non-event dates.
results = []
for event_date in ssw_event_dates:
    event_date_str = str(event_date)
    avg_lat = compute_avg_crossing_from_raw(event_date_str, raw_data, lon_of_interest)
    if avg_lat is not None:
        results.append({'date': event_date, 'avg_lat': avg_lat, 'type': 'event'})

for non_event_date in non_event_dates:
    non_event_date_str = str(non_event_date)
    avg_lat = compute_avg_crossing_from_raw(non_event_date_str, raw_data, lon_of_interest)
    if avg_lat is not None:
        results.append({'date': non_event_date, 'avg_lat': avg_lat, 'type': 'non-event'})

# Plot the average crossing latitudes as a strip plot.
fig, ax = plt.subplots(figsize=(6, 8))
jitter_strength = 0.1  # Controls horizontal spread for visualization.

for res in results:
    if res['type'] == 'event':
        x = 1 + np.random.uniform(-jitter_strength, jitter_strength)
        color = 'blue'
    else:
        x = 2 + np.random.uniform(-jitter_strength, jitter_strength)
        color = 'red'
    ax.plot(x, res['avg_lat'], 'o', markersize=8, color=color)
    ax.text(x + 0.05, res['avg_lat'], str(res['date'].astype('M8[D]')),
            fontsize=9, verticalalignment='center', color=color)

ax.set_xlim(0.5, 2.5)
ax.set_xticks([1, 2])
ax.set_xticklabels(['Event Years', 'Non-Event Years'])
if results:
    all_lats = [res['avg_lat'] for res in results]
    margin = 2
    ax.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
ax.set_ylabel("Latitude (°)")
ax.set_title("Density of Average 5400 m Crossing Latitude\n(over 15 days)")
plt.tight_layout()
plt.show()
