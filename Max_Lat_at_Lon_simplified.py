# #!/usr/bin/env python3
# import json
# import numpy as np
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature

# # ——— 1. Load the JSON data ——————————————————————
# JSON_PATH = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\all_5400m_contours_15day_raw_all.json"
# with open(JSON_PATH, "r") as f:
#     records = json.load(f)

# # ——— 2. Define SSW event & non-event dates ——————————
# EVENT_DATES = {
#     "2010-02-10","2010-03-24","2013-01-06","2018-02-12","2019-01-01",
#     "2021-01-05","2023-02-16","2000-03-20","2001-02-11","2001-12-30",
#     "2002-02-17","2003-01-18","2004-01-05","2006-01-21","2007-02-24",
#     "2008-02-22","2009-01-24","1998-12-15","1999-02-26","1980-02-29",
#     "1981-03-04","1981-12-04","1984-02-24","1985-01-01","1987-01-23",
#     "1987-12-08","1988-03-14","1989-02-21","1970-01-02","1971-01-18",
#     "1971-03-20","1973-01-31","1977-01-09","1979-02-22","1958-02-08",
#     "1960-01-17","1963-01-27","1965-12-16","1966-02-22","1968-01-07",
#     "1968-11-28","1969-03-13"
# }

# NON_EVENT_DATES = {
#     "1959-01-30","1961-01-30","1962-01-30","1964-01-30","1967-01-30",
#     "1972-01-30","1974-01-30","1975-01-30","1976-01-30","1978-01-30",
#     "1982-01-30","1983-01-30","1986-01-30","1990-01-30","1991-01-30",
#     "1992-01-30","1993-01-30","1994-01-30","1995-01-30","1996-01-30",
#     "1997-01-30","2005-01-30","2011-01-30","2012-01-30","2014-01-30",
#     "2015-01-30","2016-01-30","2017-01-30","2020-01-30","2022-01-30"
# }

# # ——— 3. Target longitudes every 10° (as floats!) —————————————
# target_lons = np.arange(-180.0, 181.0, 10.0)

# # ——— 4. Collect max‐lat per record & lon —————————————————
# event_lons, event_lats = [], []
# non_event_lons, non_event_lats = [], []

# for rec in records:
#     date_str = rec["start_date"]
#     if date_str in EVENT_DATES:
#         is_event = True
#     elif date_str in NON_EVENT_DATES:
#         is_event = False
#     else:
#         continue

#     for lon in target_lons:
#         key = f"{lon:.1f}"              # match JSON key format
#         all_cross = []
#         for day in rec.get("period", []):
#             all_cross.extend(day.get("crossings", {}).get(key, []))
#         if not all_cross:
#             continue
#         max_lat = max(all_cross)

#         if is_event:
#             event_lons.append(lon)
#             event_lats.append(max_lat)
#         else:
#             non_event_lons.append(lon)
#             non_event_lats.append(max_lat)

# # ——— 5. Plot on Plate Carrée (0–90° N) —————————————————
# fig = plt.figure(figsize=(10, 8))
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_extent([-180, 180, 0, 90], ccrs.PlateCarree())

# ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="lightgray")
# ax.add_feature(cfeature.COASTLINE.with_scale("110m"))
# ax.add_feature(cfeature.BORDERS.with_scale("110m"), linestyle=":")
# ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# ax.scatter(event_lons, event_lats,
#            c="blue", s=25, marker="o",
#            transform=ccrs.PlateCarree(),
#            label="SSW Events")
# ax.scatter(non_event_lons, non_event_lats,
#            c="red", s=25, marker="o",
#            transform=ccrs.PlateCarree(),
#            label="Non-SSW Dates")

# ax.legend(loc="lower left")
# ax.set_title(
#     "Max 5400 m Contour Latitude at 10° Lon Intervals\n"
#     "(Blue = SSW Events, Red = Non-SSW Dates)",
#     pad=12
# )
# plt.tight_layout()
# plt.show()

# #!/usr/bin/env python3
# """
# Load the precomputed 5400 m contour crossings from NetCDF and plot
# max latitude at every 10° longitude for SSW and non-SSW dates on a Plate Carrée map
# """
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature

# # 1. Load the NetCDF file
# DS_CROSS = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\nh_5400m_maxlat_10deg_per_event.nc"
# ds = xr.open_dataset(DS_CROSS)

# # Extract arrays
# lons = ds['longitude'].values
# lats = ds['latitude'].values
# types = ds['type'].values  # strings 'event' or 'non-event'

# # 2. Bin longitudes to nearest 10°
# lon_bins = np.arange(-180, 181, 10)
# binned = (np.round(lons / 10) * 10).astype(int)

# # 3. Compute max latitude per bin per type
# event_lons, event_lats = [], []
# non_event_lons, non_event_lats = [], []

# for lb in lon_bins:
#     mask_ev = (types == 'event') & (binned == lb)
#     if mask_ev.any():
#         event_lons.append(lb)
#         event_lats.append(np.max(lats[mask_ev]))

#     mask_non = (types == 'non-event') & (binned == lb)
#     if mask_non.any():
#         non_event_lons.append(lb)
#         non_event_lats.append(np.max(lats[mask_non]))

# # 4. Plot
# fig = plt.figure(figsize=(10, 8))
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_extent([-180, 180, 0, 90], ccrs.PlateCarree())

# # Basemap features
# ax.add_feature(cfeature.LAND.with_scale('110m'), facecolor='lightgray')
# ax.add_feature(cfeature.COASTLINE.with_scale('110m'))
# ax.add_feature(cfeature.BORDERS.with_scale('110m'), linestyle=':')
# ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# # Scatter the max lat points
# ax.scatter(event_lons, event_lats,
#            c='blue', s=50, marker='o',
#            transform=ccrs.PlateCarree(),
#            label='SSW Events')
# ax.scatter(non_event_lons, non_event_lats,
#            c='red', s=50, marker='o',
#            transform=ccrs.PlateCarree(),
#            label='Non-SSW Dates')

# ax.legend(loc='lower left')
# ax.set_title('Max 5400 m Contour Latitude at 10° Lon Intervals\n(Blue = SSW Events, Red = Non-SSW Dates)')
# plt.tight_layout()
# plt.show()

"""
Load the precomputed 5400 m contour crossings from NetCDF and plot
max latitude at every 10° longitude for SSW and non‑SSW dates on a Plate Carrée map
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# 1. Load the NetCDF file
DS_CROSS = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\nh_5400m_maxlat_10deg_per_event.nc"
ds = xr.open_dataset(DS_CROSS)

# 2. Extract arrays
lons = ds['longitude'].values                         # shape (num_lon_bins,)
max_lat = ds['max_latitude'].values                   # shape (num_dates, num_lon_bins)
types = ds['type'].values                             # shape (num_dates,), entries 'event' or 'non-event'

# 3. Unroll into per‑type point lists
event_lons = []
event_lats = []
non_event_lons = []
non_event_lats = []

# find the date‑indices for each type
ev_idx = np.where(types == 'event')[0]
ne_idx = np.where(types == 'non-event')[0]

# for events
ev_matrix = max_lat[ev_idx, :]                        # (num_events, num_lon_bins)
ev_flat = ev_matrix.flatten()
# repeat lons once per event
ev_lons_flat = np.tile(lons, ev_matrix.shape[0])
# mask out missing data
mask_ev = ~np.isnan(ev_flat)
event_lons = ev_lons_flat[mask_ev].tolist()
event_lats = ev_flat[mask_ev].tolist()

# for non-events
ne_matrix = max_lat[ne_idx, :]
ne_flat = ne_matrix.flatten()
ne_lons_flat = np.tile(lons, ne_matrix.shape[0])
mask_ne = ~np.isnan(ne_flat)
non_event_lons = ne_lons_flat[mask_ne].tolist()
non_event_lats = ne_flat[mask_ne].tolist()

# 4. Plot
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, 0, 90], ccrs.PlateCarree())

# Basemap features
ax.add_feature(cfeature.LAND.with_scale('110m'), facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax.add_feature(cfeature.BORDERS.with_scale('110m'), linestyle=':')
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

import rasterio

with rasterio.open(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\PermaFrost\llipa.byte") as perma_ds:
    perma_data    = perma_ds.read(1)
    perma_mask = np.where(perma_data > 0, perma_data, np.nan)
    lon_min_p, lon_max_p = perma_ds.bounds.left, perma_ds.bounds.right
    lat_min_p, lat_max_p = perma_ds.bounds.bottom, perma_ds.bounds.top

ax.imshow(
    perma_mask,
    origin='upper',
    extent=[lon_min_p, lon_max_p, lat_min_p, lat_max_p],
    transform=ccrs.PlateCarree(),
    alpha=0.75,
    interpolation='nearest',
    zorder=1
)

# Scatter the max‑lat points for each event/non‑event
ax.scatter(event_lons, event_lats,
           c='blue', s=30, marker='o',
           transform=ccrs.PlateCarree(),
           label='SSW Events')
ax.scatter(non_event_lons, non_event_lats,
           c='red', s=30, marker='o',
           transform=ccrs.PlateCarree(),
           label='Non-SSW Dates')

ax.legend(loc='lower left')
ax.set_title('Max 5400 m Contour Latitude at 10° Lon Intervals\n'
             '(Blue = SSW Events, Red = Non-SSW Dates)')
plt.tight_layout()
plt.show()
