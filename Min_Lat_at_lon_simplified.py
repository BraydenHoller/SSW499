#!/usr/bin/env python3
"""
Load the precomputed 5400 m contour crossings from NetCDF and plot
minimum latitude at every 10° longitude for SSW and non‑SSW dates
on a Plate Carrée map.
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# 1. Load the NetCDF file (must contain a 'min_latitude' variable)
DS_CROSS = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\nh_5400m_minlat_10deg_per_event.nc"
ds = xr.open_dataset(DS_CROSS)

# 2. Extract arrays
lons      = ds['longitude'].values        # shape (num_lon_bins,)
min_lat   = ds['minimum_latitude'].values     # shape (num_dates, num_lon_bins)
types     = ds['type'].values             # shape (num_dates,), entries 'event' or 'non-event'

# 3. Unroll into per‑type point lists
event_lons, event_lats = [], []
non_event_lons, non_event_lats = [], []

# find indices of each type
ev_idx = np.where(types == 'event')[0]
ne_idx = np.where(types == 'non-event')[0]

# for SSW events
ev_matrix = min_lat[ev_idx, :]
ev_flat   = ev_matrix.flatten()
ev_lons_flat = np.tile(lons, ev_matrix.shape[0])
mask_ev   = ~np.isnan(ev_flat)
event_lons = ev_lons_flat[mask_ev].tolist()
event_lats = ev_flat[mask_ev].tolist()

# for non‑SSW dates
ne_matrix = min_lat[ne_idx, :]
ne_flat   = ne_matrix.flatten()
ne_lons_flat = np.tile(lons, ne_matrix.shape[0])
mask_ne   = ~np.isnan(ne_flat)
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


# Scatter the min‑lat points
ax.scatter(event_lons, event_lats,
           c='blue', s=30, marker='o',
           transform=ccrs.PlateCarree(),
           label='SSW Events')
ax.scatter(non_event_lons, non_event_lats,
           c='red', s=30, marker='o',
           transform=ccrs.PlateCarree(),
           label='Non‑SSW Dates')

ax.legend(loc='lower left')
ax.set_title('Minimum 5400 m Contour Latitude at 10° Lon Intervals\n'
             '(Blue = SSW Events, Red = Non‑SSW Dates)')
plt.tight_layout()
plt.show()
