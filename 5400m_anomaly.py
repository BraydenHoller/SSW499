# Makes a line-type anomaly plot of the 5400 m contour
# #!/usr/bin/env python3
# """
# Plot the average SSW-induced poleward shift of the 5400 m contour
# as a function of longitude, relative to the non-SSW baseline.
# """
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt

# # ─── USER: point this at your NetCDF of crossings ───────────────────────────
# DS_CROSS = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\nh_5400m_contours_crossings.nc"
# # ───────────────────────────────────────────────────────────────────────────────

# # 1) load dataset
# ds = xr.open_dataset(DS_CROSS)

# # 2) pull out flat arrays
# lons  = ds['longitude'].values         # e.g. [ -180, -180, -170, -170, ... ]
# lats  = ds['latitude'].values          # same shape as lons
# types = ds['type'].astype(str).values  # 'event' vs 'non-event'

# # 3) find the unique longitudes (should be your 10° bins)
# unique_lons = np.sort(np.unique(lons))

# # prepare arrays
# mean_non = np.full_like(unique_lons, np.nan, dtype=float)
# mean_ssw = np.full_like(unique_lons, np.nan, dtype=float)
# anom     = np.full_like(unique_lons, np.nan, dtype=float)

# # 4) for each lon, compute mean non-SSW and mean SSW lat, then anomaly
# for i, lonval in enumerate(unique_lons):
#     m_non = (lons == lonval) & (types == 'non-event')
#     m_ssw = (lons == lonval) & (types == 'event')
#     if np.any(m_non):
#         mean_non[i] = np.nanmean(lats[m_non])
#     if np.any(m_ssw):
#         mean_ssw[i] = np.nanmean(lats[m_ssw])
#     # only compute anomaly where both exist
#     if np.isfinite(mean_non[i]) and np.isfinite(mean_ssw[i]):
#         anom[i] = mean_ssw[i] - mean_non[i]

# # 5) plot
# fig, ax = plt.subplots(figsize=(14,5))

# # plot mean anomaly
# ax.plot(unique_lons, anom, '-o', color='C0', lw=2, label='Mean SSW anomaly')

# # zero line
# ax.axhline(0, color='k', lw=1, ls='--', label='No shift')

# # shade
# ax.fill_between(unique_lons, anom, 0,
#                 where=anom>0, interpolate=True,
#                 color='C0', alpha=0.3, label='Poleward shift')
# ax.fill_between(unique_lons, anom, 0,
#                 where=anom<0, interpolate=True,
#                 color='C3', alpha=0.3, label='Equatorward shift')

# # format
# ax.set_xlim(-180, 180)
# ax.set_xticks(np.arange(-180, 181, 30))
# ax.set_xlabel('Longitude (°)')
# ax.set_ylabel('Latitude anomaly (°) relative to non-SSW')
# ax.set_title('Average Poleward Shift of 5400 m Contour during SSWs')
# ax.grid(True, ls=':')
# ax.legend(loc='upper right')

# plt.tight_layout()
# plt.show()

# Poifect Plot
#!/usr/bin/env python3
"""
Plot the average poleward (or equatorward) shift
of the 5400 m contour during SSW events, relative
to the non-SSW baseline, on a Plate Carrée map—
covering the full Northern Hemisphere.
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import TwoSlopeNorm

# ─── USER: point this at your average-latitudes NetCDF ───────
DS_AVG = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\nh_5400m_avglat_10deg_per_event.nc"
# ───────────────────────────────────────────────────────────────

# 1) load dataset
ds = xr.open_dataset(DS_AVG)

# 2) pull out variables
avglat = ds['average_latitude'].values    # shape (n_events, n_lon)
lons    = ds['longitude'].values           # shape (n_lon,)
types   = ds['type'].values.astype(str)    # e.g. 'event' vs 'non-event'
types_l = np.char.lower(types)

# 3) build masks
mask_non = (types_l == 'non-event') | (types_l == 'non-ssw')
mask_ssw = (types_l == 'event')     | (types_l == 'ssw')

# 4) compute means
mean_non = np.nanmean(avglat[mask_non, :], axis=0)  # (n_lon,)
mean_ssw = np.nanmean(avglat[mask_ssw, :], axis=0)  # (n_lon,)

# 5) anomaly = SSW minus non-SSW
anom = mean_ssw - mean_non

# 6) set up map (full NH)
fig = plt.figure(figsize=(12,6))
ax  = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-180, 180,   0,  90], crs=ccrs.PlateCarree())

# add background
ax.add_feature(cfeature.LAND,      facecolor='lightgray', zorder=0)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

# optional: draw gridlines with labels
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels    = False
gl.right_labels  = False
gl.xlabel_style = {'size':10}
gl.ylabel_style = {'size':10}

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

# 7) scatter at each lon, colored by anomaly
norm = TwoSlopeNorm(vmin=anom.min(), vcenter=0, vmax=anom.max())
sc   = ax.scatter(lons, mean_non,
                  c=anom, cmap='RdBu_r', norm=norm,
                  s=80, edgecolor='k',
                  transform=ccrs.PlateCarree())

# 8) add a horizontal colorbar
cbar = plt.colorbar(sc, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Mean SSW – non-SSW latitude (°)')

# 9) finishing touches
ax.set_title('Average Poleward Shift of 5400 m Contour during SSWs\n(Full Northern Hemisphere)')
plt.tight_layout()
plt.show()
