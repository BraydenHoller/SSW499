# #!/usr/bin/env python3
# """
# Extract 5400 m thickness contour crossings across the Northern Hemisphere
# for SSW and non-SSW dates, then export the results to NetCDF4 for efficient plotting.
# """
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt  # used for contour extraction

# # 1. Load and prepare the dataset
# DS_PATH = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All"  # adjust to your dataset path

# ds = xr.open_dataset(DS_PATH)

# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(
#             longitude=((ds.longitude + 180) % 360) - 180
#         ).sortby('longitude')
#     return ds

# ds = adjust_longitude(ds)

# # Subset for the two pressure levels; keep full spatial domain
# dz = ds['z'].sel(pressure_level=[500.0, 1000.0])

# # 2. Define SSW event & non-event dates
# ssw_event_dates = [
#     np.datetime64(d, 'D') for d in [
#         "2010-02-10","2010-03-24","2013-01-06","2018-02-12",
#         "2019-01-01","2021-01-05","2023-02-16","2000-03-20",
#         "2001-02-11","2001-12-30","2002-02-17","2003-01-18",
#         "2004-01-05","2006-01-21","2007-02-24","2008-02-22",
#         "2009-01-24","1998-12-15","1999-02-26","1980-02-29",
#         "1981-03-04","1981-12-04","1984-02-24","1985-01-01",
#         "1987-01-23","1987-12-08","1988-03-14","1989-02-21",
#         "1970-01-02","1971-01-18","1971-03-20","1973-01-31",
#         "1977-01-09","1979-02-22","1958-02-08","1960-01-17",
#         "1963-01-27","1965-12-16","1966-02-22","1968-01-07",
#         "1968-11-28","1969-03-13"
#     ]
# ]
# non_event_dates = [
#     np.datetime64(d, 'D') for d in [
#         "1959-01-30","1961-01-30","1962-01-30","1964-01-30",
#         "1967-01-30","1972-01-30","1974-01-30","1975-01-30",
#         "1976-01-30","1978-01-30","1982-01-30","1983-01-30",
#         "1986-01-30","1990-01-30","1991-01-30","1992-01-30",
#         "1993-01-30","1994-01-30","1995-01-30","1996-01-30",
#         "1997-01-30","2005-01-30","2011-01-30","2012-01-30",
#         "2014-01-30","2015-01-30","2016-01-30","2017-01-30",
#         "2020-01-30","2022-01-30"
#     ]
# ]
# all_dates = [(d, 'event') for d in ssw_event_dates] + [(d, 'non-event') for d in non_event_dates]

# # 3. Contour extraction helper
# def extract_contour_points(valid_date):
#     """
#     Return an (N,2) array of [lon, lat] for the 5400 m contour
#     across the full grid, then filter to lat ≥ 0 (NH).
#     """
#     data = dz.sel(valid_time=valid_date, method='nearest')
#     geo = data / 9.80665
#     thickness = geo.sel(pressure_level=500.0) - geo.sel(pressure_level=1000.0)

#     lon_arr = thickness.longitude.values
#     lat_arr = thickness.latitude.values
#     Z = thickness.values

#     fig, ax = plt.subplots()
#     cs = ax.contour(lon_arr, lat_arr, Z, levels=[5400])
#     plt.close(fig)

#     segs = cs.allsegs[0]
#     if not segs:
#         return np.empty((0,2))
#     pts = np.vstack(segs)
#     return pts[pts[:,1] >= 0]

# # 4. Flatten all crossings into a 1D table
# obs_dates = []
# obs_types = []
# obs_lons = []
# obs_lats = []
# for date_obj, dtype in all_dates:
#     pts = extract_contour_points(date_obj)
#     for lon, lat in pts:
#         obs_dates.append(date_obj)
#         obs_types.append(dtype)
#         obs_lons.append(lon)
#         obs_lats.append(lat)

# # 5. Build an xarray Dataset for NetCDF export
# dates_arr = np.array(obs_dates, dtype='datetime64[D]')
# types_arr = np.array(obs_types, dtype='U')

# ds_out = xr.Dataset({
#     'longitude': ('obs', np.array(obs_lons)),
#     'latitude': ('obs', np.array(obs_lats)),
#     'type': ('obs', types_arr)
# }, coords={
#     'date': ('obs', dates_arr)
# })

# # 6. Export to NetCDF4
# OUT_PATH = 'nh_5400m_contours_crossings.nc'
# ds_out.to_netcdf(OUT_PATH, format='NETCDF4')
# print(f"Exported {len(dates_arr)} contour points to {OUT_PATH}")

# #!/usr/bin/env python3
# """
# Extract the maximum latitude of the 5400 m thickness contour at each 10° longitude bin
# for each SSW and non-SSW date, then export the results to a NetCDF4 file for plotting.
# """
# import xarray as xr
# import numpy as np
# import matplotlib.pyplot as plt

# # 1. Load and prepare the dataset
# DS_PATH = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All"
# ds = xr.open_dataset(DS_PATH)

# def adjust_longitude(ds):
#     if ds.longitude.max() > 180:
#         ds = ds.assign_coords(
#             longitude=((ds.longitude + 180) % 360) - 180
#         ).sortby('longitude')
#     return ds

# ds = adjust_longitude(ds)

# # Subset geopotential height for two pressure levels
# dz = ds['z'].sel(pressure_level=[500.0, 1000.0])

# # 2. Define event & non-event dates
# ssw_dates = [
#     "2010-02-10","2010-03-24","2013-01-06","2018-02-12",
#     "2019-01-01","2021-01-05","2023-02-16","2000-03-20",
#     "2001-02-11","2001-12-30","2002-02-17","2003-01-18",
#     "2004-01-05","2006-01-21","2007-02-24","2008-02-22",
#     "2009-01-24","1998-12-15","1999-02-26","1980-02-29",
#     "1981-03-04","1981-12-04","1984-02-24","1985-01-01",
#     "1987-01-23","1987-12-08","1988-03-14","1989-02-21",
#     "1970-01-02","1971-01-18","1971-03-20","1973-01-31",
#     "1977-01-09","1979-02-22","1958-02-08","1960-01-17",
#     "1963-01-27","1965-12-16","1966-02-22","1968-01-07",
#     "1968-11-28","1969-03-13"
# ]
# non_dates = [
#     "1959-01-30","1961-01-30","1962-01-30","1964-01-30",
#     "1967-01-30","1972-01-30","1974-01-30","1975-01-30",
#     "1976-01-30","1978-01-30","1982-01-30","1983-01-30",
#     "1986-01-30","1990-01-30","1991-01-30","1992-01-30",
#     "1993-01-30","1994-01-30","1995-01-30","1996-01-30",
#     "1997-01-30","2005-01-30","2011-01-30","2012-01-30",
#     "2014-01-30","2015-01-30","2016-01-30","2017-01-30",
#     "2020-01-30","2022-01-30"
# ]
# all_dates = [(np.datetime64(d, 'D'), 'event') for d in ssw_dates] + [(np.datetime64(d, 'D'), 'non-event') for d in non_dates]

# # 3. Define longitude bins (every 10°)
# lon_bins = np.arange(-180, 181, 10)

# # 4. Helper to extract raw crossings
# def extract_contour_points(date):
#     # Select nearest time slice
#     data = dz.sel(valid_time=date, method='nearest')
#     # Convert geopotential to geometric height (m)
#     geo = data / 9.80665
#     # Compute thickness (500 hPa minus 1000 hPa)
#     thickness = geo.sel(pressure_level=500.0) - geo.sel(pressure_level=1000.0)

#     # Arrays for contour
#     lon_arr = thickness.longitude.values
#     lat_arr = thickness.latitude.values
#     Z = thickness.values

#     # Extract 5400 m contour segments
#     fig, ax = plt.subplots()
#     cs = ax.contour(lon_arr, lat_arr, Z, levels=[5400])
#     plt.close(fig)

#     segs = cs.allsegs[0]
#     if not segs:
#         return np.empty((0, 2))
#     pts = np.vstack(segs)
#     # Keep only Northern Hemisphere
#     return pts[pts[:,1] >= 0]

# # 5. Compute max latitude per date and lon bin
# dates_list = []
# types_list = []
# max_lat_matrix = np.full((len(all_dates), len(lon_bins)), np.nan)

# for i, (date_obj, dtype) in enumerate(all_dates):
#     pts = extract_contour_points(date_obj)
#     # Bin longitudes
#     if pts.size == 0:
#         continue
#     binned = (np.round(pts[:,0] / 10) * 10).astype(int)
#     for j, lb in enumerate(lon_bins):
#         mask = binned == lb
#         if mask.any():
#             max_lat_matrix[i, j] = np.max(pts[mask, 1])
#     dates_list.append(date_obj)
#     types_list.append(dtype)

# # Convert to arrays
# date_arr = np.array(dates_list, dtype='datetime64[D]')
# type_arr = np.array(types_list, dtype='U')

# # 6. Build xarray Dataset
# ds_out = xr.Dataset(
#     {
#         'max_latitude': (('date', 'longitude'), max_lat_matrix)
#     },
#     coords={
#         'date': date_arr,
#         'longitude': lon_bins,
#         'type': ('date', type_arr)
#     }
# )

# # 7. Export to NetCDF
# OUT_PATH = 'nh_5400m_maxlat_10deg_per_event.nc'
# ds_out.to_netcdf(OUT_PATH, format='NETCDF4')
# print(f"Exported max latitude per event to {OUT_PATH}")

#!/usr/bin/env python3
"""
Extract the average latitude of the 5400 m thickness contour at each 10° longitude bin
for each SSW and non-SSW date, then export to a NetCDF4 file for plotting.
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# 1. Load and prepare the dataset
DS_PATH = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All"
ds = xr.open_dataset(DS_PATH)

def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(
            longitude=((ds.longitude + 180) % 360) - 180
        ).sortby('longitude')
    return ds

ds = adjust_longitude(ds)

# Subset geopotential height for two pressure levels
dz = ds['z'].sel(pressure_level=[500.0, 1000.0])

# 2. Define event & non-event dates
ssw_dates = [
    "2010-02-10","2010-03-24","2013-01-06","2018-02-12",
    "2019-01-01","2021-01-05","2023-02-16","2000-03-20",
    "2001-02-11","2001-12-30","2002-02-17","2003-01-18",
    "2004-01-05","2006-01-21","2007-02-24","2008-02-22",
    "2009-01-24","1998-12-15","1999-02-26","1980-02-29",
    "1981-03-04","1981-12-04","1984-02-24","1985-01-01",
    "1987-01-23","1987-12-08","1988-03-14","1989-02-21",
    "1970-01-02","1971-01-18","1971-03-20","1973-01-31",
    "1977-01-09","1979-02-22","1958-02-08","1960-01-17",
    "1963-01-27","1965-12-16","1966-02-22","1968-01-07",
    "1968-11-28","1969-03-13"
]
non_dates = [
    "1959-01-30","1961-01-30","1962-01-30","1964-01-30",
    "1967-01-30","1972-01-30","1974-01-30","1975-01-30",
    "1976-01-30","1978-01-30","1982-01-30","1983-01-30",
    "1986-01-30","1990-01-30","1991-01-30","1992-01-30",
    "1993-01-30","1994-01-30","1995-01-30","1996-01-30",
    "1997-01-30","2005-01-30","2011-01-30","2012-01-30",
    "2014-01-30","2015-01-30","2016-01-30","2017-01-30",
    "2020-01-30","2022-01-30"
]
all_dates = [(np.datetime64(d, 'D'), 'event') for d in ssw_dates] + [(np.datetime64(d, 'D'), 'non-event') for d in non_dates]

# 3. Define longitude bins (every 10°)
lon_bins = np.arange(-180, 181, 10)

# 4. Helper to extract raw contour crossings
def extract_contour_points(date):
    data = dz.sel(valid_time=date, method='nearest')
    geo = data / 9.80665
    thickness = geo.sel(pressure_level=500.0) - geo.sel(pressure_level=1000.0)

    lon_arr = thickness.longitude.values
    lat_arr = thickness.latitude.values
    Z = thickness.values

    fig, ax = plt.subplots()
    cs = ax.contour(lon_arr, lat_arr, Z, levels=[5400])
    plt.close(fig)

    segs = cs.allsegs[0]
    if not segs:
        return np.empty((0, 2))
    pts = np.vstack(segs)
    return pts[pts[:,1] >= 0]

# 5. Compute average latitude per date and lon bin

dates_list = []
types_list = []
avg_lat_matrix = np.full((len(all_dates), len(lon_bins)), np.nan)

for i, (date_obj, dtype) in enumerate(all_dates):
    pts = extract_contour_points(date_obj)
    if pts.size == 0:
        continue
    binned = (np.round(pts[:,0] / 10) * 10).astype(int)
    for j, lb in enumerate(lon_bins):
        mask = binned == lb
        if mask.any():
            avg_lat_matrix[i, j] = np.mean(pts[mask, 1])
    dates_list.append(date_obj)
    types_list.append(dtype)

# Convert to arrays
date_arr = np.array(dates_list, dtype='datetime64[D]')
type_arr = np.array(types_list, dtype='U')

# 6. Build xarray Dataset
ds_out = xr.Dataset(
    {
        'average_latitude': (('date', 'longitude'), avg_lat_matrix)
    },
    coords={
        'date': date_arr,
        'longitude': lon_bins,
        'type': ('date', type_arr)
    }
)

# 7. Export to NetCDF
OUT_PATH = 'nh_5400m_avglat_10deg_per_event.nc'
ds_out.to_netcdf(OUT_PATH, format='NETCDF4')
print(f"Exported average latitude per event to {OUT_PATH}")


# -----------------------------------------------------------------------------

#!/usr/bin/env python3
"""
Extract the minimum latitude of the 5400 m thickness contour at each 10° longitude bin
for each SSW and non-SSW date, then export to a NetCDF4 file for plotting.
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# 1. Load and prepare the dataset
DS_PATH = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\SSW_5400Z_All"
ds = xr.open_dataset(DS_PATH)

def adjust_longitude(ds):
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(
            longitude=((ds.longitude + 180) % 360) - 180
        ).sortby('longitude')
    return ds

ds = adjust_longitude(ds)

# Subset geopotential height for two pressure levels
dz = ds['z'].sel(pressure_level=[500.0, 1000.0])

# 2. Define event & non-event dates
# (same lists as above)
# ... (reuse ssw_dates, non_dates, all_dates definitions)

# 3. Define longitude bins (every 10°)
lon_bins = np.arange(-180, 181, 10)

# 4. Helper to extract raw contour crossings (reuse extract_contour_points)
# ... (reuse function definition)

# 5. Compute minimum latitude per date and lon bin

dates_list = []
types_list = []
min_lat_matrix = np.full((len(all_dates), len(lon_bins)), np.nan)

for i, (date_obj, dtype) in enumerate(all_dates):
    pts = extract_contour_points(date_obj)
    if pts.size == 0:
        continue
    binned = (np.round(pts[:,0] / 10) * 10).astype(int)
    for j, lb in enumerate(lon_bins):
        mask = binned == lb
        if mask.any():
            min_lat_matrix[i, j] = np.min(pts[mask, 1])
    dates_list.append(date_obj)
    types_list.append(dtype)

# Convert to arrays
date_arr = np.array(dates_list, dtype='datetime64[D]')
type_arr = np.array(types_list, dtype='U')

# 6. Build xarray Dataset
ds_out = xr.Dataset(
    {
        'minimum_latitude': (('date', 'longitude'), min_lat_matrix)
    },
    coords={
        'date': date_arr,
        'longitude': lon_bins,
        'type': ('date', type_arr)
    }
)

# 7. Export to NetCDF
OUT_PATH = 'nh_5400m_minlat_10deg_per_event.nc'
ds_out.to_netcdf(OUT_PATH, format='NETCDF4')
print(f"Exported minimum latitude per event to {OUT_PATH}")
