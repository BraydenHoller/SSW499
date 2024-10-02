# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:44:40 2024

@author: Kevin.Burris
"""

# from netCDF4 import Dataset
import xarray as xr
import matplotlib.pyplot as plt

from pathlib import Path

file = Path('C:/Users/Kevin.Burris/OneDrive - afacademy.af.edu/Documents/github/class_examples/met499/data/SSWC_v1.0_geopClimMean_ERAi_s19790101_e20140630_c20160701.nc')
ds = xr.open_dataset(file)

# check out https://docs.xarray.dev/en/stable/user-guide/plotting.html#maps
import cartopy.crs as ccrs
# The cartopy coordinate reference system provides the tools to reproject data
# onto different maps

# Here I specify what projection I want to use for my final map
map_proj = ccrs.LambertConformal(central_longitude=-95, central_latitude=45)

print(ds['geopClimMean']) # What does my data look like
'''
<xarray.DataArray 'geopClimMean' (timeClim: 366, pres: 37, lat: 73, lon: 144)> Size: 1GB
[142353504 values with dtype=float64]
Coordinates:
  * timeClim  (timeClim) object 3kB 1979-07-01 00:00:00 ... 1980-06-30 00:00:00
  * lon       (lon) float32 576B -180.0 -177.5 -175.0 ... 172.5 175.0 177.5
  * lat       (lat) float32 292B -90.0 -87.5 -85.0 -82.5 ... 82.5 85.0 87.5 90.0
  * pres      (pres) float32 148B 1.0 2.0 3.0 5.0 ... 925.0 950.0 975.0 1e+03
Attributes:
    long_name:      Geopotential height climatological mean
    standard_name:  geopotential_height
    units:          m
    valid_range:    [-32768  32766]
    cell_methods:   area: point (interval: 0.70312 degrees_east 0.70167 degre...
    comment:        Climatology runs from July 01 - June 30
'''
# The first line tells me there is a 4 dimensional array, and the dimensions are
# timeClim, pressure, latitude, and longitude

print(ds['geopClimMean'][0,30,:,:])
# This will show me the Geopotential height climatological mean during the 
# first time on the 30th pressure level for all latitudes and longitudes
'''
<xarray.DataArray 'geopClimMean' (lat: 73, lon: 144)> Size: 84kB
[10512 values with dtype=float64]
Coordinates:
    timeClim  object 8B 1979-07-01 00:00:00
  * lon       (lon) float32 576B -180.0 -177.5 -175.0 ... 172.5 175.0 177.5
  * lat       (lat) float32 292B -90.0 -87.5 -85.0 -82.5 ... 82.5 85.0 87.5 90.0
    pres      float32 4B 850.0
Attributes:
    long_name:      Geopotential height climatological mean
    standard_name:  geopotential_height
    units:          m
    valid_range:    [-32768  32766]
    cell_methods:   area: point (interval: 0.70312 degrees_east 0.70167 degre...
    comment:        Climatology runs from July 01 - June 30
'''

# Assign some data to a placeholder. In this case, the first time and the last 
# pressure level for all lat/lon
clim_mean_test = ds['geopClimMean'][0,-1,:,:]

# Call the plot method for that dataset
p = clim_mean_test.plot(
    transform=ccrs.PlateCarree(),  # the data's projection
    # check out https://en.wikipedia.org/wiki/Equirectangular_projection
    subplot_kws={"projection": map_proj}, # The projection I want the map in
)
# I want to add features to the map, so get the current drawing axes, gca()
ax = plt.gca()
# Add coastlines
ax.coastlines('50m')