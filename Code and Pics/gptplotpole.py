import xarray as xr
import pandas as pd
import plotly.express as px

# Load the NetCDF datafile without decoding times
file_path = 'SSWC_v1.0_geopClimMean_ERAi_s19790101_e20140630_c20160701.nc'
data = xr.open_dataset(file_path, decode_times=False)

# Extract the pressure data over time and reduce the dataset size for visualization
pressure_data = data.geopClimMean.isel(lon=slice(None, None, 5), lat=slice(None, None, 5)).mean(dim='pres')
pressure_data_df = pressure_data.to_dataframe().reset_index()

# Create an interactive scatter geo plot with a North Pole centered projection
fig = px.scatter_geo(
    pressure_data_df, 
    lat='lat', 
    lon='lon', 
    color='geopClimMean',
    animation_frame='timeClim',
    projection="orthographic",
    title="Geopotential Height Changes Over Time (North Pole Centered)"
)
fig.update_geos(
    center={"lat": 90, "lon": 0},
    projection_scale=1  # Adjust the scale for better visibility
)
fig.show()

# Save the figu
