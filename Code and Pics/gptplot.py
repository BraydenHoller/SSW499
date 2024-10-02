import xarray as xr
import pandas as pd
import plotly.express as px

# Load the NetCDF datafile without decoding times
data = xr.open_dataset("SSWC_v1.0_geopClimMean_ERAi_s19790101_e20140630_c20160701.nc", decode_times=False)

# Extract the pressure data over time and reduce the dataset size for visualization
pressure_data = data.geopClimMean.isel(lon=slice(None, None, 5), lat=slice(None, None, 5)).mean(dim='pres')
pressure_data_df = pressure_data.to_dataframe().reset_index()

# Create an interactive scatter geo plot
fig = px.scatter_geo(
    pressure_data_df, 
    lat='lat', 
    lon='lon', 
    color='geopClimMean',
    animation_frame='timeClim',
    projection="natural earth",
    title="Geopotential Height Changes Over Time"
)
fig.update_layout(
    geo=dict(
        showland=True,
        landcolor="rgb(217, 217, 217)",
    )
)
fig.show()
# # Save the figure to an HTML file
# html_file_path = "geopotential_height_changes_over_time.html"
# fig.write_html(html_file_path)

# # Display the saved file path (optional)
# print(f"Interactive plot saved as {html_file_path}")
