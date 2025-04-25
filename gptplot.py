# import xarray as xr
# import pandas as pd
# import plotly.express as px
# import numpy as np
# from datetime import datetime
# import plotly.graph_objects as go

# # Step 1: Load the dataset
# file_path = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Snow_2223.nc"
# dataset = xr.open_dataset(file_path)

# # Step 2: Verify the 'valid_time' coordinate and use it directly if it’s already in datetime format
# times_converted = dataset['valid_time'].values  # Get the time values directly

# # Step 3: Select the data between 2022-07-01 and 2023-07-01
# start_date = "2022-07-01"
# end_date = "2023-07-01"
# start_datetime = np.datetime64(start_date)
# end_datetime = np.datetime64(end_date)

# # Slice the dataset to the desired time range for snowfall data ('sf')
# swe_data = dataset['sf'].sel(valid_time=slice(start_datetime, end_datetime))

# # Step 4: Subset the data to the continental U.S. lat/lon ranges
# us_swe_data = swe_data.sel(latitude=slice(49.50, 24.25), longitude=slice(235, 293), method=None)

# # Step 5: Convert SWE to Snowfall (using a ratio of 1:10)
# us_snowfall_data = us_swe_data / 0.1  # Assuming 1mm of SWE equals 10mm of snowfall

# # Step 6: Calculate cumulative snowfall at each grid point
# cumulative_snowfall = us_snowfall_data.cumsum(dim='valid_time')

# # Step 7: Convert the cumulative snowfall data into a Pandas DataFrame for plotting
# # snowfall_df = cumulative_snowfall.to_dataframe().reset_index()

# # Debugging: Print out some info to check the data
# print("Data summary:")
# print(cumulative_snowfall.head())  # Print the first few rows of the DataFrame

# # Check for missing values in important columns
# print("\nChecking for missing data:")
# print(cumulative_snowfall.isnull().sum())  # Check for NaN values in the DataFrame

# # # Step 8: Create an interactive scatter geo plot with Plotly to visualize cumulative snowfall
# # fig = go.Figure(data=go.Scattergeo(
# #     cumulative_snowfall, 
# #     lat = cumulative_snowfall['latitude'], 
# #     lon = cumulative_snowfall['longitude'], 
# #     color = cumulative_snowfall,  # Color by cumulative snowfall
# #     animation_frame = cumulative_snowfall['valid_time'],  # Animate over time
# #     projection="natural earth",
# #     title=f"Cumulative Snowfall in the Continental U.S. from {start_date} to {end_date}",
# #     labels={'sf': 'Cumulative Snowfall (m)'}
# # ))

# # # Customize the layout
# # fig.update_layout(
# #     geo=dict(
# #         showland=True,
# #         landcolor="rgb(217, 217, 217)"
# #     ),
# #     coloraxis_colorbar=dict(
# #         title="Cumulative Snowfall (m)"
# #     )
# # )

# # # Show plot
# # fig.show()
# # Extract the latitude, longitude, and data values from the DataArray
# latitudes = cumulative_snowfall.coords['latitude'].values
# longitudes = cumulative_snowfall.coords['longitude'].values
# values = cumulative_snowfall.values

# lat_min, lat_max = 24, 50
# lon_min, lon_max = 235, 293

# # Generate latitude and longitude arrays at every 1 degree
# lats = np.arange(lat_min, lat_max + 1, 1)
# lons = np.arange(lon_min, lon_max + 1, 1)

# # Create a meshgrid for plotting
# lon_grid, lat_grid = np.meshgrid(lons, lats, indexing='ij')
# # Flatten the latitude, longitude, and values arrays for plotting)
# lats_flat = lat_grid.flatten()
# lons_flat = lon_grid.flatten()
# values_flat = values.flatten()  # Assuming you want to plot the first 'valid_time' slice

# # Create the Scattergeo plot
# fig = go.Figure(data=go.Scattergeo(
#     lon=lons_flat,
#     lat=lats_flat,
#     text=values_flat,
#     mode='markers',
#     marker=dict(
#         size=8,
#         color=values_flat,
#         colorscale='Viridis',
#         colorbar=dict(title='Values')
#     )
# ))

# # Update the layout to show the map
# fig.update_layout(
#     title='Scattergeo Plot',
#     geo=dict(
#         scope='north america',
#         projection_type='natural earth',
#         showland=True
#     )
# )

# fig.show()

import xarray as xr
import pandas as pd
import plotly.express as px
import numpy as np

# Load the NetCDF datafile without decoding times
file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Snow_2223.nc'
data = xr.open_dataset(file_path, decode_times=False)

times_converted = data['valid_time'].values  # Get the time values directly

# Step 3: Select the data between 2022-07-01 and 2023-07-01
start_date = "2022-07-01"
end_date = "2023-07-01"
start_datetime = np.datetime64(start_date)
end_datetime = np.datetime64(end_date)

# Slice the dataset to the desired time range for snowfall data ('sf')
sf_data = data['sf'].sel(valid_time=slice(start_datetime, end_datetime)).isel(lon=slice(50, 24, -4), lat=slice(235, 293, 4))

snowfall_at_point = sf_data / 0.1
sf_data_df = snowfall_at_point.to_dataframe().reset_index()

# Create an interactive scatter geo plot
fig = px.scatter_geo(
    sf_data_df, 
    lat='latitude', 
    lon='longitude', 
    color='sf',
    animation_frame='valid_time',
    projection="natural earth",
    title="Snowfall over time"
)
fig.update_layout(
    geo=dict(
        showland=True,
        landcolor="rgb(217, 217, 217)",
    )
)

# Save the figure to an HTML file
html_file_path = "Snowfallovertime.html"
fig.write_html(html_file_path)

# Display the saved file path (optional)
print(f"Interactive plot saved as {html_file_path}")
