import xarray as xr
import os
# Define the paths to the directories containing NetCDF files
folder1 = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\globsnow_swe_data_2022"
folder2 = r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\globsnow_swe_data_2023"
# Get all NetCDF files in both directories
files_folder1 = [os.path.join(folder1, f) for f in os.listdir(folder1) if f.endswith(".nc.gz")]
files_folder2 = [os.path.join(folder2, f) for f in os.listdir(folder2) if f.endswith(".nc.gz")]

# Combine both lists
all_files = sorted(files_folder1 + files_folder2)
# Use xarray to open and concatenate all files along the time dimension
# The combine='by_coords' argument will concatenate along existing coordinates (like time)
combined_dataset = xr.open_mfdataset(all_files, combine='by_coords')
# Sort the combined dataset by the time dimension
sorted_dataset = combined_dataset.sortby('data_date')
# Save the combined and sorted dataset to a new NetCDF file
output_file = "combined_sorted_globsnow22-23.nc"
sorted_dataset.to_netcdf(output_file)
# Load the combined file to verify its contents
combined_data = xr.open_dataset(output_file)
print(combined_data)
