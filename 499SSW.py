from netCDF4 import Dataset
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path

file_path = Path(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\Since_2010_Temps_1000_500.nc")
data = xr.open_dataset(file_path, decode_times=False)
print(data) # displays the entire dataset broken into  summaries for dimensions, 
            #   coordinates, variables, attributes
print(data.attrs) # displays all attributes in a dictionary
# print(data.data_vars['t'].values)
# for attribute, value in data.attrs.items():
#     print(attribute, value, '\n') # displays each attribute on a new line

# print(data.attrs['Conventions'])    # lays out attribute conventions (only given CF 
#                                     # conventions, not attribute conventions)

dimensions =  data.dims # accesses dimensions
print(dimensions) # prints dimensions as a dict
# print(dimensions['geopClimMean']) # prints indiv dimension

coords = data.coords # accesses coordinates
print(coords)

data_vars = data.data_vars # accesses all variables
print(data_vars)
# climMean = pres = data.variables['valid_time'][:] # gets an indiv variable
# print(climMean)

# climMean = data.data_vars['lsp'].attrs # gets an indiv variable vals
# print(climMean)

# mean_var_attributes = data.data_vars['latitude'].attrs # looks at a vars attrs
# print(mean_var_attributes)

# mean_var_attributes = data.data_vars['geopClimMean'].attrs['standard_name'] # gets an indiv attribute of
#                                                                             # an indiv variable
# print(mean_var_attributes)