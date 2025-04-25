import rasterio
import matplotlib.pyplot as plt
import numpy as np

# Define the path to your .byte file
byte_file_path = r'C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\llipa.byte'

# Open the .byte file using Rasterio
with rasterio.open(byte_file_path) as dataset:
    # Read the data into a numpy array (first band)
    data = dataset.read(1)
    
    # Determine the pixel size in the y-direction.
    # For a north-up image, dataset.transform.e is negative, so we take its absolute value.
    pixel_size_y = abs(dataset.transform.e)
    
    # Create an array of row indices (each row corresponds to a latitude)
    nrows = data.shape[0]
    rows = np.arange(nrows)
    
    # Compute the center latitude for each row.
    # For a north-up raster, the top boundary is dataset.bounds.top.
    center_lats = dataset.bounds.top - (rows + 0.5) * pixel_size_y
    
    # Find the indices of rows with center latitude above 25°
    valid_rows = np.where(center_lats > 25)[0]
    if valid_rows.size == 0:
        raise ValueError("No rows found with latitude above 25°.")
    
    # Because the rows are ordered from north to south, these valid rows are contiguous from the top.
    # We'll take all rows from the top down to the last row that meets the condition.
    last_valid_row = valid_rows[-1]
    data_subset = data[:last_valid_row + 1, :]
    
    # Update the plotting extent so that the axes show the correct longitudes and latitudes.
    # The full dataset's horizontal extent remains the same.
    new_left = dataset.bounds.left
    new_right = dataset.bounds.right
    new_top = dataset.bounds.top
    # The new bottom is the bottom edge of the last valid row.
    # For row index i, the bottom edge is: top - (i+1)*pixel_size_y.
    new_bottom = dataset.bounds.top - (last_valid_row + 1) * pixel_size_y
    new_extent = (new_left, new_right, new_bottom, new_top)
    
    # Plot the subset of the data
    plt.figure(figsize=(10, 8))
    plt.imshow(data_subset, cmap='viridis', extent=new_extent, origin='upper')
    plt.colorbar(label='Permafrost Classification')
    plt.title('Circum-Arctic Map of Permafrost and Ground-Ice Conditions (Above 25° Latitude)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()
