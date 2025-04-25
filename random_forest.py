import xarray as xr
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt

# =============================================================================
# Step 1: Load ERA5 Data
# =============================================================================
# Assume the ERA5 data is stored in a NetCDF file named 'era5_data.nc'.
# You can download and process ERA5 data following guidelines from the ERA5 documentation:
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#heading-Parameterlistings :contentReference[oaicite:1]{index=1}
data = xr.open_dataset("era5_data.nc")

# =============================================================================
# Step 2: Feature Engineering
# =============================================================================
# For demonstration, assume our dataset includes 'temperature' and 'zonal_wind'
# variables with dimensions (time, level, lat, lon). We select data at a specific pressure level,
# e.g., 10 hPa, which can be important for stratospheric studies.
# We also focus on the polar region (e.g., latitudes 60°N to 90°N).

# Select variables at a chosen pressure level (adjust 'level' accordingly)
# and average over lat and lon.
temp = data['temperature'].sel(level=10, lat=slice(60, 90)).mean(dim=["lat", "lon"])
wind = data['zonal_wind'].sel(level=10, lat=slice(60, 90)).mean(dim=["lat", "lon"])

# Create a DataFrame from the extracted time series
df = pd.DataFrame({
    'temp': temp.values,
    'wind': wind.values
}, index=pd.to_datetime(temp.time.values))

# Optional: visualize the time series
plt.figure(figsize=(10, 4))
plt.plot(df.index, df['temp'], label='Temperature')
plt.xlabel('Time')
plt.ylabel('Temperature (K)')
plt.title('Stratospheric Temperature Time Series')
plt.legend()
plt.show()

# =============================================================================
# Step 3: Define the Target Variable (SSW Event)
# =============================================================================
# In a simplified example, we define an SSW event as occurring when the temperature exceeds a
# certain threshold. Here, we use the 95th percentile of the temperature time series as the threshold.
# In practice, SSW definitions might also incorporate wind reversal or other criteria.
threshold = df['temp'].quantile(0.95)
df['ssw_event'] = (df['temp'] > threshold).astype(int)

# Display the event counts
print("Number of SSW events:", df['ssw_event'].sum())
print("Total number of samples:", len(df))

# =============================================================================
# Step 4: Prepare Data for the Random Forest Model
# =============================================================================
# Define features (predictors) and target variable
X = df[['temp', 'wind']]
y = df['ssw_event']

# Split data into training and testing sets (using stratification due to class imbalance)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# =============================================================================
# Step 5: Build and Train the Random Forest Classifier
# =============================================================================
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# =============================================================================
# Step 6: Model Evaluation
# =============================================================================
y_pred = rf.predict(X_test)
print("\nClassification Report:")
print(classification_report(y_test, y_pred))
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))

# Optionally, display feature importances
importances = rf.feature_importances_
feature_names = X.columns
plt.figure(figsize=(6, 4))
plt.bar(feature_names, importances)
plt.xlabel('Feature')
plt.ylabel('Importance')
plt.title('Feature Importances from the Random Forest Model')
plt.show()
