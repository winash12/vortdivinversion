import numpy as np
import xarray as xr
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import metpy.calc as mpcalc
from core import BoundingBox,GridHandlerFactory, find_bounding_box_indices,bounding_box_mask
from dynamics import invert_vorticity_os21,invert_divergence_os21
import sys
import os
from  datetime import datetime
import time
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.cm import get_cmap

def get_gfs_variable(run_date, forecast_hour=0, level=200, var_name='u', handler_type="xarray"):
    """
    Downloads and normalizes a SINGLE GFS variable using the Factory pattern.
    """
    # 1. Build Dynamic Paths
    date_str = run_date.strftime("%Y%m%d")
    cycle_str = run_date.strftime("%H")
    f_hr = f"{forecast_hour:03d}"
    
    s3_key = f"gfs.{date_str}/{cycle_str}/atmos/gfs.t{cycle_str}z.pgrb2.0p25.f{f_hr}"
    local_filename = f"gfs_{date_str}_{cycle_str}z_f{f_hr}.grib2"

    # 2. Download (Boto3)
    if not os.path.exists(local_filename):
        client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        client.download_file('noaa-gfs-bdp-pds', s3_key, local_filename)

    # 3. Open with specific filters (200 hPa, U or V, etc.)
    ds_raw = xr.open_dataset(
        local_filename, 
        engine='cfgrib', 
        backend_kwargs={'filter_by_keys': {
            'typeOfLevel': 'isobaricInhPa', 
            'level': level, 
            'shortName': var_name
        }}
    )

    # 4. Use Factory to create the handler (Passes the string 'xarray' or 'numpy')
    handler = GridHandlerFactory.create(handler_type)

    # 5. Normalize and return based on the handler type
    if handler_type == "xarray":
        # Returns a normalized Xarray DataArray
        ds_normalized = handler.normalize(ds_raw)
        return ds_normalized[var_name]
    
    elif handler_type == "numpy":
        # Returns (lats, lons, field_array)
        return handler.normalize(
            ds_raw.latitude.values, 
            ds_raw.longitude.values, 
            ds_raw[var_name].values
        )

# --- Usage ---
run_time = datetime(2025, 12, 1, 6)

# Get U at 200 hPa
u200 = get_gfs_variable(run_time, level=200, var_name='u', handler_type="xarray")

# Get V at 200 hPa
v200 = get_gfs_variable(run_time, level=200, var_name='v', handler_type="xarray")
print(f"First Lat: {u200.latitude.values[0]}")  # Should be -90.0
print(f"Last Lat: {u200.latitude.values[-1]}")  # Should be 90.0

div200 = mpcalc.divergence(u200, v200)

# 1. Expand the 'Source' box (Vortex ROI)
# This captures the full vorticity/divergence 'charge' to avoid the 22N clipping.
vortex_roi = BoundingBox(
    min_lat=3.0, 
    max_lat=35.0,  # Expanded from 22.0
    min_lon=75.0,  # Slightly wider to capture the Western Ghats influence
    max_lon=95.0   # Captures the full Bay of Bengal outflow
)

# 2. Expand the 'Plotting' box (Domain ROI)
# This ensures the wind field is calculated for the entire map view.
domain_roi = BoundingBox(
    min_lat=3.0, 
    max_lat=36.0,  # Expanded to match your map's top edge
    min_lon=65.0,  # Shows more of the Arabian Sea for context
    max_lon=100.0  # Shows the full Andaman Sea/Myanmar coast
)

vortex_indices = find_bounding_box_indices(div200, vortex_roi)
domain_indices = find_bounding_box_indices(div200, domain_roi)

# 2. ISOLATE THE ANOMALY (The 'Source')
# Ensure bounding_box_mask returns the FULL grid size, not a crop.
# This keeps index [100, 200] in the same physical place.
div_anomaly = bounding_box_mask(div200, vortex_indices, fill_value=0.0)


# 3. CALCULATE METRICS
# Use the full coordinates so dx/dy match the full grid indices
dx, dy = mpcalc.lat_lon_grid_deltas(div200.longitude, div200.latitude)

# 4. INVOKE INVERSION
# Inside invert_vorticity_os21:
# - extract_source_data uses 'vortex_indices' to slice the anomaly
# - the loops use 'domain_indices' to calculate the wind in the plot area
starttime = time.time()
upsi, vpsi = invert_divergence_os21(
    divmask=div_anomaly, 
    dx=dx, 
    dy=dy, 
    target=vortex_indices, 
    domain=domain_indices
)
print(f"Vectorization done in {time.time()-starttime:.2f}s")
# 6. QUANTIFY & ANALYZE
# Now u_psi and v_psi are the cloud-induced wind components in m/s
nd_spd_raw = np.sqrt(upsi**2 + vpsi**2)

valid_time = pd.to_datetime(v200.valid_time.values)
date_str = valid_time.strftime('%H UTC %d %B %Y')

# Calculate dynamic levels for the colorbar
max_val = np.ceil(float(nd_spd_raw.max()))
plot_levels = np.linspace(0, max_val, 21)

# --- 2. FIGURE SETUP ---
fig = plt.figure(figsize=(12, 10), dpi=100)
ax = plt.axes(projection=crs.PlateCarree())

nd_spd = xr.DataArray(
    data=nd_spd_raw,
    coords=div200.coords,
    dims=div200.dims,
    name="irrotational_speed"
)

# --- 3. BACKGROUND PLOT (Wind Magnitude) ---
# We use xarray's plot functionality with our dynamic levels
plot = nd_spd.plot(
    ax=ax,
    levels=plot_levels,
    cmap='YlGnBu',
    transform=crs.PlateCarree(),
    add_colorbar=True,
    cbar_kwargs={
        'label': 'Irrotational Wind Speed ($m\ s^{-1}$)',
        'shrink': 0.8,
        'pad': 0.05,
        'aspect': 30,
        'orientation': 'vertical'
    }
)

# --- 4. OVERLAY WIND VECTORS (Quivers) ---
# We skip every 6th data point so the arrows don't overlap too much
skip = 6
# Slicing the 2D arrays: [y_slice, x_slice]
sub_u = upsi[::skip, ::skip]
sub_v = vpsi[::skip, ::skip]
lats = div200.latitude.values
lons = div200.longitude.values
sub_lons = lons[::skip]
sub_lats = lats[::skip]

# Add the arrows
q = ax.quiver(
    sub_lons, sub_lats, sub_u, sub_v,
    transform=crs.PlateCarree(),
    color='black',
    alpha=0.6,    # Slightly transparent so we can see the colors underneath
    scale=150,    # Adjust this if arrows look too long or short
    width=0.002,
    headwidth=3
)

# Add a legend for the arrows (Top Right)
ax.quiverkey(q, X=0.9, Y=1.03, U=10, label='10 $m\ s^{-1}$', labelpos='E')

# --- 5. GEOGRAPHIC FEATURES ---
# Standard coastlines
ax.coastlines('50m', linewidth=1.0, color='black')

# Country borders (Cultural features)
countries = NaturalEarthFeature(
    category="cultural", scale="50m", 
    facecolor="none", name="admin_0_countries"
)
ax.add_feature(countries, linewidth=0.6, edgecolor="#444444")

# Set the map view to your specific region
ax.set_extent([70, 100., 6., 36.], crs=crs.PlateCarree())

# --- 6. GRIDLINES & LABELS ---
gl = ax.gridlines(
    draw_labels=True, 
    color="grey", 
    linestyle="--", 
    alpha=0.4, 
    linewidth=0.5
)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 11, 'color': 'black'}
gl.ylabel_style = {'size': 11, 'color': 'black'}

# --- 7. DYNAMIC TITLE & OUTPUT ---
plt.title(
    f"GFS 200 hPa Divergent Wind Field\nValid: {date_str}", 
    fontsize=15, 
    fontweight='bold', 
    pad=25
)

# Save with a high DPI for professional quality
plt.savefig('divergent_wind_final.png', dpi=300, bbox_inches='tight')
plt.show()


