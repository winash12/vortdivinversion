import numpy as np
import xarray as xr
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import metpy.calc as mpcalc
from core import BoundingBox,GridHandlerFactory, find_bounding_box_indices,bounding_box_mask
from dynamics import invert_vorticity_os21
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
def get_gfs_variable(run_date, forecast_hour=0, level=700, var_name='u', handler_type="xarray"):
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

    # 3. Open with specific filters (700 hPa, U or V, etc.)
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
run_time = datetime(2026, 3, 18, 18)

# Get U at 700 hPa
u700 = get_gfs_variable(run_time, level=700, var_name='u', handler_type="xarray")

# Get V at 700 hPa
v700 = get_gfs_variable(run_time, level=700, var_name='v', handler_type="xarray")
print(f"First Lat: {u700.latitude.values[0]}")  # Should be -90.0
print(f"Last Lat: {u700.latitude.values[-1]}")  # Should be 90.0

vort700 = mpcalc.vorticity(u700, v700)

vortex_roi = BoundingBox(
    min_lat=3.0, 
    max_lat=22.0, 
    min_lon=80.0, 
    max_lon=100.0
)
domain_roi = BoundingBox(
    min_lat=3.0, 
    max_lat=22.0, 
    min_lon=70.0, 
    max_lon=100.0
)

vortex_indices = find_bounding_box_indices(vort700, vortex_roi)
domain_indices = find_bounding_box_indices(vort700, domain_roi)

# 2. ISOLATE THE ANOMALY (The 'Source')
# Ensure bounding_box_mask returns the FULL grid size, not a crop.
# This keeps index [100, 200] in the same physical place.
vort_anomaly = bounding_box_mask(vort700, vortex_indices, fill_value=0.0)


# 3. CALCULATE METRICS
# Use the full coordinates so dx/dy match the full grid indices
dx, dy = mpcalc.lat_lon_grid_deltas(vort700.longitude, vort700.latitude)

# 4. INVOKE INVERSION
# Inside invert_vorticity_os21:
# - extract_source_data uses 'vortex_indices' to slice the anomaly
# - the loops use 'domain_indices' to calculate the wind in the plot area
starttime = time.time()
upsi, vpsi = invert_vorticity_os21(
    vortmask=vort_anomaly, 
    dx=dx, 
    dy=dy, 
    target=vortex_indices, 
    domain=domain_indices
)
print(f"Vectorization done in {time.time()-starttime:.2f}s")
# 6. QUANTIFY & ANALYZE
# Now u_psi and v_psi are the cloud-induced wind components in m/s
nd_spd_raw = np.sqrt(upsi**2 + vpsi**2)

    

valid_time = pd.to_datetime(v700.valid_time.values)
date_str = valid_time.strftime('%H UTC %d %B %Y')

# Calculate dynamic levels for the colorbar
max_val = np.ceil(float(nd_spd_raw.max()))
plot_levels = np.linspace(0, max_val, 21)

# --- 2. FIGURE SETUP ---
fig = plt.figure(figsize=(12, 10), dpi=100)
ax = plt.axes(projection=crs.PlateCarree())

nd_spd = xr.DataArray(
    data=nd_spd_raw,
    coords=vort700.coords,
    dims=vort700.dims,
    name="non_divergent_speed"
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
        'label': 'Non-divergent Wind Speed ($m\ s^{-1}$)',
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
lats = vort700.latitude.values
lons = vort700.longitude.values
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
    f"GFS 700 hPa Non-Divergent Wind Field\nValid: {date_str}", 
    fontsize=15, 
    fontweight='bold', 
    pad=25
)

# Save with a high DPI for professional quality
plt.savefig('non_divergent_wind_final.png', dpi=300, bbox_inches='tight')
plt.show()


