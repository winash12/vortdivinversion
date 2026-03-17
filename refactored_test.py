import pytest
import xarray as xr
import metpy.calc as mpcalc
import numpy as np
from .core import GridHandlerFactory, find_bounding_box_indices, get_vectorized_array_indices
from .dynamics import prepare_inversion_inputs, invert_vorticity_to_wind, invert_divergence_to_wind

def test_gfs_reconstruction_integrity():
    """
    End-to-end verification using real GFS data.
    The reconstructed wind must converge to the original GFS wind.
    """
    # 1. LOAD & NORMALIZE
    # Use a small local GFS grib/nc file for the CI/CD pipeline
    ds_raw = xr.open_dataset("tests/data/gfs_sample.nc")
    handler = GridHandlerFactory.create("xarray")
    ds = handler.normalize(ds_raw)

    # 2. PHYSICS DERIVATION
    dx, dy = mpcalc.lat_lon_grid_deltas(ds.longitude, ds.latitude)
    vort = mpcalc.vorticity(ds.ugrd, ds.vgrd)
    div = mpcalc.divergence(ds.ugrd, ds.vgrd)

    # 3. DOMAIN SUBSETTING
    min_lat, max_lat = 10.0, 20.0
    min_lon, max_lon = 200.0, 210.0
    bb = find_bounding_box_indices(ds, min_lat, max_lat, min_lon, max_lon)
    x_src, y_src = get_vectorized_array_indices(bb)

    # 4. PREPARE & INVERT
    area_vort, dx_m, dy_m = prepare_inversion_inputs(vort, dx, dy, bb)
    area_div, _, _ = prepare_inversion_inputs(div, dx, dy, bb)

    u_psi, v_psi = invert_vorticity_to_wind(area_vort, bb, x_src, y_src, dx_m, dy_m)
    u_chi, v_chi = invert_divergence_to_wind(area_div, bb, x_src, y_src, dx_m, dy_m)

    # 5. VERIFICATION (The Architect's Final Check)
    u_recon = u_psi + u_chi
    v_recon = v_psi + v_chi

    # Extract original winds for the same subset to compare
    y_s, x_s = slice(bb.y_ll, bb.y_ur + 1), slice(bb.x_ll, bb.x_ur + 1)
    u_orig = ds.ugrd.values[y_s, x_s]
    
    # Assert that the error is within a reasonable scientific tolerance (e.g., 5%)
    # Inversions aren't 100% perfect at the boundaries, but the core should match
    rmse = np.sqrt(np.mean((u_orig - u_recon)**2))
    assert rmse < 0.5  # Max allowable error in m/s
