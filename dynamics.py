import numpy as np

# Import the architectural components from your core.py
# We use BoundingBoxIndices for type-hinting the arguments
from core import BoundingBoxIndices,extract_source_data,extract_source_data_ms


def invert_vorticity_os21(vortmask, dx, dy, target, domain):
    # 1. Helper to strip units/metadata (The "Agnostic Utility")
    def to_raw(obj):
        if hasattr(obj, "magnitude"): return obj.magnitude # Strips MetPy units
        if hasattr(obj, "values"): return obj.values       # Strips Xarray metadata
        return obj

    # 2. Convert inputs to raw NumPy arrays immediately
    vort_raw = to_raw(vortmask)
    dx_raw = to_raw(dx)
    dy_raw = to_raw(dy)

    # 3. Pass the RAW arrays to your extractor
    #area_vort, x_src, y_src, dx_f, dy_f = extract_source_data(vort_raw, dx_raw, dy_raw,
    #                                                          target)
    area_vort, x_src, y_src, dx_f, dy_f = extract_source_data_ms(vort_raw, dx_raw, dy_raw,
                                                                 target,ds_meta=vortmask)
    # 4. Pre-allocate results using the raw shape
    upsi = np.zeros_like(vort_raw)
    vpsi = np.zeros_like(vort_raw)

    # 5. Fast Vectorized Loop
    for i in range(domain.x_ll, domain.x_ur):
        for j in range(domain.y_ll, domain.y_ur):
            
            xdiff = (i - x_src) * dx_f
            ydiff = (j - y_src) * dy_f
            rsq = xdiff**2 + ydiff**2
            
            # This comparison now works because rsq is raw float, not Pint Quantity
            mask = (rsq > 1e-6) 
            inv_rsq = 1.0 / rsq[mask]
            
            upsi[j, i] = np.sum(area_vort[mask] * -ydiff[mask] * inv_rsq)
            vpsi[j, i] = np.sum(area_vort[mask] * xdiff[mask] * inv_rsq)

    return upsi / (2 * np.pi), vpsi / (2 * np.pi)

def invert_divergence_os21(divmask, dx, dy, target, domain):
    """
    Performs OS2021 Green's Function Inversion for Divergence (Eq 10 & 11).
    Agnostic to Xarray/NumPy and unit-safe.
    """
    # 1. Helper to strip units/metadata (Same as we used for Vorticity)
    def to_raw(obj):
        if hasattr(obj, "magnitude"): return obj.magnitude 
        if hasattr(obj, "values"): return obj.values       
        return obj

    # 2. Pre-process inputs to raw NumPy
    div_raw = to_raw(divmask)
    dx_raw = to_raw(dx)
    dy_raw = to_raw(dy)

    # 3. Get pre-calculated strengths (using Divergence)
    area_div, x_src, y_src, dx_f, dy_f = extract_source_data(div_raw, dx_raw, dy_raw, target)

    area_vort, x_src, y_src, dx_f, dy_f = extract_source_data_ms(div_raw, dx_raw, dy_raw,
                                                                 target,ds_meta=divmask)

    
    # 4. Pre-allocate results using the raw shape
    uchi = np.zeros_like(div_raw)
    vchi = np.zeros_like(div_raw)

    # 5. Vectorized Calculation Loop
    for i in range(domain.x_ll, domain.x_ur):
        for j in range(domain.y_ll, domain.y_ur):
            
            # Displacement vectors (Target point [i,j] - Source meshgrid [x_src, y_src])
            xdiff = (i - x_src) * dx_f
            ydiff = (j - y_src) * dy_f
            rsq = xdiff**2 + ydiff**2
            
            # Singularity mask
            mask = (rsq > 1e-6)
            inv_rsq = 1.0 / rsq[mask]
            
            # --- Oertel & Schemm 2021 Eq 10 & 11 ---
            # Irrotational wind (chi) uses POSITIVE signs for both components:
            # u_chi = d(chi)/dx,  v_chi = d(chi)/dy
            uchi[j, i] = np.sum(area_div[mask] * xdiff[mask] * inv_rsq)
            vchi[j, i] = np.sum(area_div[mask] * ydiff[mask] * inv_rsq)

    scaling = 1 / (2 * np.pi)
    return uchi * scaling, vchi * scaling

def reconstruct_total_induced_wind(vortmask, divmask, dx, dy, target, domain):
    """
    Coordinates the dual inversion (Vorticity & Divergence) to get the 
    total induced wind field within the calculation domain.
    """
    # 1. Perform Rotational Inversion (Eq 8 & 9)
    u_psi, v_psi = invert_vorticity_os21(vortmask, dx, dy, target, domain)
    
    # 2. Perform Divergent Inversion (Eq 10 & 11)
    u_chi, v_chi = invert_divergence_os21(divmask, dx, dy, target, domain)
    
    # 3. Sum the components for the Total Induced Wind
    u_total = u_psi + u_chi
    v_total = v_psi + v_chi
    
    # 4. Wrap back into Xarray for plotting/saving
    # We use the original coordinates from the vortmask
    results = xr.Dataset(
        data_vars={
            "u_psi": (("latitude", "longitude"), u_psi),
            "v_psi": (("latitude", "longitude"), v_psi),
            "u_chi": (("latitude", "longitude"), u_chi),
            "v_chi": (("latitude", "longitude"), v_chi),
            "u_total": (("latitude", "longitude"), u_total),
            "v_total": (("latitude", "longitude"), v_total),
        },
        coords=vortmask.coords
    )
    
    # Slice the final dataset to the domain for clean output
    return results.isel(
        latitude=slice(domain.y_ll, domain.y_ur),
        longitude=slice(domain.x_ll, domain.x_ur)
    )
