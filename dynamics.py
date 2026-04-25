import numpy as np

# Import the architectural components from your core.py
# We use BoundingBoxIndices for type-hinting the arguments
from core import BoundingBoxIndices,extract_source_data,extract_source_data_ms,coriolis_parameter


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


    for i in range(domain.x_ll, domain.x_ur):
        for j in range(domain.y_ll, domain.y_ur):

            xdiff = (i - x_src) * dx_f
            ydiff = (j - y_src) * dy_f
            rsq = xdiff**2 + ydiff**2
            # This comparison now works because rsq is raw float, not Pint Quantity

            mask = (rsq > 1e-6) 
            inv_rsq = 1.0 / rsq[mask]

            upsi[j, i] = -np.sum(area_vort[mask] * ydiff[mask] * inv_rsq)
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
    #area_div, x_src, y_src, dx_f, dy_f = extract_source_data(div_raw, dx_raw, dy_raw, target)

    area_div, x_src, y_src, dx_f, dy_f = extract_source_data_ms(div_raw, dx_raw, dy_raw,
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
    return utotal,vtotal

def invert_scalars_os21(field_mask, dx, dy, target, domain, scalar_type="streamfunction"):
    """
    Performs OS2021 Green's Function Inversion for Scalars (Eq 6 & 7).
    Equation 6: Streamfunction (psi) from Vorticity
    Equation 7: Velocity Potential (chi) from Divergence
    """
    f_raw = field_mask.values if hasattr(field_mask, "values") else field_mask
    dx_raw = dx.values if hasattr(dx, "values") else dx
    dy_raw = dy.values if hasattr(dy, "values") else dy

    area_f, x_src, y_src, dx_f, dy_f = extract_source_data_ms(f_raw, dx_raw, dy_raw, 
                                                              target, ds_meta=field_mask)
    
    # Pre-allocate scalar field
    scalar_field = np.zeros_like(f_raw)

    for i in range(domain.x_ll, domain.x_ur):
        for j in range(domain.y_ll, domain.y_ur):
            xdiff = (i - x_src) * dx_f
            ydiff = (j - y_src) * dy_f
            rsq = xdiff**2 + ydiff**2
            mask = (rsq > 1e-6)
            # --- Oertel & Schemm 2021 Eq 12 & 13 ---
            # Both psi and chi use the log of the distance for the 2D Green's function
            # G = 1/(2*pi) * ln(r)
            # Result = sum( Strength * G )
            scalar_field[j, i] = np.sum(area_f[mask] * 0.5 * np.log(rsq[mask]))

    scaling = 1 / (2 * np.pi)
    return scalar_field * scaling

def reconstruct_total_induced_wind(vortmask, divmask, dx, dy, target, domain, u_obs=None, v_obs=None):
    """
    Coordinates the triple-decomposition: Rotational, Divergent, and Environmental.
    Reference: Oertel & Schemm (2021) Eq 8-12.
    """
    # 1. Induced Rotational Wind (Eq 8 & 9)
    u_psi, v_psi = invert_vorticity_os21(vortmask, dx, dy, target, domain)
    
    # 2. Induced Divergent Wind (Eq 10 & 11)
    u_chi, v_chi = invert_divergence_os21(divmask, dx, dy, target, domain)
    
    # 3. Sum local induced components
    u_induced = u_psi + u_chi
    v_induced = v_psi + v_chi
    
    # 4. Calculate Environmental/Harmonic Wind (Eq 12)
    # If observed total winds are provided (e.g. from GFS/WRF)
    u_env, v_env = None, None
    if u_obs is not None and v_obs is not None:
        # Subtract induced components from the total observed field
        # Ensure indices match the calculation 'domain'
        u_raw = u_obs.values if hasattr(u_obs, "values") else u_obs
        v_raw = v_obs.values if hasattr(v_obs, "values") else v_obs
        
        # Slicing the observed wind to match the domain result
        u_env = u_raw[domain.y_ll:domain.y_ur, domain.x_ll:domain.x_ur] - u_induced
        v_env = v_raw[domain.y_ll:domain.y_ur, domain.x_ll:domain.x_ur] - v_induced

    return u_induced, v_induced, u_env, v_env

def calculate_full_vertical_vorticity_tendency_ms(vortmask, divmask, dx, dy, target, 
                                                  u_obs, v_obs, lats_1d, ds_meta=None):
    """
    Calculates the full horizontal vorticity tendency (Eq 13) for GFS data.
    """
    
    def to_raw(obj):
        if hasattr(obj, "magnitude"): return obj.magnitude
        if hasattr(obj, "values"): return obj.values
        return obj

    vort_raw = to_raw(vortmask)
    div_raw = to_raw(divmask)
    u_raw = to_raw(u_obs)
    v_raw = to_raw(v_obs)
    dx_raw = to_raw(dx)
    dy_raw = to_raw(dy)
    
    # 1. Slice horizontal wind and mask arrays according to target domain
    u_full = u_raw[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    v_full = v_raw[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    vort_domain = vort_raw[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    div_domain = div_raw[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    
    dx_domain = dx_raw[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    dy_domain = dy_raw[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    
    # 2. CALCULATE CUSTOM CORIOLIS PARAMETER (f)
    # Slice the 1D GFS latitudes matching the Y-axis calculation domain
    lats_sliced = lats_1d[target.y_ll:target.y_ur]
    
    # Use our unit-free custom function
    f_1d = coriolis_parameter(lats_sliced)
    
    # Broadcast the 1D 'f' vector into a matching 2D array grid
    _, f_grid = np.meshgrid(np.arange(target.x_ll, target.x_ur), f_1d)
    
    # 3. GET MAP FACTORS
    meta_source = ds_meta if ds_meta is not None else vortmask
    mx, my = get_projection_scaling(meta_source, target)
    
    # 4. Compute Cartesian gradients of vorticity
    vort_grad_y_raw, vort_grad_x_raw = np.gradient(vort_domain, dy_domain, dx_domain)
    vort_grad_x = vort_grad_x_raw / mx
    vort_grad_y = vort_grad_y_raw / my
    
    # 5. Component A: Horizontal Advection of Vorticity
    advection = -(u_full * vort_grad_x + v_full * vort_grad_y)
    
    # 6. Component B: Vortex Stretching (Now dynamically accurate per latitude)
    absolute_vorticity = vort_domain + f_grid
    stretching = -(absolute_vorticity * div_domain)
    
    total_tendency = advection + stretching
    
    return total_tendency
