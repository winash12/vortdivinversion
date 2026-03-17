import numpy as np

# Import the architectural components from your core.py
# We use BoundingBoxIndices for type-hinting the arguments
from .core import BoundingBoxIndices

def prepare_inversion_inputs(masked_field, dx, dy, bb):
    # Slice the already-masked field using the indices found in step 4
    y_slice = slice(bb.y_ll, bb.y_ur + 1)
    x_slice = slice(bb.x_ll, bb.x_ur + 1)
    
    f_vals = masked_field.values[y_slice, x_slice]
    dx_vals = dx.magnitude[y_slice, x_slice]
    dy_vals = dy.magnitude[y_slice, x_slice]
    
    return f_vals * dx_vals * dy_vals, np.mean(dx_vals), np.mean(dy_vals)

def invert_vorticity_to_wind(vort_field, bb_indices, x_indices, y_indices, dx, dy):
    """
    Performs the Biot-Savart inversion of relative vorticity into rotational wind.
    
    Implements Oertel & Schemm (2021) Equations 8 & 9.
    Assumes input data is normalized (South-to-North) via GridHandler.
    
    Parameters
    ----------
    vort_field : ndarray
        The 2D area-weighted vorticity field (zeta * dA).
    bb_indices : BoundingBoxIndices
        The indices of the target domain subset.
    x_indices, y_indices : ndarray
        The 2D meshgrids representing the source grid (from core.py).
    dx, dy : float
        Grid spacing in meters.
        
    Returns
    -------
    tuple : (u_psi, v_psi) 2D arrays of rotational wind components.
    """
    # Initialize output arrays with the same shape as the target subset
    # We use a subset-sized array to save memory
    shape = (bb_indices.y_ur - bb_indices.y_ll + 1, 
             bb_indices.x_ur - bb_indices.x_ll + 1)
    
    u_psi = np.zeros(shape)
    v_psi = np.zeros(shape)

    # The Nested Inversion (The "Real Deal")
    # i and j are local indices for the output arrays
    for i_idx, i_glob in enumerate(range(bb_indices.x_ll, bb_indices.x_ur + 1)):
        for j_idx, j_glob in enumerate(range(bb_indices.y_ll, bb_indices.y_ur + 1)):
            
            # Calculate distance: (Target Global Index - Source Mesh) * Spacing
            xdiff = (i_glob - x_indices) * dx
            ydiff = (j_glob - y_indices) * dy
            rsq = xdiff**2 + ydiff**2
            
            # Singularity mask
            mask = (rsq > 1e-6)
            inv_rsq = 1.0 / rsq[mask]
            
            # Eq 8: u_psi = -sum(vort * dy / r^2)
            u_psi[j_idx, i_idx] = -np.sum(vort_field[mask] * ydiff[mask] * inv_rsq)
            
            # Eq 9: v_psi = sum(vort * dx / r^2)
            v_psi[j_idx, i_idx] =  np.sum(vort_field[mask] * xdiff[mask] * inv_rsq)

    # Global Scaling (Efficient Vectorization)
    scale = 1.0 / (2.0 * np.pi)
    return u_psi * scale, v_psi * scale

def invert_divergence_to_wind(div_field, bb_indices, x_indices, y_indices, dx, dy):
    """
    Performs the Biot-Savart inversion of divergence into divergent wind.
    
    Implements Oertel & Schemm (2021) Equations 10 & 11.
    Assumes input data is normalized (South-to-North) via GridHandler.
    
    Parameters
    ----------
    div_field : ndarray
        The 2D area-weighted divergence field (delta * dA).
    bb_indices : BoundingBoxIndices
        The indices of the target domain subset.
    x_indices, y_indices : ndarray
        The 2D meshgrids representing the source grid (from core.py).
    dx, dy : float
        Grid spacing in meters.
        
    Returns
    -------
    tuple : (u_chi, v_chi) 2D arrays of divergent wind components.
    """
    # Initialize output arrays (matching subset dimensions)
    shape = (bb_indices.y_ur - bb_indices.y_ll + 1, 
             bb_indices.x_ur - bb_indices.x_ll + 1)
    
    u_chi = np.zeros(shape)
    v_chi = np.zeros(shape)

    # i_glob/j_glob are global indices for distance; i_idx/j_idx for local array storage
    for i_idx, i_glob in enumerate(range(bb_indices.x_ll, bb_indices.x_ur + 1)):
        for j_idx, j_glob in enumerate(range(bb_indices.y_ll, bb_indices.y_ur + 1)):
            
            # Calculate distance (Target Global Index - Source Mesh)
            xdiff = (i_glob - x_indices) * dx
            ydiff = (j_glob - y_indices) * dy
            rsq = xdiff**2 + ydiff**2
            
            mask = (rsq > 1e-6)
            inv_rsq = 1.0 / rsq[mask]
            
            # EQ 10: u_chi = sum(div * dx / r^2)
            u_chi[j_idx, i_idx] = np.sum(div_field[mask] * xdiff[mask] * inv_rsq)
            
            # EQ 11: v_chi = sum(div * dy / r^2)
            v_chi[j_idx, i_idx] = np.sum(div_field[mask] * ydiff[mask] * inv_rsq)

    # Global Scaling (Efficient Vectorization)
    scale = 1.0 / (2.0 * np.pi)
    return u_chi * scale, v_chi * scale
