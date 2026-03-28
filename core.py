import numpy as np
import xarray as xr
from pyproj import Proj,CRS,Geod
from dataclasses import dataclass
from abc import ABC, abstractmethod

@dataclass
class GridMetadata:
    """Metadata to handle map projections without hardcoding."""
    dx_metric: np.ndarray  # True distance in meters
    dy_metric: np.ndarray
    map_factor_x: np.ndarray = None
    map_factor_y: np.ndarray = None

@dataclass
class BoundingBox:
    """Geographic coordinates for a region of interest."""
    min_lat: float
    max_lat: float
    min_lon: float
    max_lon: float
    
@dataclass
class BoundingBoxIndices:
    """Container for array slice bounds to ensure type-safety in the inversion."""
    x_ll: int
    x_ur: int
    y_ll: int
    y_ur: int


class GridHandler(ABC):
    """Abstract Base Class defining the interface for Grid Normalization."""
    @abstractmethod
    def normalize(self, *args, **kwargs):
        pass


class XarrayGridHandler(GridHandler):
    """Strategy for Xarray objects: Leverages coordinate-aware sorting."""
    def normalize(self, obj):
        # 1. Flip Latitude to South-to-North (Ascending)
        # Use .values[0] for explicit scalar comparison
        if obj.latitude.values[0] > obj.latitude.values[-1]:
            obj = obj.sortby("latitude", ascending=True)
            
        # 2. Normalize Longitude to 0-360 and ensure it is monotonic
        obj = obj.assign_coords(longitude=(obj.longitude % 360))
        return obj.sortby("longitude", ascending=True)


class NumpyGridHandler(GridHandler):
    """Strategy for raw NumPy arrays: Manually synchronizes lats and data."""
    def normalize(self, lats, lons, *fields):
        # 1. Flip Lats/Fields to South-to-North if descending
        if lats[0] > lats[-1]:
            lats = lats[::-1]
            fields = [np.flip(f, axis=0) for f in fields]
            
        # 2. Normalize Longitude to 0-360
        lons = lons % 360
        
        # 3. Ensure monotonic longitude (West-to-East)
        if not np.all(np.diff(lons) > 0):
            idx = np.argsort(lons)
            lons = lons[idx]
            # Use ellipsis to handle 2D or 3D field structures
            fields = [f[..., idx] for f in fields]
            
        return lats, lons, *fields


class GridHandlerFactory:
    """Dynamic Factory to instantiate the correct GridHandler by name."""
    _handlers = {
        "xarray": XarrayGridHandler,
        "numpy": NumpyGridHandler
    }

    @staticmethod
    def create(handler_type: str) -> GridHandler:
        handler_class = GridHandlerFactory._handlers.get(handler_type.lower())
        if not handler_class:
            raise ValueError(f"Unknown handler type: {handler_type}. Use 'xarray' or 'numpy'.")
        return handler_class()


def bounding_box_mask(data, indices: BoundingBoxIndices, fill_value=0.0):
    """
    TRULY AGNOSTIC: Handles Xarray DataArrays or NumPy ndarrays.
    Ensures the vortex is pasted into a full grid of zeros.
    """
    # 1. Create the 'Empty' Background Grid
    # xr.zeros_like preserves coordinates/metadata for publication plots
    if hasattr(data, "zeros_like"):
        mask_result = xr.zeros_like(data)
    else:
        mask_result = np.full_like(data, fill_value)

    # 2. Define the Slice for Readability
    # Slicing is [South_Index : North_Index, West_Index : East_Index]
    # Indices are guaranteed Small-to-Large by your find_bounding_box_indices
    y_slice = slice(indices.y_ll, indices.y_ur)
    x_slice = slice(indices.x_ll, indices.x_ur)

    # 3. 'Paste' the Data
    # We use .values for Xarray to ensure we are modifying the underlying NumPy buffer
    if hasattr(mask_result, "values"):
        mask_result.values[..., y_slice, x_slice] = data.values[..., y_slice, x_slice]
    else:
        mask_result[..., y_slice, x_slice] = data[..., y_slice, x_slice]
            
    return mask_result


def find_bounding_box_indices(data, bbox: BoundingBox) -> BoundingBoxIndices:
    """
    TRULY AGNOSTIC: Handles xarray objects OR a tuple of (lats, lons) numpy arrays.
    """
    # 1. Extract Coordinate Values
    if hasattr(data, "latitude") and hasattr(data, "longitude"):
        # Case: Xarray Object
        lats = data.latitude.values
        lons = data.longitude.values
    elif isinstance(data, (tuple, list)) and len(data) >= 2:
        # Case: NumPy Tuple (lats, lons)
        lats, lons = data[0], data[1]
    else:
        raise TypeError("Input must be an Xarray object or a tuple of (lats, lons) arrays.")

    # 2. Find closest index (Nearest Neighbour)
    # This works regardless of whether lats are [90...0] or [0...90]
    def get_idx(arr, val):
        return int(np.argmin(np.abs(arr - val)))

    iy1, iy2 = get_idx(lats, bbox.min_lat), get_idx(lats, bbox.max_lat)
    ix1, ix2 = get_idx(lons, bbox.min_lon), get_idx(lons, bbox.max_lon)

    # 3. Force Small-to-Large Slicing
    # This fixes the (0,0) shape error and the Evans N-S assumption
    return BoundingBoxIndices(
        y_ll = min(iy1, iy2),
        y_ur = max(iy1, iy2) + 1,
        x_ll = min(ix1, ix2),
        x_ur = max(ix1, ix2) + 1
    )

def get_projection_scaling(ds, indices: BoundingBoxIndices):
    """
    Independent Projection Handler.
    Detects CRS from CF-attributes (no MetPy) and calculates scale factors.
    """
    ny = indices.y_ur - indices.y_ll
    nx = indices.x_ur - indices.x_ll

    # 1. Try to find CF-Convention grid_mapping
    # Most professional datasets use a variable to define the CRS
    crs_var = None
    for var in ds.variables:
        if 'grid_mapping_name' in ds[var].attrs:
            crs_var = ds[var]
            break

    if crs_var is not None:
        # Build the Proj object directly from CF attributes
        cf_params = crs_var.attrs
        crs = CRS.from_cf(cf_params)
        p = Proj(crs)
        
        lons = ds.longitude.values[indices.y_ll:indices.y_ur, indices.x_ll:indices.x_ur]
        lats = ds.latitude.values[indices.y_ll:indices.y_ur, indices.x_ll:indices.x_ur]
        
        # Calculate Parallel (k) and Meridional (h) scale factors
        # These account for Rotated Pole (COSMO) or Lambert (WRF)
        factors = p.get_factors(lons, lats)
        return factors.parallel_scale, factors.meridional_scale

    # 2. Fallback for Standard Lat/Lon (GFS / ERA5)
    # If no CRS is found, we assume Equidistant Cylindrical
    lats = ds.latitude.values[indices.y_ll:indices.y_ur]
    mx_1d = 1.0 / np.cos(np.deg2rad(lats))
    mx = np.repeat(mx_1d[:, np.newaxis], nx, axis=1)
    my = np.ones_like(mx)
    
    return mx, my

def extract_source_data(vortmask, dx, dy, target: BoundingBoxIndices):
    """
    Final Publication Version:
    1. Agnostic to input type (Xarray/NumPy).
    2. Guaranteed South-to-North ordering.
    3. Pre-calculates 'area_vort' for Green's Summation.
    """
    # Agnostic extraction of raw numbers
    v_raw = vortmask.values if hasattr(vortmask, "values") else vortmask
    dx_raw = dx.values if hasattr(dx, "values") else dx
    dy_raw = dy.values if hasattr(dy, "values") else dy

    # Slice using the normalized South-to-North indices
    vort_f = v_raw[target.y_ll : target.y_ur, target.x_ll : target.x_ur]
    dx_f = dx_raw[target.y_ll : target.y_ur, target.x_ll : target.x_ur]
    dy_f = dy_raw[target.y_ll : target.y_ur, target.x_ll : target.x_ur]

    # Pre-calculate source strength (The 'Charge' in Green's function)
    area_vort = vort_f * dx_f * dy_f

    # Create matching coordinate grids for the inversion loop
    source_y = np.arange(target.y_ll, target.y_ur)
    source_x = np.arange(target.x_ll, target.x_ur)
    x_src, y_src = np.meshgrid(source_x, source_y)

    return area_vort, x_src, y_src, dx_f, dy_f

def extract_source_data_ms(field_mask, dx, dy, target: BoundingBoxIndices):
    """
    JOSS-compliant extractor: Uses Map Factors to ensure 
    Area Strength is projection-independent.
    """
    # 1. Get raw slices
    field_f = field_mask.values[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    dx_f = dx.magnitude[target.y_ll:target.y_ur, target.x_ll:target.x_ur]
    dy_f = dy.magnitude[target.y_ll:target.y_ur, target.x_ll:target.x_ur]

    # 2. Get Map Factors (No hardcoding)
    mx, my = get_projection_scaling(field_mask, target)

    # 3. AREA STRENGTH (The 'Physical Charge')
    # We divide by the map factors to convert 'Map Area' to 'True Earth Area'
    area_strength = field_f * (dx_f / mx) * (dy_f / my)

    # ... [Coordinate grid logic] ...
    return area_strength, x_src, y_src, dx_f, dy_f

def get_geod_from_ds(ds):
    """
    Extracts the specific ellipsoid/sphere from the dataset's CRS.
    If no CRS is found, defaults to the GFS sphere (R=6371229).
    """
    # 1. Search for CF-compliant grid mapping
    crs_var = None
    for var in ds.variables:
        if 'grid_mapping_name' in ds[var].attrs:
            crs_var = ds[var]
            break
            
    if crs_var is not None:
        # Build CRS from CF attributes (WRF/COSMO)
        crs = CRS.from_cf(crs_var.attrs)
    else:
        # Default to GFS Sphere (The standard for global models)
        crs = CRS.from_dict({'proj': 'longlat', 'a': 6371229, 'b': 6371229})
        
    # Return the specific Geodesic object for this Earth shape
    return crs.get_geod()

def calculate_metric_deltas(lons, lats, geod: Geod):
    """
    Uses the model-specific Geod to find dx and dy.
    lons/lats: 2D arrays of coordinates.
    """
    # Calculate dy (North-South) - inv(lon1, lat1, lon2, lat2)
    _, _, dy = geod.inv(lons[:-1, :], lats[:-1, :], 
                        lons[1:, :],  lats[1:, :])
    
    # Calculate dx (West-East)
    _, _, dx = geod.inv(lons[:, :-1], lats[:, :-1], 
                        lons[:, 1:],  lats[:, 1:])
    
    return dx, dy
