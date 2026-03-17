import numpy as np
import xarray as xr
from dataclasses import dataclass
from abc import ABC, abstractmethod

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


def bounding_box_mask(data_array, min_lat, max_lat, min_lon, max_lon):
    """
    Slices the data to the bounding box and masks values outside.
    Pre-requisite: data_array must be normalized (South-to-North).
    """
    # Slice first to minimize memory footprint
    subset = data_array.sel(
        latitude=slice(min_lat, max_lat), 
        longitude=slice(min_lon, max_lon)
    )
    
    mask = (
        (subset.latitude <= max_lat) & (subset.latitude >= min_lat) &
        (subset.longitude <= max_lon) & (subset.longitude >= min_lon)
    )
    return subset.where(mask).fillna(0.0)


def find_bounding_box_indices(data_array, min_lat, max_lat, min_lon, max_lon):
    """
    Returns array indices of a bounding box using nearest-neighbor lookup.
    """
    def get_idx(coords, val):
        return int(np.argmin(np.abs(coords - val)))

    return BoundingBoxIndices(
        x_ll=get_idx(data_array.longitude.values, min_lon),
        x_ur=get_idx(data_array.longitude.values, max_lon),
        y_ll=get_idx(data_array.latitude.values, min_lat),
        y_ur=get_idx(data_array.latitude.values, max_lat)
    )


def get_vectorized_array_indices(bb_indices):
    """
    Generates index meshgrids for the inversion integral.
    'ij' indexing ensures Axis 0 is Latitude (y) and Axis 1 is Longitude (x).
    """
    x_range = np.arange(bb_indices.x_ll, bb_indices.x_ur + 1)
    y_range = np.arange(bb_indices.y_ll, bb_indices.y_ur + 1)
    
    # yindex (lat) becomes the leading dimension
    yindex, xindex = np.meshgrid(y_range, x_range, indexing='ij')
    
    return [xindex, yindex]
