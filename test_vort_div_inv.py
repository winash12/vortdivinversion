import os.path
import xarray as xr
import numpy as np
import boto3
import metpy.calc as mpcalc
from botocore import UNSIGNED
from botocore.config import Config
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.ticker as mticker
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys
import time
if (not os.path.isfile('gfs.t12z.pgrb2.0p25.f000')):

    client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    client.download_file('noaa-gfs-bdp-pds', 'gfs.20230809/12/atmos/gfs.t12z.pgrb2.0p25.f000',
                         'gfs.t12z.pgrb2.0p25.f000')

u850 = xr.open_dataset('gfs.t12z.pgrb2.0p25.f000', engine='cfgrib',
                       backend_kwargs={'filter_by_keys':
                                       {'typeOfLevel': 'isobaricInhPa', 'shortName': 'u',
                                        'level': 850}})
u = u850.u

v850 = xr.open_dataset('gfs.t12z.pgrb2.0p25.f000', engine='cfgrib',
                       backend_kwargs={'filter_by_keys':
                                       {'typeOfLevel': 'isobaricInhPa', 'shortName': 'v',
                                        'level': 850}})
v = v850.v


# Compute the 850 hPa relative vorticity.


vort850 = mpcalc.vorticity(u, v)

# Compute the 850 hPa divergence.

div850 = mpcalc.divergence(u, v)

umask = mpcalc.bounding_box_mask(u, 5., 13.5, 191., 202.)

vmask = mpcalc.bounding_box_mask(v, 5., 13.5, 191., 202.)


vortmask = mpcalc.bounding_box_mask(vort850, 5., 13.5, 191., 202.)


divmask = mpcalc.bounding_box_mask(div850, 5., 13.5, 191., 202.)

i_bb_indices = mpcalc.find_bounding_box_indices(vortmask, 5., 13.5, 191., 202.)


o_bb_indices = mpcalc.find_bounding_box_indices(vortmask, 0., 30., 180., 220.)


dx, dy = mpcalc.lat_lon_grid_deltas(vortmask.longitude, vortmask.latitude)

upsi, vpsi = mpcalc.rotational_wind_from_inversion(umask, vmask, vortmask, dx, dy,
                                                   o_bb_indices, i_bb_indices)

recons_vorticity = mpcalc.vorticity(upsi,vpsi)

diffArray = abs(recons_vorticity.values.flatten()-vortmask.values.flatten())


print(np.allclose(recons_vorticity,vortmask,1e-05,1e-08))

minval = np.min(diffArray[np.nonzero(diffArray)])
maxval = np.max(diffArray[np.nonzero(diffArray)])

print(minval,maxval)



uchi, vchi = mpcalc.divergent_wind_from_inversion(umask, vmask, divmask, dx, dy,
                                                  o_bb_indices, i_bb_indices)

