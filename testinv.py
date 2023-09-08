import numpy as np

from botocore import UNSIGNED
from botocore.config import Config
import sys
import os.path
import xarray as xr
from bounding_box import bounding_box,point
from testvorinv import vortdiv_inversion

def main():
    
    u850,v850=read_file()

    outer_ll_point = point(0.,180.)
    outer_ul_point = point(30.,180.)
    outer_lr_point = point(0.,220.)
    outer_ur_point = point(30.,220.)
    list_of_points = []
    list_of_points.append(outer_ll_point)
    list_of_points.append(outer_ul_point)
    list_of_points.append(outer_lr_point)
    list_of_points.append(outer_ur_point)

    inner_ll_point = point(5.,191.)
    inner_ul_point = point(13.5,191.)
    inner_lr_point = point(5.,202.)
    inner_ur_point = point(13.5,202.)
    list_Of_points2 = []
    list_Of_points2.append(inner_ll_point)
    list_of_points2.append(inner_ul_point)
    list_of_points2.append(inner_lr_point)
    list_of_points2.append(inner_ur_point)

    inner_boundingbox = bounding_box(listOfPoints2)
    outer_boundingbox = bounding_box(listOfPoints)
    vortdivinv = vortdiv_inversion(u850,v850,outerBoundingBox,innerBoundingBox)
    vortdivinv.set_data()
    vortdivinv.calculate_vorticity()
    vortdivinv.calculate_divergence()
    vortdivinv.vorticity_mask()
    vortdivinv.divergence_mask()
    vortdivinv.calculate_distance_matrix()
    vortdivinv. initiailize_rotational_and_irrotationalwind()
    vortdivinv.define_boundingbox_indices()
    vortdivinv.calculate_rotational_wind_from_inversion()
def read_file():

    if (not os.path.isfile('gfs.t12z.pgrb2.0p50.f000')):
        
        client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        client.download_file('noaa-gfs-bdp-pds', 'gfs.20230824/12/atmos/gfs.t12z.pgrb2.0p25.f000', 'gfs.t12z.pgrb2.0p25.f000')
        
    u850 = xr.open_dataset('gfs.t12z.pgrb2.0p50.f000', engine='cfgrib',backend_kwargs={'filter_by_keys':{'typeOfLevel': 'isobaricInhPa', 'shortName': 'u', 'level': 850}})
    
    v850 = xr.open_dataset('gfs.t12z.pgrb2.0p50.f000', engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel': 'isobaricInhPa', 'shortName': 'v', 'level': 850}})
    return [u850,v850]
main()

