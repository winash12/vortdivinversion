import os.path
import xarray as xr
import numpy as np
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import sys
import time
import pandas as pd
if (not os.path.isfile('gfs.t12z.pgrb2.0p25.f000')):

    client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    client.download_file('noaa-gfs-bdp-pds', 'gfs.20230809/12/atmos/gfs.t12z.pgrb2.0p25.f000',
                         'gfs.t12z.pgrb2.0p25.f000')

u850 = xr.open_dataset('gfs.t12z.pgrb2.0p25.f000', engine='cfgrib',
                       backend_kwargs={'filter_by_keys':
                                       {'typeOfLevel': 'isobaricInhPa', 'shortName': 'u',

                                        'level': 850}})

u850 = u850.to_dataframe().unstack()

print(type(u850))

u850.to_csv('out.csv')
