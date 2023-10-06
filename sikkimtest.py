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
if (not os.path.isfile('gfs.t18z.pgrb2.0p25.f000')):

    client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    client.download_file('noaa-gfs-bdp-pds', 'gfs.20231004/18/atmos/gfs.t18z.pgrb2.0p25.f000', 'gfs.t18z.pgrb2.0p25.f000')
u850 = xr.open_dataset('gfs.t18z.pgrb2.0p25.f000', engine='cfgrib',backend_kwargs={'filter_by_keys':{'typeOfLevel': 'isobaricInhPa', 'shortName': 'u', 'level': 200}})
u = u850.u
print(u.shape)
v850 = xr.open_dataset('gfs.t18z.pgrb2.0p25.f000', engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel': 'isobaricInhPa', 'shortName': 'v', 'level': 200}})
v = v850.v

print(v.shape)

# Compute the 850 hPa relative vorticity.
vort850 = mpcalc.vorticity(u, v)


# Compute the 850 hPa divergence.
div850 = mpcalc.divergence(u, v)



mask = ((vort850.latitude <= 35) & (vort850.latitude >= 14.0) & (vort850.longitude <= 100.) & (vort850.longitude >= 70.))
vortmask = vort850.where(mask)

vortmask = vortmask.fillna(0.0)
divmask = div850.where(mask)
divmas = div850.fillna(0.0)
dx, dy = mpcalc.lat_lon_grid_deltas(vortmask.longitude, vortmask.latitude)

print(dy.shape)


upsi = xr.zeros_like(vortmask)

vpsi = xr.zeros_like(vortmask)



x_ll = list(vortmask.longitude.values).index(80.0)
x_ur = list(vortmask.longitude.values).index(93.0)
y_ll = list(vortmask.latitude.values).index(18.0)
y_ur = list(vortmask.latitude.values).index(36.0)

print(x_ll,x_ur,y_ll,y_ur)



x_ll_subset = list(vortmask.longitude.values).index(70.0)
x_ur_subset = list(vortmask.longitude.values).index(100.0)
y_ll_subset = list(vortmask.latitude.values).index(10.0)
y_ur_subset = list(vortmask.latitude.values).index(36.0)

print(x_ll_subset,x_ur_subset,y_ll_subset,y_ur_subset)

i = np.abs(x_ll_subset-x_ur_subset)
j = np.abs(y_ll_subset-y_ur_subset)



x = np.abs(x_ll-x_ur)
y = np.abs(y_ll-y_ur)


xindex = np.linspace(x_ll,x_ur,num=x*y,endpoint=False,dtype=np.int32)
yindex = np.linspace(y_ur,y_ll,num=y*x,endpoint=False,dtype=np.int32)


xindex = xindex.reshape((y,x),order='F')
yindex = yindex.reshape((y,x),order='C')




iindex = np.zeros((y,x),dtype=np.int32)
jindex = np.zeros((y,x),dtype=np.int32)


dx1 = dx.magnitude

dy1 = dy.magnitude

vortmask1 = vortmask.values

starttime = time.time()

for i in range(x_ll_subset, x_ur_subset):

    
    for j in range(y_ur_subset, y_ll_subset): 

        iindex[:,:] = i
        jindex[:,:] = j
        xdiff = (iindex-xindex)*dx1[y_ur:y_ll,x_ll:x_ur]
        ydiff = (jindex-yindex)*dy1[y_ur:y_ll,x_ll:x_ur]
        rsq = xdiff * xdiff + ydiff * ydiff
        upsi[j,i] = np.where(rsq > 0., vortmask1[y_ur:y_ll,x_ll:x_ur]*-1.0*(ydiff/rsq)*dx1[y_ur:y_ll,x_ll:x_ur]*dy1[y_ur:y_ll,x_ll:x_ur], 0.0).sum()
        vpsi[j,i] = np.where(rsq > 0., vortmask1[y_ur:y_ll,x_ll:x_ur]*-1.0*(xdiff/rsq)*dx1[y_ur:y_ll,x_ll:x_ur]*dy1[y_ur:y_ll,x_ll:x_ur], 0.0).sum()




upsi[:,:] = (1/(2*np.pi)) * upsi[:,:]
vpsi[:,:] = (1/(2*np.pi)) * vpsi[:,:]

endtime = time.time()

print(endtime-starttime)

print("done vectorizing")

# Create a set of axes for the figure and set
# its map projection to that of the input data.

ax = plt.axes(projection=crs.PlateCarree())

# Add country borders and coastlines.
countries = NaturalEarthFeature(category="cultural", scale="50m",
                                      facecolor="none",
                                      name="admin_0_countries")
ax.add_feature(countries, linewidth=.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)

# Compute the magnitude of the non-divergent component of the 850 hPa wind.
nd_spd = np.sqrt(upsi**2 + vpsi**2)

# Plot this using xarray's plot functionality.
plot = nd_spd.plot(levels=np.arange(0., 13., 1.), cmap=get_cmap('YlGnBu'), transform=crs.PlateCarree(), cbar_kwargs={'label':'non-divergent wind ($m s^{-1}$)', 'shrink': 0.98})

# Set the map's extent to match that over which we computed the non-divergent wind.
ax.set_extent([70,100.,6.,36.],crs=crs.PlateCarree())

# Add latitude/longitude gridlines.
gridlines = ax.gridlines(color="grey", linestyle="dotted", draw_labels=True)
gridlines.xlabels_top = False
gridlines.ylabels_right = False
gridlines.xlocator = mticker.FixedLocator(np.arange(70.,100.,5.))
gridlines.ylocator = mticker.FixedLocator(np.arange(6.,36.,5.))
gridlines.xlabel_style = {'size':12, 'color':'black'}
gridlines.ylabel_style = {'size':12, 'color':'black'}
gridlines.xformatter = LONGITUDE_FORMATTER
gridlines.yformatter = LATITUDE_FORMATTER

# Add a plot title, then show the image.
plt.title("GFS 0-h 200 hPa non-divergent wind magnitude ($m s^{-1}$) due to tropical low at 1800 UTC 4 October 2023")
plt.savefig('vectorized_version')
plt.show()

sys.exit()

        

     
