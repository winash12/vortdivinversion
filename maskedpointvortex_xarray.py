import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import sys

def main():
    setUpData()

def setUpData():

    
    # Define the vortex strength and position
    Gamma = 1.0
    lon0, lat0 = 0.0,0.0

    lat = np.linspace(90.,-90.,181,endpoint=True)*units.degree
    lon = np.linspace(-180.,180. ,360,endpoint=False)*units.degree

    u = np.empty((181,360))
    v = np.empty((181,360))

    longitude = xr.DataArray(lon,
                       attrs={'standard_name': 'longitude', 'units': 'degrees_east'})
    latitude = xr.DataArray(lat,
                       attrs={'standard_name': 'latitude', 'units': 'degrees_north'})
    
    lat,lon = np.meshgrid(lon,lat)
    
    # Calculate the velocity field
    x = np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon))
    y = np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon))
    x0 = np.cos(np.deg2rad(lat0)) * np.cos(np.deg2rad(lon0))
    y0 = np.cos(np.deg2rad(lat0)) * np.sin(np.deg2rad(lon0))
    u = -Gamma/(2*np.pi)*(y-y0)/((x-x0)**2 + (y-y0)**2)*units('m/s')
    v = Gamma/(2*np.pi)*(x-x0)/((x-x0)**2 + (x-x0)**2)*units('m/s')


    u1 = xr.DataArray(u ,coords=(latitude,longitude),dims=('latitude','longitude'))

    v1 = xr.DataArray(v ,coords=(latitude,longitude),dims=('latitude','longitude'))

    vort = mpcalc.vorticity(u1,v1)
    
    umask = mpcalc.bounding_box_mask(u1, -1., 1, -1., 1.)

    vmask = mpcalc.bounding_box_mask(v1, -1., 1, -1., 1.)

    vortmask = mpcalc.bounding_box_mask(vort, -1., 1., -1., 1.)

    i_bb_indices = mpcalc.find_bounding_box_indices(vortmask, -1., 1., -1., 1.)
    
    o_bb_indices = mpcalc.find_bounding_box_indices(vortmask, -2., 2., -2., 2.)


    dx, dy = mpcalc.lat_lon_grid_deltas(vortmask.longitude, vortmask.latitude)

    upsi, vpsi = mpcalc.rotational_wind_from_inversion(umask, vmask, vortmask, dx, dy,
                                                       o_bb_indices, i_bb_indices)

    upsi = upsi*units('m/s')
    vpsi = vpsi*units('m/s')
    recovered_vort = mpcalc.vorticity(upsi,vpsi)

    upsi = upsi.values
    umask = umask.values
    print(np.allclose(umask,upsi,atol=0.6))
    sys.exit()

main()
