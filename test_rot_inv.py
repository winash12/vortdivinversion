import xarray as xr
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units
import pyproj
import sys
def main():
    setUpData()
    setUpDivergence()
    
def setUpData():

    lat = np.linspace(90.,-90.,721,endpoint=True)*units.degree
    lon = np.linspace(0.,360. ,1440,endpoint=False)*units.degree
    u = np.empty((721,1440))
    v = np.empty((721,1440))
    longitude = xr.DataArray(lon,
                       attrs={'standard_name': 'longitude', 'units': 'degrees_east'})
    latitude = xr.DataArray(lat,
                       attrs={'standard_name': 'latitude', 'units': 'degrees_north'})
    u1 = -2. * np.sin(lat/2)*units('m/s')
    v1 = - np.cos(lon)*units('m/s')
    u1,v1 = np.meshgrid(v1,u1)
    u = xr.DataArray(u1 ,coords=(latitude,longitude),dims=('latitude','longitude'))

    v = xr.DataArray(v1 ,coords=(latitude,longitude),dims=('latitude','longitude'))

    u = u.metpy.assign_crs(grid_mapping_name='latitude_longitude',earth_radius=6371229.0)
    v = v.metpy.assign_crs(grid_mapping_name='latitude_longitude',earth_radius=6371229.0)    
    vort = mpcalc.vorticity(u,v)


    

    umask = mpcalc.bounding_box_mask(u, 5., 13.5, 191., 202.)

    vmask = mpcalc.bounding_box_mask(v, 5., 13.5, 191., 202.)

    
    vortmask = mpcalc.bounding_box_mask(vort, 5., 13.5, 191., 202.)

    
    i_bb_indices = mpcalc.find_bounding_box_indices(vortmask, 5., 13.5, 191., 202.)

    
    o_bb_indices = mpcalc.find_bounding_box_indices(vortmask, 0., 30., 180., 220.)



    dx, dy = mpcalc.lat_lon_grid_deltas(vortmask.longitude, vortmask.latitude)
    upsi, vpsi = mpcalc.rotational_wind_from_inversion(umask, vmask, vortmask, dx, dy,
                                                       o_bb_indices, i_bb_indices)
    upsi = upsi*units('m/s')
    vpsi = vpsi*units('m/s')
    recovered_vort = mpcalc.vorticity(upsi,vpsi)

    dx1 = np.array(dx)
    dy1 = np.array(dy)
    aa = np.max(dx1)
    bb = np.max(dy1)
    tolerance = (max(aa,bb)) *np.abs(vort.values).max()
    print(np.allclose(recovered_vort,vort,atol=tolerance))

def setUpDivergence():
    
    lat = np.linspace(90.,-90.,721,endpoint=True)*units.degree
    lon = np.linspace(0.,360. ,1440,endpoint=False)*units.degree
    u = np.empty((721,1440))
    v = np.empty((721,1440))
    longitude = xr.DataArray(lon,
                             attrs={'standard_name': 'longitude', 'units': 'degrees_east'})
    latitude = xr.DataArray(lat,
                            attrs={'standard_name': 'latitude', 'units': 'degrees_north'})


    u1 = - np.cos(lon)*units('m/s')
    v1 = -2. * np.sin(lat/2)*units('m/s')

    u1,v1 = np.meshgrid(u1,v1)

    u = xr.DataArray(u1 ,coords=(latitude,longitude),dims=('latitude','longitude'))
    
    v = xr.DataArray(v1 ,coords=(latitude,longitude),dims=('latitude','longitude'))

    u = u.metpy.assign_crs(grid_mapping_name='latitude_longitude',earth_radius=6371229.0)
    v = v.metpy.assign_crs(grid_mapping_name='latitude_longitude',earth_radius=6371229.0)    
    
    div = mpcalc.divergence(u,v)
    
    divmask = mpcalc.bounding_box_mask(div, 5., 13.5, 191., 202.)
    
    i_bb_indices = mpcalc.find_bounding_box_indices(divmask, 5., 13.5, 191., 202.)
    
    
    o_bb_indices = mpcalc.find_bounding_box_indices(divmask, 0., 30., 180., 220.)
    
    
    dx, dy = mpcalc.lat_lon_grid_deltas(divmask.longitude, divmask.latitude)

    umask = mpcalc.bounding_box_mask(u, 5., 13.5, 191., 202.)

    vmask = mpcalc.bounding_box_mask(v, 5., 13.5, 191., 202.)


    upsi, vpsi = mpcalc.divergent_wind_from_inversion(umask, vmask, divmask, dx, dy,
                                                      o_bb_indices, i_bb_indices)
    upsi = upsi*units('m/s')
    vpsi = vpsi*units('m/s')
    recovered_div = mpcalc.divergence(upsi,vpsi)
    
    dx1 = np.array(dx)
    dy1 = np.array(dy)
    aa = np.max(dx1)
    bb = np.max(dy1)
    tolerance = (max(aa,bb)) *np.abs(div.values).max()
    print(np.allclose(recovered_div,div,atol=tolerance))

main()

