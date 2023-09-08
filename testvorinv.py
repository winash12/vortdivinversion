import numpy as np
import os.path
import xarray as xr
import boto3
import metpy.calc as mpcalc
from botocore import UNSIGNED
from botocore.config import Config
import sys


class vortdiv_inversion:

    def __init__(self,u850,v850,outerboundingbox,innerboundingbox):
        
        self.u850 = u850
        self.v850 = v850
        self.u = None
        self.v = None
        self.vort850 = None
        self.div850 = None
        self.vortmask = None
        self.divmask = None
        self.dx = None
        self.dy = None
        self.upsi = None
        self.vpsi = None
        self.uchi = None
        self.vchi = None
        self.outerboundingbox = outerboundingbox
        self.innerboundingbox = innerboundingbox
        self.upper_latitude = self.innerboundingbox.get_upper_left_point().get_latitude()
        self.lower_latitude = self.innerboundingbox.get_lower_left_point().get_latitude()
        self.upper_longitude = self.innerboundingbox.get_upper_right_point().get_longitude()
        self.lower_longitude = self.innerboundingbox.get_lower_left_point().get_longitude()

        self.oupper_latitude = self.outerboundingbox.get_upper_left_point().get_latitude()
        self.olower_latitude = self.outerboundingbox.get_lower_left_point().get_latitude()

        self.oupper_longitude = self.outerboundingbox.get_upper_right_point().get_longitude()
        self.olower_longitude = self.outerboundingbox.get_lower_left_point().get_longitude()

        self.x_ll = None
        self.x_ur = None
        self.y_ll = None
        self.y_ur = None
            
        self.x_ll_subset = None
        self.x_ur_subset = None
        self.y_ll_subset = None
        self.y_ur_subset = None

        
    def set_data(self):
        self.u = self.u850.u
        self.v = self.v850.v

    def calculate_vorticity(self):
        self.vort850 = mpcalc.vorticity(self.u,self.v)
        
    def calculate_divergence(self):

        self.div850 = mpcalc.divergence(self.u,self.v)
        
    def vorticity_mask(self):

        mask = ((self.vort850.latitude <= self.upper_latitude) & (self.vort850.latitude >= self.lower_latitude) & (self.vort850.longitude <= self.upper_longitude) & (self.vort850.longitude >= self.lower_longitude))
        self.vortmask = self.vort850.where(mask)
        self.vortmask = self.vortmask.fillna(0.0)
        
    def divergence_mask(self):
        mask = ((self.div850.latitude <= self.upper_latitude) & (self.div850.latitude >= self.lower_latitude) & (self.div850.longitude <= self.upper_longitude) & (self.div850.longitude >= self.lower_longitude))
        self.divmask = self.div850.where(mask)
        self.divmask = self.div850.fillna(0.0)

        
    def calculate_distance_matrix(self):
        self.dx, self.dy = mpcalc.lat_lon_grid_deltas(self.vortmask.longitude, self.vortmask.latitude)
        self.dx = np.abs(self.dx)
        self.dy = np.abs(self.dy)

    def initiailize_rotational_and_irrotational_wind(self):
        self.upsi = xr.zeros_like(self.vortmask)
        self.vpsi = xr.zeros_like(self.vortmask)
        self.uchi = xr.zeros_like(self.divmask)
        self.vchi = xr.zeros_like(self.divmask)
        
    def define_bounding_box_indices(self):
        self.x_ll = list(self.vortmask.longitude.values).index(self.lower_longitude)
        self.x_ur = list(self.vortmask.longitude.values).index(self.upper_longitude)
        self.y_ll = list(self.vortmask.latitude.values).index(self.lower_latitude)
        self.y_ur = list(self.vortmask.latitude.values).index(self.upper_latitude)

        
        self.x_ll_subset = list(self.vortmask.longitude.values).index(self.olower_longitude)
        self.x_ur_subset = list(self.vortmask.longitude.values).index(self.oupper_longitude)



        self.y_ll_subset = list(self.vortmask.latitude.values).index(self.olower_latitude)
        self.y_ur_subset = list(self.vortmask.latitude.values).index(self.oupper_latitude)

    def calculate_rotational_wind_from_inversion(self):

        for i in range(self.x_ll_subset, self.x_ur_subset):
            for j in range(self.y_ur_subset, self.y_ll_subset): 
                # Computing the contribution to each other grid point (masked area).
                # x1 and y1 refer to x' and y' in the Oertel and Schemm (2021) equations.
                for x1 in range(self.x_ll, self.x_ur):
                    for y1 in range(self.y_ur, self.y_ll):
                        outer_point = [i,j]
                        inner_point = [x1,y1]
                        if inner_point != outer_point:

                            # Compute x-x', y-y', and r^2...
                            xdiff = (i-x1)*self.dx[y1,x1].magnitude
                            ydiff = (j-y1)*self.dy[y1,x1].magnitude
                            rsq = (xdiff*xdiff) + (ydiff*ydiff)
                            self.upsi[j,i] += self.vortmask[y1,x1].values * -1.0 * (ydiff / rsq) * self.dx[y1,x1].magnitude * self.dy[y1,x1].magnitude
                            self.vpsi[j,i] += self.vortmask[y1,x1].values * (xdiff / rsq) * self.dx[y1,x1].magnitude * self.dy[y1,x1].magnitude


        self.upsi[:,:] = (1/(2*np.pi)) * self.upsi[:,:]
        self.vpsi[:,:] = (1/(2*np.pi)) * self.vpsi[:,:]


    def calculate_divergent_wind_from_inversion(self):
        
        for i in range(x_ll_subset, x_ur_subset):
            for j in range(y_ur_subset, y_ll_subset): 
                
                # Computing the contribution to each other grid point (masked area).
                # x1 and y1 refer to x' and y' in the Oertel and Schemm (2021) equations.
                for x1 in range(x_ll, x_ur):
                    for y1 in range(y_ur, y_ll):
                        outer_point = [i,j]
                        inner_point = [x1,y1]
                        if inner_point != outer_point:
                            
                            # Compute x-x', y-y', and r^2...
                            xdiff = (i-x1)*dx[y1,x1].magnitude
                            ydiff = (j-y1)*dy[y1,x1].magnitude
                            rsq = (xdiff*xdiff) + (ydiff*ydiff)
                            
                            # Compute the irrotational flow contribution.
                            uchi[j,i] += divmask[y1,x1].values * (xdiff / rsq) * dx[y1,x1].magnitude * dy[y1,x1].magnitude
                            vchi[j,i] += divmask[y1,x1].values * (ydiff / rsq) * dx[y1,x1].magnitude * dy[y1,x1].magnitude
                    

        uchi[:,:] = (1/(2*np.pi)) * uchi[:,:]
        vchi[:,:] = (1/(2*np.pi)) * vchi[:,:]


