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
        self.upperlatitude = self.innerboundingbox.getUpperLeftPoint().getLatitude()
        self.lowerlatitude = self.innerboundingbox.getLowerLeftPoint().getLatitude()
        self.upperlongitude = self.innerboundingbox.getUpperRightPoint().getLongitude()
        self.lowerlongitude = self.innerboundingbox.getLowerLeftPoint().getLongitude()

        self.oupperlatitude = self.outerboundingbox.getUpperLeftPoint().getLatitude()
        self.olowerlatitude = self.outerboundingbox.getLowerLeftPoint().getLatitude()

        self.oupperlongitude = self.outerboundingbox.getUpperRightPoint().getLongitude()
        self.olowerlongitude = self.outerboundingbox.getLowerLeftPoint().getLongitude()

        self.x_ll = None
        self.x_ur = None
        self.y_ll = None
        self.y_ur = None
            
        self.x_ll_subset = None
        self.x_ur_subset = None
        self.y_ll_subset = None
        self.y_ur_subset = None

        
    def setData(self):
        self.u = self.u850.u
        self.v = self.v850.v

    def calculateVorticity(self):
        self.vort850 = mpcalc.vorticity(self.u,self.v)
        
    def calculateDivergence(self):

        self.div850 = mpcalc.divergence(self.u,self.v)
        
    def vorticityMask(self):

        mask = ((self.vort850.latitude <= self.upperlatitude) & (self.vort850.latitude >= self.lowerlatitude) & (self.vort850.longitude <= self.upperlongitude) & (self.vort850.longitude >= self.lowerlongitude))
        self.vortmask = self.vort850.where(mask)
        self.vortmask = self.vortmask.fillna(0.0)
        
    def divergenceMask(self):
        mask = ((self.div850.latitude <= self.upperlatitude) & (self.div850.latitude >= self.lowerlatitude) & (self.div850.longitude <= self.upperlongitude) & (self.div850.longitude >= self.lowerlongitude))
        self.divmask = self.div850.where(mask)
        self.divmask = self.div850.fillna(0.0)

        
    def calculateDistanceMatrix(self):
        self.dx, self.dy = mpcalc.lat_lon_grid_deltas(self.vortmask.longitude, self.vortmask.latitude)
        self.dx = np.abs(self.dx)
        self.dy = np.abs(self.dy)

    def initiailizeRotationalAndIrrotationalWind(self):
        self.upsi = xr.zeros_like(self.vortmask)
        self.vpsi = xr.zeros_like(self.vortmask)
        self.uchi = xr.zeros_like(self.divmask)
        self.vchi = xr.zeros_like(self.divmask)
        
    def defineBoundingBoxIndices(self):
        self.x_ll = list(self.vortmask.longitude.values).index(self.lowerlongitude)
        self.x_ur = list(self.vortmask.longitude.values).index(self.upperlongitude)
        self.y_ll = list(self.vortmask.latitude.values).index(self.lowerlatitude)
        self.y_ur = list(self.vortmask.latitude.values).index(self.upperlatitude)

        
        self.x_ll_subset = list(self.vortmask.longitude.values).index(self.olowerlongitude)
        self.x_ur_subset = list(self.vortmask.longitude.values).index(self.oupperlongitude)



        self.y_ll_subset = list(self.vortmask.latitude.values).index(self.olowerlatitude)
        self.y_ur_subset = list(self.vortmask.latitude.values).index(self.oupperlatitude)

    def calculateRotationalWindFromInversion(self):

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


    def calculateDivergentWindFromInversion(self):
        
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


