import numpy as np
import os.path
import xarray as xr
import boto3
import metpy.calc as mpcalc
from botocore import UNSIGNED
from botocore.config import Config
import sys


class vortdiv_inversion:

    def __init__(self,u850,v850,outerBoundingBox,innerBoundingBox):
        
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
        self.outerBoundingBox = outerBoundingBox
        self.innerBoundingBox = innerBoundingBox
        
    def setData(self):
        self.u = self.u850.u
        self.v = self.v850.v

    def calculateVorticity(self):
        self.vort850 = mpcalc.vorticity(self.u,self.v)
        
    def calculateDivergence(self):

        self.div850 = mpcalc.divergence(self.u,self.v)
        
    def vorticityMask(self):
        upperlatitude =     self.innerboundingbox.getUpperLeftPoint().getLatitude()
        lowerlatitude = self.innerboundingbox.getLowerLeftPoint().getLatitude()
        upperlongitude = self.innerboundingbox.getUpperLeftPoint().getLongitude()
        lowerlongitude = self.innerboundingbox.getLowerLeftPoint().getLongitude()
        mask = ((self.vort850.latitude <= upperlatitude) & (self.vort850.latitude >= lowerlatitude) & (self.vort850.longitude <= upperlongitude) & (self.vort850.longitude >= lowerlongitude))
        self.vortmask = self.vort850.where(mask)
        self.vortmask = self.vortmask.fillna(0.0)
        
    def divergenceMask(self):
        mask = ((self.div850.latitude <= 13.5) & (self.div850.latitude >= 5.0) & (self.div850.longitude <= 202.) & (self.div850.longitude >= 191.))
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
        x_ll = list(self.vortmask.longitude.values).index(191.0)
        x_ur = list(self.vortmask.longitude.values).index(202.0)
        y_ll = list(self.vortmask.latitude.values).index(5.0)
        y_ur = list(self.vortmask.latitude.values).index(13.5)
            
        x_ll_subset = list(self.vortmask.longitude.values).index(180.0)
        x_ur_subset = list(self.vortmask.longitude.values).index(220.0)
        y_ll_subset = list(self.vortmask.latitude.values).index(0.0)
        y_ur_subset = list(self.vortmask.latitude.values).index(30.0)

    def calculateRotationalWindFromInversion(self):

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
                            upsi[j,i] += vortmask[y1,x1].values * -1.0 * (ydiff / rsq) * dx[y1,x1].magnitude * dy[y1,x1].magnitude
                            vpsi[j,i] += vortmask[y1,x1].values * (xdiff / rsq) * dx[y1,x1].magnitude * dy[y1,x1].magnitude
                    
        upsi[:,:] = (1/(2*np.pi)) * upsi[:,:]
        vpsi[:,:] = (1/(2*np.pi)) * vpsi[:,:]


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


