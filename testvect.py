import xarray as xr
import numpy as np
import numpy.random as rnd
import sys
import metpy.calc as mpcalc

dx = rnd.randn(361, 719)
dy = rnd.randn(361, 719)
u850 = xr.open_dataset('gfs.t12z.pgrb2.0p50.f000', engine='cfgrib',backend_kwargs={'filter_by_keys':{'typeOfLevel': 'isobaricInhPa', 'shortName': 'u', 'level': 850}})
u = u850.u
v850 = xr.open_dataset('gfs.t12z.pgrb2.0p50.f000', engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel': 'isobaricInhPa', 'shortName': 'v', 'level': 850}})
v = v850.v
vort850 = mpcalc.vorticity(u, v)

mask = ((vort850.latitude <= 13.5) & (vort850.latitude >= 5.0) & (vort850.longitude <= 202.) & (vort850.longitude >= 191.))
vortmask = vort850.where(mask)
vortmask = vortmask.fillna(0.0)

#dx, dy = mpcalc.lat_lon_grid_deltas(vortmask.longitude, vortmask.latitude)

#dx = np.abs(dx)
#dy = np.abs(dy)

#dx = dx.magnitude
#dy = dy.magnitude

x_ll_subset = 360
x_ur_subset = 440
y_ll_subset = 180
y_ur_subset = 120

#
x_ll = 382
x_ur = 404
y_ll = 170
y_ur = 153

x = np.abs(x_ll-x_ur)
y = np.abs(y_ll-y_ur)


xstart = np.linspace(x_ll,x_ur,num = x,endpoint=False,dtype=np.int32)
ystart = np.linspace(y_ll,y_ur,num=y,endpoint=False,dtype=np.int32)

xindex,yindex = np.meshgrid(xstart,ystart)

iindex = np.zeros((y,x),dtype=np.int32)
jindex = np.zeros((y,x),dtype=np.int32)
#vectorized code

clist = []
clist1 = []
for i in range(x_ll_subset, x_ur_subset):

    for j in range(y_ur_subset, y_ll_subset): 

        iindex[:,:] = i
        jindex[:,:] = j
        xdiff = (iindex-xindex)*dx[y_ur:y_ll,x_ll:x_ur]
        ydiff = (jindex-yindex)*dy[y_ur:y_ll,x_ll:x_ur]
    clist.append(xdiff)
    clist1.append(ydiff)
cxdiff = np.concatenate(clist,axis=0)
cydiff = np.concatenate(clist1,axis=0)

print(cxdiff.shape)
print(cydiff.shape)

#scalar code

xdiff1 = np.zeros((y,x))
ydiff1 = np.zeros((y,x))
clistsc = []
clistsc1 = []
for i in range(x_ll_subset, x_ur_subset):
    for j in range(y_ur_subset, y_ll_subset): 
        x11 = 0
        for x1 in range(x_ll, x_ur):
            y11 = 0
            for y1 in range(y_ur, y_ll):

                outer_point = [i,j]

                inner_point = [x1,y1]

                if inner_point != outer_point:
                    xdiff1[y11,x11] = (i-x1)*dx[y1,x1]
                    ydiff1[y11,x11] = (j-y1)*dy[y1,x1]
                y11 += 1
            x11 += 1
    clistsc.append(xdiff1)
    clistsc1.append(ydiff1)
cxdiff1 = np.concatenate(clistsc,axis=0)
cydiff1 = np.concatenate(clistsc1,axis=0)
print(cxdiff1.shape)
print(np.array_equal(cxdiff,cxdiff1))
#print(np.array_equal(ydiff,ydiff1))
