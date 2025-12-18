import xarray as xr
import numpy as np
import numpy.random as rnd
import sys

dx = rnd.randn(361,719)
dy = rnd.randn(360,720)
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


xindex = np.linspace(x_ll,x_ur,num = x*y,endpoint=False,dtype=np.int32)
yindex = np.linspace(y_ll,y_ur,num=y*x,endpoint=False,dtype=np.int32)



xindex = xindex.reshape((y,x),order='F')
yindex = yindex.reshape((y,x),order='C')

iindex = np.zeros((y,x),dtype=np.int32)
jindex = np.zeros((y,x),dtype=np.int32)
#vectorized code

cxdiff = []
cydiff = []
for i in range(x_ll_subset, x_ur_subset):

    for j in range(y_ur_subset, y_ll_subset): 

        iindex[:,:] = i
        jindex[:,:] = j
        xdiff = (iindex-xindex)#*dx[y_ur:y_ll,x_ll:x_ur]
        ydiff = (jindex-yindex)#*dy[y_ur:y_ll,x_ll:x_ur]
        cxdiff.append(xdiff)
        cydiff.append(ydiff)

cxdiffarr = np.concatenate(cxdiff,axis=0)
cydiffarr = np.concatenate(cydiff,axis=0)

cxdiffarr = cxdiffarr.ravel()
print(cxdiffarr.shape)

#scalar code

xdiff1 = np.zeros((y,x))
ydiff1 = np.zeros((y,x))
cunvec1 = []
cunvec2 = []
wow = 0
for i in range(x_ll_subset, x_ur_subset):
    for j in range(y_ur_subset, y_ll_subset): 
        x11 = 0
        for x1 in range(x_ll, x_ur):
            y11 = 0
            for y1 in range(y_ur, y_ll):

                outer_point = [i,j]

                inner_point = [x1,y1]

                if inner_point != outer_point:
                    xdiff1[y11,x11] = (i-x1)#*dx[y1,x1]
                    ydiff1[y11,x11] = (j-y1)#*dy[y1,x1]
        
                y11 += 1
            x11 += 1
        cunvec1.append(xdiff1)
        cunvec2.append(ydiff1)
cxdiff1 = np.concatenate(cunvec1,axis=0)
cydiff1 = np.concatenate(cunvec2,axis=0)

cxdiff1 = cxdiff1.ravel()

print(np.array_equal(cxdiffarr,cxdiff1))

