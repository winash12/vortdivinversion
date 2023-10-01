import numpy as np
import numpy.random as rnd
import sys

dx = rnd.randn(361, 719)



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

#vectorized code

clist = []
for i in range(x_ll_subset, x_ur_subset):

    for j in range(y_ur_subset, y_ll_subset): 

        iindex[:,:] = i
        xdiff = (iindex-xindex)*dx[y_ur:y_ll,x_ll:x_ur]
    clist.append(xdiff)

cxdiff = np.concatenate(clist,axis=0)

print(cxdiff.shape)


#scalar code

xdiff1 = np.zeros((y,x))
clist1 = []
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

                y11 += 1
            x11 += 1
    clist1.append(xdiff1)

cxdiff1 = np.concatenate(clist1,axis=0)
print(cxdiff1.shape)
print(np.array_equal(xdiff,xdiff1))
