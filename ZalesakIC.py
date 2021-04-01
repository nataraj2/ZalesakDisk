import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
from numpy import *

xdata = loadtxt('xvec.txt')
xc = xdata[:]
nx = size(xdata)
ydata = loadtxt('yvec.txt')
yc = ydata[:]
ny = size(ydata)
data = loadtxt('zalesakdiskIC.txt')
n = size(data[:,0])

VOF = zeros([nx,ny])

for i in arange(0,n,1):
	VOF[data[i,0],data[i,1]] = data[i,2]


Xc, Yc = np.meshgrid(xc, yc)

plt.figure(1)
for i in arange(0,nx,1):
        plt.plot(Xc[:,i],Yc[:,0],'k')
for j in arange(0,ny,1):
        plt.plot(Xc[0,:],Yc[j,:],'k')

plt.contourf(Xc, Yc, transpose(VOF), 20, colormap='hot');
plt.tight_layout()
plt.axis('equal')
plt.show()
