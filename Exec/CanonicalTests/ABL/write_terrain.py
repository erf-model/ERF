import numpy as np

target=open("erf_terrain_def","w")

xc=1024
yc=1024
x=np.arange(xc-256,xc+256,16)
y=np.arange(yc-256,yc+256,16)

X,Y=np.meshgrid(x,y)

for i in range(0,X.shape[0]):
    for j in range(0,X.shape[1]):
        target.write("%g %g %g\n"%(X[i,j],Y[i,j],108))
target.close()
