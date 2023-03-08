import matplotlib.pyplot as plt
from perlin_noise import PerlinNoise
from scipy.io import savemat
import numpy as np

dx = 20
dy = 20
dz = 20
xmin = -300 - dx//2
ymin = -300 - dy//2
zmin = 0 - dz//2
xmax = 300 + dx//2
ymax = 300 + dy//2
zmax = 0 + dz//2
octaves = 5
seed = 1
substrateNum = 1

#Generate Perlin noise matrix and output as MATLAB file
noise = PerlinNoise(octaves, seed)

#Convert actual dimensions to indices
xpix, ypix, zpix = int((xmax - xmin) / dx), int((ymax - ymin) / dy), int((zmax - zmin) / dz)

#Zero substrate matrix
substrate_1 = np.zeros((xpix, ypix, zpix))

#Array to hold noise values
substrate_2 = [[[(noise([k / zpix, j / ypix, i / xpix]) + 1) / 2 for k in range(zpix)] for j \
	             in range(ypix)] for i in range(xpix)]

#Display generate noise array 
plt.imshow(substrate_2, cmap='gray')
plt.colorbar()
print(np.min(substrate_2))
plt.show()

#Generate data matrix for output
dataArray = np.zeros((3+1+2, xpix * ypix * zpix))

#Insert substrate, noise information into data matrix
for k in range(zpix): 
	for j in range(ypix): 
		for i in range(xpix):
			n = k * (zpix*ypix) + j*ypix + i
			dataArray[0,n] = int((i * dx + xmin) + dx // 2)
			dataArray[1,n] = int((j * dy + ymin) + dy // 2)
			dataArray[2,n] = int((k * dz + zmin) + dz // 2)
			dataArray[3,n] = dx * dy * dz 
			#Substrate information 
			dataArray[5,n] = substrate_1[i][j][k]
			dataArray[4,n] = substrate_2[i][j][k]
				
savemat('initialConditions.mat', {'dataArray': dataArray}, format='4')
			
			