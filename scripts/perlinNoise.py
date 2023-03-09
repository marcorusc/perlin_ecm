import matplotlib.pyplot as plt
from perlin_noise import PerlinNoise
from scipy.io import savemat
import numpy as np

'''
Script currently returns one uniform substrate matrix and one noise 
matrix. Future versions will increase the adaptability of substrates.
'''

'''
Parameters
'''
dx = 20                     #x-spacing
dy = 20                     #y-spacing 
dz = 20                     #z-spacing
xmin = -300 - dx//2         #minimum x-value
ymin = -300 - dy//2         #minimum y-value
zmin = 0 - dz//2            #minimum z-value
xmax = 300 + dx//2          #maximum x-value
ymax = 300 + dy//2          #maximum y-value
zmax = 0 + dz//2            #maximum z-value
octaves = 5                 #peak/trough coarseness parameter
seed = 1                    #random seed
substrateNum = 2            #number of substrate matrices required 
uniformVal = 0.5            #concentration for uniform substrate

#Rotate the generated array to match PhysiCell input
def rotated(array):
    listOfTuples = zip(*array[::-1])
    return [list(elem) for elem in listOfTuples]

#Generate Perlin noise matrix 
noise = PerlinNoise(octaves, seed)

#Convert actual positions to indices
xpix, ypix, zpix = int((xmax - xmin) / dx), int((ymax - ymin) / dy), int((zmax - zmin) / dz)

#Zero substrate matrix
substrate_1 = np.full((xpix, ypix, zpix), uniformVal)

#Array to hold noise values
substrate_2 = [[[(noise([k / zpix, j / ypix, i / xpix]) + 1) / 2 for k in range(zpix)] for j \
	             in range(ypix)] for i in range(xpix)]
	             
#Rotate substrates to match PhysiCell
substrate_2 = rotated(substrate_2)
substrate_1 = rotated(substrate_1)

#Display generate noise array 
plt.imshow(substrate_2, cmap='gray')
plt.colorbar()
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
				
savemat(r"../config/initialConditions.mat", {'dataArray': dataArray}, format='4')
			
			