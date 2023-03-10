import matplotlib.pyplot as plt
from perlin_noise import PerlinNoise
from scipy.io import savemat
import numpy as np
import subprocess
import sys

#Install non-installed package
def install(package):
	subprocess.check_call([sys.executable, "-m", "pip", "install", package])

#Install perlin_noise package
if not ("perlin_noise" in sys.modules): 
    print('installing perlin_noise')
    install("perlin_noise")
else: 
    print('perlin_noise already installed')

'''
Script currently returns one uniform substrate matrix and one noise 
matrix. Future versions will increase the adaptability of substrates.
'''

import xml.etree.ElementTree as ET
tree = ET.parse("../config/PhysiCell_settings.xml")
root = tree.getroot()

'''
Parameters
'''
xmin = int(root[0][0].text) 
xmax = int(root[0][1].text) 
ymin = int(root[0][2].text) 
ymax = int(root[0][3].text) 
zmin = int(root[0][4].text) 
zmax = int(root[0][5].text) 
dx = int(root[0][6].text) 
dy = int(root[0][7].text) 
dz = int(root[0][8].text) 
sub_str = str( (root[5][-2]).attrib)
for c in sub_str:
    if c.isdigit():
        substrate_max_index = c
substrateNum = int(substrate_max_index) + 1

#collagen_string = "ECM_collagen"
fibrin_string = "ECM_fibrin" 
oxygen_string = "oxygen" 
for child in root[5]:
    c = str(child.attrib)
    if oxygen_string in c:
        oxy_val = float(child[1].text)
    elif fibrin_string in c: 
        fibrin_val = float(child[1].text)
        
    
octaves = int(root[-1][-2].text)
seed = int(root[-1][-1].text)


#Rotate the generated array to match PhysiCell input
def rotated(array):
    listOfTuples = zip(*array[::-1])
    return [list(elem) for elem in listOfTuples]

#Generate Perlin noise matrix 
noise = PerlinNoise(octaves, seed)

#Convert actual positions to indices
xpix, ypix, zpix = int((xmax - xmin) / dx), int((ymax - ymin) / dy), int((zmax - zmin) / dz)

#Zero substrate matrix
substrate_1 = np.full((xpix, ypix, zpix), fibrin_val)

#Array to hold noise values
#substrate_2 = [[[(noise([k / zpix, j / ypix, i / xpix]) + 1) / 2 for k in range(zpix)] for j \
	            # in range(ypix)] for i in range(xpix)]
substrate_2 = [[[(noise([i / xpix, j / ypix, k / zpix]) + 1) / 2 for k in range(zpix)] for j in range(ypix)] for i in range(xpix)]
	             
#Array to hold noise values
substrate_3 = np.full((xpix, ypix, zpix), oxy_val)
	             
#Rotate substrates to match PhysiCell
substrate_2 = rotated(substrate_2)
substrate_1 = rotated(substrate_1)
substrate_3 = rotated(substrate_3)

#Display generate noise array 
plt.imshow(substrate_2[:][:][0], cmap='gray')
plt.title("Perlin noise at $z=0$")
plt.colorbar()
plt.show()

#Generate data matrix for output
dataArray = np.zeros((3+1+substrateNum, xpix * ypix * zpix))

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
			dataArray[6,n] = substrate_3[i][j][k]
				
savemat(r"../config/initialConditions.mat", {'dataArray': dataArray}, format='4')
			
			
