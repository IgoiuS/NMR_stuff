#Shileding Tensor Diagonalizator

import numpy as np
from scipy.linalg import block_diag

#path of the output file from Turbomole
#path = 'C:/Users/igordas/OneDrive/Bachelors thesis Wiebke/DFT/Ru_triphos/shifts/shifts.out'
path = 'C:/Users/igordas/OneDrive/Bachelors thesis Wiebke/DFT/wo_Ru/shifts.out'

fhandle = open(path)

data = []

for line in fhandle:
    line = line.strip()
    if len(line) == 0:
        continue
    else:
        data.append(line)
        
tensor = []
name = []
for i in range(len(data)):
    if data[i][0:7] == 'ATOM  p':
        name.append(data[i].strip())
        tensor.append(data[i+22:i+25])

tensor_num = []

for i in range(3):
    for j in range(3):
        for k in range(3):
            tensor_num.append(float(tensor[i][j].split()[k]))

shield_block = block_diag(np.asarray(tensor_num[0:9]).reshape(3,3), 
                          np.asarray(tensor_num[9:18]).reshape(3,3), 
                          np.asarray(tensor_num[18:27]).reshape(3,3))


diagonal = np.linalg.eigvals(shield_block)

iso = []
aniso = []
indices = [3, 6]
prev = 0
for l in indices:
    s = np.sum(diagonal[prev:l])
    m = np.max(diagonal[prev:l])
    iso.append(s/3)
    aniso.append(m-(s-m)/2)
    prev = l
iso.append(np.sum(diagonal[prev:])/3)
aniso.append(np.max(diagonal[prev:])-(np.sum(diagonal[prev:])-np.max(diagonal[prev:]))/2)

print('\n Isotropic Shielding values are: ', iso, '\n Anisotropic Shielding values are: ', aniso)