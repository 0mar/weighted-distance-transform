#!/usr/bin/env python3
import numpy as np
from fortran.wdt import weighted_distance_transform
from scipy.misc import imread
import time
from script import WDT
import matplotlib.pyplot as plt

filename = 'images/ex2.png'
data = imread(filename, mode='RGB')
cost_field = WDT.convert_image_to_cost_field(data)
time1 = time.time()
phi = weighted_distance_transform(cost_field,*cost_field.shape,50)
phi[phi>=49]=np.inf
time2 = time.time()
print(time2-time1)
plt.imshow(phi)
plt.colorbar()
plt.show()
