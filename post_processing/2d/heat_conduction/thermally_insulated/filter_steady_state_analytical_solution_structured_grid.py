import numpy as np

l = 1.
h = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

T = np.cos(np.pi*x/l) *  np.cos(np.pi*y/h)

output.PointData.append(T, f'Ta')
output.PointData.append(inputs[0].PointData['T'], f'T')
