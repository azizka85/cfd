import numpy as np

alpha = 1.

l = 1.
h = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

t = inputs[0].FieldData['Time']

T = np.cos(np.pi*x/l) *  np.cos(np.pi*y/h) * np.exp(-np.pi**2*alpha*(1/l/l + 1/h/h)*t)

output.PointData.append(T, f'Ta')
output.PointData.append(inputs[0].PointData['T'], f'T')
