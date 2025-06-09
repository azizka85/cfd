import numpy as np

Gx = 1.
H = 1.

nu = 1.

y = inputs[0].Points[:,1]
t = inputs[0].FieldData['Time']

u = 0

for n in range(1,1001):
    pn = np.pi*n
    pnh = pn / H

    p1 = 2*H*H*(1. - cos(pn)) / nu / pn**3
    p2 = sin(pnh*y)
    p3 = 1. - exp(-pnh*pnh*nu*t)

    u += p1 * p2 * p3 * Gx

output.PointData.append(u, f'u1000')
output.PointData.append(inputs[0].PointData['u'], f'u')
