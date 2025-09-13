import numpy as np

ul = -1.
ur = 1.

x = inputs[0].Points[:,0]

t = inputs[0].FieldData['Time']

nx = len(x)

x0 = x[(nx+1)//2]

xl = x0 + ul*t
xr = x0 + ur*t

u = (ur - ul)*x/(xr - xl) + (ul*xr - ur*xl)/(xr - xl)

u[x <= xl] = ul
u[x >= xr] = ur

output.PointData.append(u, f'ua')
output.PointData.append(inputs[0].PointData['u'], f'u')
