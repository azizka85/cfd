import numpy as np

A0 = 1.
A1 = 1.

L = 1.

alpha = 1.

x = inputs[0].Points[:, 0]
t = inputs[0].FieldData['Time']

PL = np.pi / L
PS = PL**2

T = A0 + A1 * exp(-alpha*PS*t) * cos(PL * x)

output.PointData.append(T, f'T_Analytical')
output.PointData.append(inputs[0].PointData['T'], f'T')
