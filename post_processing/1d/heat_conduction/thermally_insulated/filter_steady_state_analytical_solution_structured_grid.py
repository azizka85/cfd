import numpy as np

A0 = 1.

T = A0

output.PointData.append(T, f'T_Analytical')
output.PointData.append(inputs[0].PointData['T'], f'T')
