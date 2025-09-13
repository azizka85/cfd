import numpy as np

UL = 1.
UR = 0.
dy = 0.01

output.FieldData.append(UR, 'UR')
output.FieldData.append(UL, 'UL')
output.FieldData.append(dy, 'dy')

y = np.arange(-2, 2 + dy, dy) 

output.RowData.append(y, 'y')

def get_solution(nu):
	return (UL + UR)/2 - (UL - UR)*np.tanh(y*(UL - UR)/4/nu)/2

viscosities = [0.25, 0.1, 0.02]

for viscosity in viscosities:
	w = get_solution(viscosity)
	
	output.RowData.append(w, f'w_{viscosity}')
