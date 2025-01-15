import numpy as np

def create_convergence_over_time(time_name: str, value_name: str, source_name: str, data_name: str):
    input_source = FindSource(source_name)
    
    timesteps = input_source.TimestepValues

    forceTime = ForceTime(Input=input_source)

    v = None

    points = []

    for timestep in timesteps:          
        forceTime.ForcedTime = timestep

        data = servermanager.Fetch(forceTime)

        time = data.GetFieldData().GetArray(time_name).GetValue(0)
        w = np.array([data.GetPointData().GetArray(value_name).GetValue(j) for j in range(data.GetNumberOfPoints())])

        if v is None:
            v = w
        else:
            value = np.abs(v - w).max()

            points.extend([time, value, 0.])

            v = w

    polyLine = PolyLineSource()

    polyLine.Points = points

    RenameSource(data_name, polyLine)
    Delete(forceTime)
