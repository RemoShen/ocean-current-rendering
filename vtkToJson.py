import json
import vtk

reader = vtk.vtkDataSetReader()
reader.SetFileName("./public/ocean_data.vtk")
reader.Update()
grid = reader.GetOutput()

points = grid.GetPoints()
salinity = grid.GetPointData().GetArray("salinity")
temperature = grid.GetPointData().GetArray("temperature")
velocity = grid.GetPointData().GetArray("velocity")

data = {
    "points": [],
    "salinity": [],
    "temperature": [],
    "velocity": []
}

for i in range(points.GetNumberOfPoints()):
    x, y, z = points.GetPoint(i)
    data["points"].append([x, y, z])

    data["salinity"].append(salinity.GetTuple1(i))
    data["temperature"].append(temperature.GetTuple1(i))
    velocity_tuple = velocity.GetTuple3(i)
    data["velocity"].append(
        [velocity_tuple[0], velocity_tuple[1], velocity_tuple[2]])


with open("./public/ocean_data.json", "w") as outfile:
    json.dump(data, outfile)
