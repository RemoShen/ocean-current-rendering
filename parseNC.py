import xarray
import vtk

if __name__ == "__main__":
    nc_file = r"./public/source.nc"
    ds = xarray.open_dataset(nc_file)
    print(ds)

    dims = [ds.dims['longitude'], ds.dims['latitude'], ds.dims['depth']]
    salt = ds.salt.data.tolist()[0]
    temp = ds.temp.data.tolist()[0]
    u = ds.u.data.tolist()[0]
    v = ds.v.data.tolist()[0]
    w = ds.w.data.tolist()[0]

    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(dims[0], dims[1], dims[2])

    points = vtk.vtkPoints()
    points.Allocate(dims[0] * dims[1] * dims[2])
    for depth in ds.depth.data.tolist():
        for lat in ds.latitude.data.tolist():
            for lon in ds.longitude.data.tolist():
                points.InsertNextPoint(lon, lat, depth)

    # zeta = vtk.vtkFloatArray()
    salinity = vtk.vtkFloatArray()
    temperature = vtk.vtkFloatArray()
    velocity = vtk.vtkFloatArray()

    # zeta.SetNumberOfComponents(1)
    salinity.SetNumberOfComponents(1)
    temperature.SetNumberOfComponents(1)
    velocity.SetNumberOfComponents(3)

    # zeta.SetNumberOfTuples(dims[0] * dims[1] * dims[2])
    salinity.SetNumberOfTuples(dims[0] * dims[1] * dims[2])
    temperature.SetNumberOfTuples(dims[0] * dims[1] * dims[2])
    velocity.SetNumberOfTuples(dims[0] * dims[1] * dims[2])

    # zeta.SetName("zeta")
    salinity.SetName("salinity")
    temperature.SetName("temperature")
    velocity.SetName("velocity")
    for z in range(dims[2]):
        for y in range(dims[1]):
            for x in range(dims[0]):
                index = z * dims[1] * dims[0] + y * dims[0] + x
                # print(index)
                salinity.InsertTuple(index,
                                     [salt[z][y][x] if salt[z][y][x] == salt[z][y][x] else 0.0])
                temperature.InsertTuple(index,
                                        [temp[z][y][x] if temp[z][y][x] == temp[z][y][x] else 0.0])
                velocity.InsertTuple(index,
                                     [u[z][y][x] if u[z][y][x] == u[z][y][x] else 0.0,  #
                                      v[z][y][x] if v[z][y][x] == v[z][y][x] else 0.0,  #
                                      w[z][y][x] if w[z][y][x] == w[z][y][x] else 0.0]) #
    print("Done")
    grid.SetPoints(points)
    grid.GetPointData().AddArray(salinity)
    grid.GetPointData().AddArray(temperature)
    grid.GetPointData().AddArray(velocity)

    writer = vtk.vtkDataSetWriter()
    writer.SetFileName("./public/ocean_data.vtk")
    writer.SetInputData(grid)
    writer.Write()