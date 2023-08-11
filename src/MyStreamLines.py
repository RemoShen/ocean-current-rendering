import vtk
import xarray
import math
import random
import csv
import numpy as np
from enum import Enum
from copy import deepcopy, copy


class Cell:
    def __init__(self, cell_points, cell_data):
        self.CellPoints = cell_points
        self.CellData = cell_data
        self.Bounds = (self.CellPoints[0], self.CellPoints[-1])

    def GetLength2(self):
        l = 0.0

        for i in range(3):
            diff = self.Bounds[1][i] - self.Bounds[0][i]
            l += diff * diff

        return l

    def trilinear_interpolate(self, pos):
        xd = (pos[0] - self.Bounds[0][0]) / (self.Bounds[1][0] - self.Bounds[0][0])
        yd = (pos[1] - self.Bounds[0][1]) / (self.Bounds[1][1] - self.Bounds[0][1])
        zd = (pos[2] - self.Bounds[0][2]) / (self.Bounds[1][2] - self.Bounds[0][2])

        res = np.zeros(self.CellData[0].shape)
        res += self.CellData[0] * (1 - xd) * (1 - yd) * (1 - zd)
        res += self.CellData[1] * xd * (1 - yd) * (1 - zd)
        res += self.CellData[2] * (1 - xd) * yd * (1 - zd)
        res += self.CellData[3] * xd * yd * (1 - zd)
        res += self.CellData[4] * (1 - xd) * (1 - yd) * zd
        res += self.CellData[5] * xd * (1 - yd) * zd
        res += self.CellData[6] * (1 - xd) * yd * zd
        res += self.CellData[7] * xd * yd * zd
        return res


class StructureGrid:
    dims = [0.0, 0.0, 0.0]
    points = None
    pointData = {}
    bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    def xyz2index(self, x, y, z):
        return z * self.dims[1] * self.dims[0] + y * self.dims[0] + x

    def __init__(self, dataset: xarray.Dataset):
        self.dims = [dataset.dims['longitude'], dataset.dims['latitude'], dataset.dims['depth']]
        # salt = dataset.salt.data.tolist()[0]
        # temp = dataset.temp.data.tolist()[0]
        u = dataset.u.data.tolist()[0]
        v = dataset.v.data.tolist()[0]
        w = dataset.w.data.tolist()[0]

        self.bounds = [np.inf, np.inf, np.inf, -np.inf, -np.inf, -np.inf]

        self.points = np.empty(shape=[self.dims[2], self.dims[1], self.dims[0], 3])

        for z, depth in enumerate(dataset.depth.data.tolist()):
            for y, lat in enumerate(dataset.latitude.data.tolist()):
                for x, lon in enumerate(dataset.longitude.data.tolist()):
                    self.points[z, y, x] = [lon, lat, depth]
                    if depth < self.bounds[2]:
                        self.bounds[2] = depth
                    if depth > self.bounds[5]:
                        self.bounds[5] = depth
                    if lat < self.bounds[1]:
                        self.bounds[1] = lat
                    if lat > self.bounds[4]:
                        self.bounds[4] = lat
                    if lon < self.bounds[0]:
                        self.bounds[0] = lon
                    if lon > self.bounds[3]:
                        self.bounds[3] = lon

        velocity = []
        for z in range(self.dims[2]):
            for y in range(self.dims[1]):
                for x in range(self.dims[0]):
                    index = self.xyz2index(x, y, z)
                    velocity.append(np.array([u[z][y][x] if u[z][y][x] == u[z][y][x] else 0.0,
                                              v[z][y][x] if v[z][y][x] == v[z][y][x] else 0.0,
                                              w[z][y][x] if w[z][y][x] == w[z][y][x] else 0.0]))

        self.pointData["velocity"] = velocity

    def inBounds(self, pos):
        if self.bounds[0] <= pos[0] <= self.bounds[3] and \
                self.bounds[1] <= pos[1] <= self.bounds[4] and \
                self.bounds[2] <= pos[2] <= self.bounds[5]:
            return True
        else:
            return False

    def getCell(self, pos, arrayName="velocity"):
        assert len(pos) >= 3
        if not self.inBounds(pos):
            return None
        xIndex = np.searchsorted(self.points[0, 0, :, 0], pos[0])
        yIndex = np.searchsorted(self.points[0, :, 0, 1], pos[1])
        zIndex = np.searchsorted(self.points[:, 0, 0, 2], pos[2])
        assert 0 < xIndex < self.dims[0]
        assert 0 < yIndex < self.dims[1]
        assert 0 < zIndex < self.dims[2]
        xIndex = xIndex - 1
        yIndex = yIndex - 1
        zIndex = zIndex - 1

        cell_points = []
        cell_data = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    cell_points.append(self.points[zIndex + i, yIndex + i, xIndex + i])
                    cell_data.append(
                        np.array(self.pointData[arrayName][self.xyz2index(xIndex, yIndex, zIndex)]))
        return Cell(cell_points, cell_data)

    def sample(self, pos: np.array, vec: np.array, cell=None, arrayName="velocity", useNormalize=False):
        if not self.inBounds(pos):
            return 0
        if cell is None:
            cell = self.getCell(pos, arrayName)
        inter = cell.trilinear_interpolate(pos)
        if useNormalize:
            norm_inter = np.linalg.norm(inter)
            if norm_inter != 0:
                for i in range(len(inter)):
                    vec[i] = inter[i] / norm_inter
        else:
            for i in range(len(inter)):
                vec[i] = inter[i]
        return 1


def GenerateSeeds(field: StructureGrid, num=1000):
    seeds = []
    while len(seeds) < num:
        seed = np.zeros(3)
        seed[0] = random.uniform(field.bounds[0], field.bounds[3])
        seed[1] = random.uniform(field.bounds[1], field.bounds[4])
        seed[2] = random.uniform(field.bounds[2], field.bounds[5])
        cell = field.getCell(seed)
        if np.linalg.norm(cell.trilinear_interpolate(seed)) > 0:
            # if field.inBounds(seed):
            seeds.append(seed)

    return seeds


class IntegratorRK45:
    A = np.zeros(5)
    B = np.zeros([5, 5])
    C = np.zeros(6)
    DC = np.zeros(6)

    NextDerivs = [np.zeros(4)] * 6
    Vals = np.zeros(4)
    Derivs = np.zeros(3)

    def __init__(self):
        self.A = np.array([1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0])
        self.B = np.array([[1.0 / 5.0, 0, 0, 0, 0],
                           [3.0 / 40.0, 9.0 / 40.0, 0, 0, 0],
                           [3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0, 0, 0],
                           [-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0, 0],
                           [1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0]])
        self.C = [37.0 / 378.0, 0, 250.0 / 621.0, 125.0 / 594.0, 0, 512.0 / 1771.0]
        self.DC = [37.0 / 378.0 - 2825.0 / 27648.0,
                   0,
                   250.0 / 621.0 - 18575.0 / 48384.0,
                   125.0 / 594.0 - 13525.0 / 55296.0,
                   -277.0 / 14336.0,
                   512.0 / 1771.0 - 1.0 / 4.0]

    def ComputeAStep(self, xprev, dxprev, xnext, t, delT, field: StructureGrid):
        delTActual = 0
        numDerivs = 3
        numVals = numDerivs + 1
        error = 0.0
        for i in range(numVals - 1):
            self.Vals[i] = xprev[i]
        self.Vals[numVals - 1] = t
        if dxprev:
            for i in range(numDerivs):
                self.NextDerivs[0][i] = dxprev[i]
        elif not field.sample(self.Vals, self.NextDerivs[0], useNormalize=True):
            for i in range(numVals - 1):
                xnext[i] = self.Vals[i]
            return ReasonForTermination.OUT_OF_DOMAIN, delTActual, error

        summary = 0.0
        for i in range(1, 6):
            for j in range(numVals - 1):
                summary = 0
                for k in range(i):
                    summary += self.B[i - 1][k] * self.NextDerivs[k][j]
                self.Vals[j] = xprev[j] + delT * summary
            self.Vals[numVals - 1] = t + delT * self.A[i - 1]
            if not field.sample(self.Vals, self.NextDerivs[0], useNormalize=True):
                for l in range(numVals - 1):
                    xnext[l] = self.Vals[l]
                delTActual = delT * self.A[i - 1]
                return ReasonForTermination.OUT_OF_DOMAIN, delTActual, error
        for i in range(numDerivs):
            summary = 0
            for j in range(6):
                summary += self.C[j] * self.NextDerivs[j][i]
            xnext[i] = xprev[i] + delT * summary
        delTActual = delT

        for i in range(numDerivs):
            summary = 0
            for j in range(6):
                summary += self.DC[j] * self.NextDerivs[j][i]
            error += delT * summary * delT * summary
        error = math.sqrt(error)

        numZero = 0
        for i in range(numDerivs):
            if xnext[i] == xprev[i]:
                numZero += 1
        if numZero == numDerivs:
            return ReasonForTermination.UNEXPECTED_VALUE, delTActual, error
        return 0, delTActual, error

    def ComputeNextStep(self, xprev, dxprev, xnext, t, delT, minStep, maxStep, maxError, field: StructureGrid):
        estError = float("inf")
        if minStep < 0:
            minStep = -minStep
        if maxStep < 0:
            maxStep = -maxStep

        delTActual = 0

        absDT = abs(delT)
        if ((minStep == absDT) and (maxStep == absDT)) or (maxError <= 0.0):
            retVal, delTActual, estError = self.ComputeAStep(xprev, dxprev, xnext, t, delT, field)
            return retVal, delTActual
        elif minStep > maxStep:
            return ReasonForTermination.UNEXPECTED_VALUE, delTActual

        while estError > maxError:
            retVal, delTActual, estError = self.ComputeAStep(xprev, dxprev, xnext, t, delT, field)
            if retVal:
                return retVal, delTActual
            absDT = abs(delT)
            if absDT == minStep:
                break
            errRatio = estError / maxError
            if errRatio == 0.0:
                tmp = -minStep if delT < 0 else minStep
            elif errRatio > 1:
                tmp = 0.9 * delT * math.pow(errRatio, -0.25)
            else:
                tmp = 0.9 * delT * math.pow(errRatio, -0.2)
            tmp2 = abs(tmp)
            shouldBreak = False
            if tmp2 > maxStep:
                delT = maxStep * delT / abs(delT)
                shouldBreak = True
            elif tmp2 < minStep:
                delT = minStep * delT / abs(delT)
                shouldBreak = True
            else:
                delT = tmp

            tmp2 = t + delT
            if tmp2 == t:
                return ReasonForTermination.UNEXPECTED_VALUE
            if shouldBreak:
                retVal, delTActual, estError = self.ComputeAStep(xprev, dxprev, xnext, t, delT, field)
                if retVal:
                    return retVal, delTActual
                break
        return 0, delTActual


class Line:
    p1 = 0  # point's index
    p2 = 0

    def __init__(self, p1: int, p2: int):
        self.p1 = p1
        self.p2 = p2


class Units(Enum):
    LENGTH_UNIT = 1
    CELL_LENGTH_UNIT = 2


class ReasonForTermination(Enum):
    OUT_OF_DOMAIN = 1
    NOT_INITIALIZED = 2
    UNEXPECTED_VALUE = 3
    OUT_OF_LENGTH = 4
    OUT_OF_STEPS = 5
    STAGNATION = 6
    FIXED_REASONS_FOR_TERMINATION_COUNT = 7


class Direction(Enum):
    FORWARD = 1
    BACKWARD = 2
    BOTH = 3


def ConvertToLength(interval, unit, cellLength):
    retValue = 0.0
    if unit == Units.LENGTH_UNIT:
        retValue = interval
    elif unit == Units.CELL_LENGTH_UNIT:
        retValue = interval * cellLength
    return retValue


class StreamTracer:
    points = []
    polyLines = []
    pointData = {}
    Seeds = []
    eps = 1e-7
    # StartPos = [0.0, 0.0, 0.0]
    TerminalSpeed = 0
    LastUsedStepSize = 0.0
    MaximumPropagation = 1.0
    MinimumIntegrationStep = 0.01
    MaximumIntegrationStep = 1.0
    InitialIntegrationStep = 0.5
    IntegrationStepUnit = Units.CELL_LENGTH_UNIT
    IntegrationDirection = Direction.FORWARD
    Integrator = None
    MaximumError = 1.0e-6
    MaximumNumberOfStep = 2000

    ComputeVorticity = False
    RotationScale = 1.0

    field = None

    def __init__(self,
                 field: StructureGrid,
                 integrator: IntegratorRK45,
                 seeds=None,
                 terminalSpeed=1.0e-12,
                 maximumPropagation=1.0,
                 minimumIntegrationStep=0.01,
                 maximumIntegrationStep=0.5,
                 initialIntegrationStep=0.2,
                 integrationStepUnit=Units.CELL_LENGTH_UNIT,
                 integrationDirection=Direction.FORWARD,
                 maximumError=1.0e-6,
                 maximumNumberOfStep=2000,
                 computeVorticity=False,
                 rotationScale=1.0
                 ):
        assert seeds is not None
        self.field = field
        self.Integrator = integrator
        self.Seeds = seeds
        self.TerminalSpeed = terminalSpeed
        self.MaximumPropagation = maximumPropagation
        self.MinimumIntegrationStep = minimumIntegrationStep
        self.MaximumIntegrationStep = maximumIntegrationStep
        self.InitialIntegrationStep = initialIntegrationStep
        self.IntegrationStepUnit = integrationStepUnit
        self.IntegrationDirection = integrationDirection
        self.MaximumError = maximumError
        self.MaximumNumberOfStep = maximumNumberOfStep
        self.ComputeVorticity = computeVorticity
        self.RotationScale = rotationScale
        for array in self.field.pointData.keys():
            self.pointData[array] = []

    def ConvertIntervals(self, direction, cellLength):
        minStep = maxStep = step = \
            direction * ConvertToLength(self.InitialIntegrationStep, self.IntegrationStepUnit, cellLength)
        if self.MinimumIntegrationStep > 0.0:
            minStep = ConvertToLength(self.MinimumIntegrationStep, self.IntegrationStepUnit, cellLength)
        if self.MaximumIntegrationStep > 0.0:
            maxStep = ConvertToLength(self.MaximumIntegrationStep, self.IntegrationStepUnit, cellLength)
        return step, minStep, maxStep

    def addLine(self, points, pointData):
        origin_pts = len(self.points)
        polyLine = [i + origin_pts for i in range(len(points))]
        self.points.extend(points)
        self.pointData["velocity"].extend(pointData)
        self.polyLines.append(polyLine)

    def integrate(self):
        assert self.Integrator is not None
        forward_points = []
        forward_pointData = []
        backward_points = []
        backward_pointData = []
        for i, seed in enumerate(self.Seeds):
            forward_points.clear()
            forward_pointData.clear()
            backward_points.clear()
            backward_points.clear()
            if self.IntegrationDirection == Direction.FORWARD:
                self.integrate_dir(seed, forward_points, forward_pointData, direction=1)
                if len(forward_points) > 1:
                    self.addLine(forward_points, forward_pointData)
            elif self.IntegrationDirection == Direction.BACKWARD:
                self.integrate_dir(seed, backward_points, backward_pointData, direction=-1)
                if len(backward_points) > 1:
                    self.addLine(backward_points, backward_pointData)
            elif self.IntegrationDirection == Direction.BOTH:
                self.integrate_dir(seed, forward_points, forward_pointData, direction=1)
                self.integrate_dir(seed, backward_points, backward_pointData, direction=-1)
                if len(forward_points) > 1:
                    self.addLine(forward_points, forward_pointData)
                if len(backward_points) > 1:
                    self.addLine(backward_points, backward_pointData)

    def integrate_dir(self, start: np.array, points: [], pointData: [], direction=1):
        propagation = 0
        numSteps = 0
        integrationTime = 0
        retVal = 0
        # velocity = [0.0, 0.0, 0.0]
        # for seed in self.Seeds:
        point1 = deepcopy(start)
        point2 = deepcopy(point1)
        velocity = [0.0, 0.0, 0.0]
        # pcoords = [0.0, 0.0, 0.0]
        # vort = [0.0, 0.0, 0.0]
        # omega = 0.0
        # index = 0
        cell = self.field.getCell(point1)
        if not self.field.sample(point1, velocity, cell):
            retVal = ReasonForTermination.OUT_OF_DOMAIN
            return
        points.append(deepcopy(point1))
        pointData.append(deepcopy(velocity))

        stepSize = 0
        minStep = 0
        maxStep = 0

        cellLength = math.sqrt(cell.GetLength2())
        speed = np.linalg.norm(velocity)
        if speed != 0.0:
            stepSize, minStep, maxStep = self.ConvertIntervals(direction, cellLength)
        if self.ComputeVorticity:
            print("TODO")  # TODO
        while propagation < self.MaximumPropagation:
            numSteps += 1
            if numSteps > self.MaximumNumberOfStep:
                break
            if speed == 0 or speed <= self.TerminalSpeed:
                break
            aStep = abs(stepSize)
            if propagation + aStep > self.MaximumPropagation:
                aStep = self.MaximumPropagation
                if stepSize >= 0:
                    stepSize = ConvertToLength(aStep, Units.LENGTH_UNIT, cellLength)
                else:
                    stepSize = ConvertToLength(aStep, Units.LENGTH_UNIT, cellLength) * (-1.0)
                maxStep = stepSize
            tmp, stepTaken = self.Integrator.ComputeNextStep(xprev=point1, dxprev=None, xnext=point2, t=0,
                                                             delT=stepSize, minStep=minStep, maxStep=maxStep,
                                                             maxError=self.MaximumError, field=self.field)
            if tmp != 0:
                retVal = tmp
                break
            for i in range(3):
                point1[i] = point2[i]

            cell = self.field.getCell(point2)
            if not self.field.sample(point2, velocity, cell):
                retVal = ReasonForTermination.OUT_OF_DOMAIN
                break

            speed2 = np.linalg.norm(velocity)
            if (speed2 + speed) / 2 <= self.TerminalSpeed:
                retVal = ReasonForTermination.STAGNATION
                break

            integrationTime += stepTaken / speed
            propagation += abs(stepSize)
            cellLength = math.sqrt(cell.GetLength2())
            speed = speed2

            # convertedPoint = np.zeros(3)
            # for i in range(3):
            #     convertedPoint[i] = point1[i]
            if abs(points[-1][0] - point1[0]) > self.eps or \
                    abs(points[-1][1] - point1[1]) > self.eps or \
                    abs(points[-1][2] - point1[2]) > self.eps:
                points.append(deepcopy(point1))
                pointData.append(deepcopy(velocity))
                # self.lines.append(Line(len(self.points) - 2, len(self.points) - 1))
                if self.ComputeVorticity:
                    print("TODO")  # TODO
            if speed == 0 or speed <= self.TerminalSpeed:
                retVal = ReasonForTermination.STAGNATION
                break
            step, minStep, maxStep = self.ConvertIntervals(direction, cellLength)
            if abs(stepSize) < abs(minStep):
                stepSize = abs(minStep) * stepSize / abs(stepSize)
            elif abs(stepSize) > abs(maxStep):
                stepSize = abs(maxStep) * stepSize / abs(stepSize)
        return retVal


if __name__ == "__main__":
    nc_file = r"E:\\BaiduNetdiskDownload\\ocean_data\\ECS_SCS_WYJ_Coarse_skp_hcst20103_1242.nc"
    ds = xarray.open_dataset(nc_file)
    print(ds)

    field = StructureGrid(ds)

    seeds = GenerateSeeds(field, num=100)

    # with open("random_seed.csv", "w") as f:
    #     writer = csv.writer(f)
    #     writer.writerow(["x", "y", "z"])
    #     writer.writerows(seeds)

    streamTracer = StreamTracer(field, IntegratorRK45(), seeds,
                                integrationDirection=Direction.BOTH,
                                integrationStepUnit=Units.CELL_LENGTH_UNIT,
                                initialIntegrationStep=0.2,
                                minimumIntegrationStep=0.01,
                                maximumIntegrationStep=0.5,
                                maximumNumberOfStep=2000,
                                maximumPropagation=31,
                                terminalSpeed=1e-12,
                                maximumError=1e-6)
    print("Integrate Start!")
    streamTracer.integrate()
    print("Integrate Finished!")

    streamLines = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    # points.SetNumberOfPoints(len(streamTracer.points))
    for p in streamTracer.points:
        points.InsertNextPoint(p)
    streamLines.SetPoints(points)
    for l in streamTracer.polyLines:
        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(len(l))
        for i in range(len(l)):
            polyLine.GetPointIds().SetId(i, l[i])
        streamLines.InsertNextCell(polyLine.GetCellType(), polyLine.GetPointIds())
    pointData = vtk.vtkFloatArray()
    pointData.SetName("velocity")
    pointData.SetNumberOfTuples(len(streamTracer.pointData["velocity"]))
    pointData.SetNumberOfComponents(3)
    for i, v in enumerate(streamTracer.pointData["velocity"]):
        pointData.InsertTuple3(i, v[0], v[1], v[2])
    streamLines.GetPointData().AddArray(pointData)

    writer = vtk.vtkDataSetWriter()
    writer.SetFileName("test_stream.vtk")
    writer.SetInputData(streamLines)
    writer.Write()
