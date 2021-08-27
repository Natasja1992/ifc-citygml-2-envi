from osgeo.osr import SpatialReference, CoordinateTransformation
import numpy as np


class Transform(object):
    def __init__(self):
        pass

    def transform_point(self, coord):
        raise NotImplementedError

    def transform_points(self, coord_list):
        assert type(coord_list) is list, "expected a list of coordinates"

        return [self.transform_point(c) for c in coord_list]

    def transform_building(self, building):
        assert type(building) is list, "expected a list of polygons in building"

        return [self.transform_points(polygon) for polygon in building]


class TransformEPSG2EPSG(Transform):
    def __init__(self, source_srs_number, target_srs_number):
        super().__init__()
        source_srs = SpatialReference()
        source_srs.ImportFromEPSG(source_srs_number)
        target_srs = SpatialReference()
        target_srs.ImportFromEPSG(target_srs_number)

        self._source2target = CoordinateTransformation(source_srs, target_srs)

    # transform point from surroundings to target srs
    def transform_point(self, coord):
        assert type(coord) is tuple, "expected a tuple as coordinate"
        assert len(coord) == 3, "expected 3 elements in coordinate"
        for v in coord:
            assert type(v) is float, "expected type float in elements of coordinate"

        return self._source2target.TransformPoint(*coord)


class TransformCustom2EPSG(Transform):
    # with reference point as longitude, lattitude, elevation in decimals and true north in degrees
    def __init__(self, ref_point_in_target, true_north):
        super().__init__()

        # define transform matrix based on true north, units and reference point
        mm2m_factor = 0.001

        self._T = np.identity(4)

        i = 0
        while i < 3:
            self._T[i][3] = ref_point_in_target[i]
            i += 1

        rotation_rad = np.radians(true_north)
        self._T[0][0] = np.cos(rotation_rad) * mm2m_factor
        self._T[0][1] = -np.sin(rotation_rad) * mm2m_factor
        self._T[1][0] = np.sin(rotation_rad) * mm2m_factor
        self._T[1][1] = np.cos(rotation_rad) * mm2m_factor
        self._T[2][2] = mm2m_factor

    def transform_point(self, coord):
        assert type(coord) is tuple, "expected a tuple as coordinate"
        assert len(coord) == 3, "expected 3 elements in coordinate"
        for v in coord:
            assert type(v) is float, "expected type float in elements of coordinate"

        point = np.array([coord[0], coord[1], coord[2], 1])
        transformed_point = self._T.dot(point)
        return transformed_point[0], transformed_point[1], transformed_point[2]


class TransformEPSG2Custom(Transform):
    # assuming reference point srs is same as target srs (which is the only case necessary in this tool)
    # with reference point as longitude, lattitude, elevation in decimals and true north in degrees
    # and resolution is resolution of custom system in meters,
    def __init__(self, ref_point_in_target, true_north_target, resolution):
        super().__init__()

        # define transform matrix based on true north, units and reference point
        m2voxel_factor_x = resolution[0]
        m2voxel_factor_y = resolution[1]
        m2voxel_factor_z = resolution[2]

        # generate 4x4 identity matrix
        self._T = np.identity(4)

        # add transpose to matrix
        i = 0
        while i < 3:
            self._T[i][3] = ref_point_in_target[i]
            i += 1

        # add rotation and scale factor to matrix
        rotation_rad = np.radians(true_north_target)
        self._T[0][0] = np.cos(rotation_rad) * m2voxel_factor_x
        self._T[0][1] = -np.sin(rotation_rad) * m2voxel_factor_y
        self._T[1][0] = np.sin(rotation_rad) * m2voxel_factor_x
        self._T[1][1] = np.cos(rotation_rad) * m2voxel_factor_y
        self._T[2][2] = m2voxel_factor_z

        self._T_inv = np.linalg.inv(self._T)

    def transform_point(self, coord):
        assert type(coord) is tuple, "expected a tuple as coordinate"
        assert len(coord) == 3, "expected 3 elements in coordinate"
        for v in coord:
            assert type(v) is float or type(v) is np.float64, "expected type float in elements of coordinate"

        point = np.array([coord[0], coord[1], coord[2], 1])
        transformed_point = self._T_inv.dot(point)
        return transformed_point[0], transformed_point[1], transformed_point[2]

    def transform_point_reverse(self, coord):
        assert type(coord) is tuple, "expected a tuple as coordinate"
        assert len(coord) == 3, "expected 3 elements in coordinate"
        for v in coord:
            assert type(v) is float or type(v) is np.float64, "expected type float in elements of coordinate"

        point = np.array([coord[0], coord[1], coord[2], 1])
        transformed_point = self._T.dot(point)
        return transformed_point[0], transformed_point[1], transformed_point[2]
