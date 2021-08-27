from ifcParser import IfcParser
from citygmlParser import CitygmlParser
import georeferencing
from elements import *
from shapely.geometry import Point, Polygon, LineString
from shapely import speedups
import reverse_geocode
from timezonefinder import TimezoneFinder

import os
import math
import numpy as np

speedups.disable()


class ConversionTool(object):
    def __init__(self, input_ifc, input_citygml, variables):
        # initiate classes
        print('loading ifc [%s]...' % input_ifc)
        self._ifc = IfcParser(input_ifc)
        print('loading citygml [%s]...' % input_citygml)
        self._cg = CitygmlParser(input_citygml)
        print('input files loaded')
        self._voxel_grid = VoxelGrid()

        print('calculation parameters...')

        # variables and values
        self.variables = variables
        self._orientation = self._ifc.model_orientation()
        self._ref_point_ifc = self._ifc.reference_point()
        self._bbox_citygml = self._cg.get_bbox()

        self._srs = self._bbox_citygml[2]
        assert self._srs.split(':')[0].lower() == 'epsg', "expected epsg as srs"
        self.variables["reference_system_EPSG"] = int(self._srs.split(':')[1])

        ## size of bbox site in target data structure, in voxels and in meters ceiled based on grid size ##

        # calculate size of bbox site in target data structure, ceiled based on grid size
        (res_x, res_y, res_z) = self.variables["resolution"]
        ((ifc_min_x, ifc_min_y, ifc_min_z), (ifc_max_x, ifc_max_y, ifc_max_z)) = self._ifc.bbox_site()
        length_x = (ifc_max_x - ifc_min_x) / 1000.0
        length_y = (ifc_max_y - ifc_min_y) / 1000.0
        length_x_voxels = math.ceil(length_x / res_x)
        length_y_voxels = math.ceil(length_y / res_y)
        length_x_resceil = length_x_voxels * res_x
        length_y_resceil = length_y_voxels * res_y

        # ceil the borderwidth based on the resolution (grid size)
        borderwidth = self.variables["citygml_borderwidth"]
        borderwidth_x_voxels = math.ceil(borderwidth / res_x)
        borderwidth_y_voxels = math.ceil(borderwidth / res_y)
        borderwidth_x_resceil = borderwidth_x_voxels * res_x
        borderwidth_y_resceil = borderwidth_y_voxels * res_y

        # calculate size of border grid in meters
        border_grid = self.variables["border_grid"]
        border_grid_x = border_grid * res_x
        border_grid_y = border_grid * res_y

        # generate bounding boxes/rectangles of site, surrounding and envi-met model (x and y) in METERS
        self.bbox_surrounding = ((border_grid_x, border_grid_y),
                                 (border_grid_x + 2.0 * borderwidth_x_resceil + length_x_resceil,
                                  border_grid_y + 2.0 * borderwidth_y_resceil + length_y_resceil))
        self.bbox_site = ((border_grid_x + borderwidth_x_resceil, border_grid_y + borderwidth_y_resceil),
                          (border_grid_x + borderwidth_x_resceil + length_x_resceil,
                           border_grid_y + borderwidth_y_resceil + length_y_resceil))
        self.bbox_envi = ((0, 0), (2.0 * border_grid_x + 2.0 * borderwidth_x_resceil + length_x_resceil,
                                   2.0 * border_grid_y + 2.0 * borderwidth_y_resceil + length_y_resceil))

        # generate bounding boxes/rectangles of site, surrounding and envi-met model (x and y) in VOXELS
        self.bbox_surrounding_v = ((float(border_grid), float(border_grid), 0.0),
                                   (border_grid + 2.0 * borderwidth_x_voxels + length_x_voxels,
                                    border_grid + 2.0 * borderwidth_y_voxels + length_y_voxels, 0.0))  #todo: z
        self.bbox_site_v = ((border_grid + borderwidth_x_voxels, border_grid + borderwidth_y_voxels),
                          (border_grid + borderwidth_x_voxels + length_x_voxels,
                           border_grid + borderwidth_y_voxels + length_y_voxels))
        self.bbox_envi_v = ((0.0, 0.0), (2.0 * border_grid + 2.0 * borderwidth_x_voxels + length_x_voxels,
                                   2.0 * border_grid + 2.0 * borderwidth_y_voxels + length_y_voxels))

        # calculate real world location of 0,0 in target (voxelized) data structure
        m2mm = 1000.0
        z = ifc_min_z   #todo: compare with min_z citygml, get lowest point in total model (and floor based on grid)
        ref_target_in_ifc = (ifc_min_x - borderwidth_x_resceil * m2mm - border_grid_x * m2mm,
                             ifc_min_y - borderwidth_y_resceil * m2mm - border_grid_y * m2mm,
                             z)

        # create transformations
        ifc2voxel = georeferencing.TransformEPSG2Custom(ref_target_in_ifc, 0,
                                                        tuple(m2mm * x for x in self.variables['resolution']))

        wgs842target = georeferencing.TransformEPSG2EPSG(4326, self.variables["reference_system_EPSG"])
        ref_point_switched = (self._ref_point_ifc[1], self._ref_point_ifc[0], self._ref_point_ifc[2])  # switched lon lat because of library change
        self._ifc_ref_point_in_target = wgs842target.transform_point(ref_point_switched)
        ifc2refsrs = georeferencing.TransformCustom2EPSG(self._ifc_ref_point_in_target, -self._orientation)
        self.ref_point_target = ifc2refsrs.transform_point(ref_target_in_ifc)   # in target reference system

        cg2voxel = georeferencing.TransformEPSG2Custom(self.ref_point_target, self._orientation, self.variables['resolution'])

        # real world coordinates bbox/rectangle surroundings (4 points because of rotation)
        self.define_realworld_bboxes(cg2voxel)

        # check if ifc model is inside bounding box citygml
        assert self.ref_point_in_bbox(), "ifc model not within citygml model boundries"

        # extract building from input files
        print('extracting buildings from input files...')
        cg_buildings = self._cg.get_buildings()
        ifc_slabs = self._ifc.get_slabs()

        # create building number counter and convert IFC building
        print('converting IFC building...')
        building_no = 1
        self.convert_ifc_building(ifc2voxel, ifc_slabs, building_no)

        # convert CityGML building
        building_no += 1
        self.convert_cg_buildings(cg2voxel, cg_buildings, building_no)

        # extracting and converting terrain model from CityGML
        print()
        print('extracting and generating dem...')
        self.convert_terrain(cg2voxel)

        # extracting and converting trees from CityGML
        print('extracting and converting trees...')
        self.convert_trees(cg2voxel)

    # calculate bounding boxes with real world coordinates for surrounding area and envi-met model area
    def define_realworld_bboxes(self, cg2voxel):
        self._sur_bbox_rw = (cg2voxel.transform_point_reverse(self.bbox_surrounding_v[0]),
                             cg2voxel.transform_point_reverse(
                                 (self.bbox_surrounding_v[1][0], self.bbox_surrounding_v[0][1], 0.0)),
                             cg2voxel.transform_point_reverse(self.bbox_surrounding_v[1]),
                             cg2voxel.transform_point_reverse(
                                 (self.bbox_surrounding_v[0][0], self.bbox_surrounding_v[1][1], 0.0)))
        self._envi_bbox_rw = (cg2voxel.transform_point_reverse((self.bbox_envi_v[0][0], self.bbox_envi_v[0][1], 0.0)),
                              cg2voxel.transform_point_reverse((self.bbox_envi_v[1][0], self.bbox_envi_v[0][1], 0.0)),
                              cg2voxel.transform_point_reverse((self.bbox_envi_v[1][0], self.bbox_envi_v[1][1], 0.0)),
                              cg2voxel.transform_point_reverse((self.bbox_envi_v[0][0], self.bbox_envi_v[1][1], 0.0)))

    # convert buildings from citygml within area of interest to custom coordinate system and update corresponding voxels
    # in the voxel grid
    def convert_cg_buildings(self, cg2voxel, cg_buildings, building_no):
        buildings = []
        for building in cg_buildings:
            building_found = False
            for polygon in building:
                for point in polygon:
                    # check if building within area of interest
                    if self.point_in_bbox(self._sur_bbox_rw, point):
                        buildings.append(Building(cg2voxel.transform_building(building), building_no, 'citygml'))
                        building_found = True
                        building_no += 1
                        break
                if building_found:
                    break
        cg_building_ctr = 1
        for building in buildings:
            print('\rconverting CityGML building %d of %d...' % (cg_building_ctr, len(buildings)), end="")
            cg_building_ctr += 1
            geom = building.get_geom()
            bbox_building = building.get_bbox()
            bbox = [list(c) for c in bbox_building]
            for axis in range(2):
                bbox[0][axis] = max(bbox[0][axis], int(self.bbox_surrounding_v[0][axis]))
                bbox[1][axis] = min(bbox[1][axis], int(self.bbox_surrounding_v[1][axis]))

            for x in range(bbox[0][0], bbox[1][0] + 1):
                for y in range(bbox[0][1], bbox[1][1] + 1):
                    for z in range(bbox[0][2], bbox[1][2] + 1):
                        voxel = (x, y, z)
                        line = [(bbox_building[0][0] - 1, y + 0.5, z + 0.5), (x + 0.5, y + 0.5, z + 0.5)]
                        count = 0
                        for polygon in geom:
                            if self.line_crosses_polygon(line, polygon) is not None:
                                count += 1
                        if count % 2 == 1:
                            self._voxel_grid.get_voxel(voxel).set_building(building, 100)

    # convert building from ifc to custom coordinate system and update corresponding voxels in the voxel grid
    def convert_ifc_building(self, ifc2voxel, slabs, building_no):
        ifc_building = Building(ifc2voxel.transform_building(slabs), building_no, 'ifc')
        geom = ifc_building.get_geom()
        bbox_building = ifc_building.get_bbox()
        bbox = [list(c) for c in bbox_building]
        for axis in range(2):
            bbox[0][axis] = max(bbox[0][axis], int(self.bbox_surrounding_v[0][axis]))
            bbox[1][axis] = min(bbox[1][axis], int(self.bbox_surrounding_v[1][axis]))
        for x in range(bbox[0][0], bbox[1][0] + 1):
            for y in range(bbox[0][1], bbox[1][1] + 1):
                line = [(x + 0.5, y + 0.5, bbox_building[0][2] - 1), (x + 0.5, y + 0.5, bbox_building[1][2] + 1)]

                collision_points = []
                for polygon in geom:
                    collision_point = self.line_crosses_polygon(line, polygon)
                    if collision_point is not None:
                        collision_points.append(collision_point)

                if len(collision_points) >= 2:  # found collisions, adding voxels
                    collision_z = [p[2] for p in collision_points]
                    min_collision = min(collision_z)
                    max_collision = max(collision_z)

                    for z in range(bbox[0][2], bbox[1][2] + 1):
                        if min_collision <= z + 0.5 <= max_collision:
                            voxel = (x, y, z)
                            self._voxel_grid.get_voxel(voxel).set_building(ifc_building, 100)

    # convert terrain from citygml to custom coordinate system and update corresponding voxels in the voxel grid
    def convert_terrain(self, cg2voxel):
        filtered_tin = []
        tin_list = self._cg.get_tin()
        for triangle in tin_list:
            for point in triangle:
                if self.point_in_bbox(self._envi_bbox_rw, point):
                    filtered_tin.append(cg2voxel.transform_points(triangle))
                    break
        for triangle in filtered_tin:
            bbox = [[int(math.floor(min([c[0] for c in triangle]))),
                     int(math.floor(min([c[1] for c in triangle]))),
                     int(math.floor(min([c[2] for c in triangle])))],
                    [int(math.ceil(max([c[0] for c in triangle]))),
                     int(math.ceil(max([c[1] for c in triangle]))),
                     int(math.ceil(max([c[2] for c in triangle])))]]

            for axis in range(2):
                bbox[0][axis] = max(bbox[0][axis], int(self.bbox_envi_v[0][axis]))
                bbox[1][axis] = min(bbox[1][axis], int(self.bbox_envi_v[1][axis]))

            for x in range(bbox[0][0], bbox[1][0]):
                for y in range(bbox[0][1], bbox[1][1]):
                    line = [(x + 0.5, y + 0.5, bbox[0][2] - 1), (x + 0.5, y + 0.5, bbox[1][2] + 1)]

                    collision_point = self.line_crosses_polygon(line, triangle)
                    if collision_point is not None:
                        for z in range(0, int(round(collision_point[2]))):
                            voxel = (x, y, z)
                            self._voxel_grid.get_voxel(voxel).set_terrain()

    # convert trees from citygml to custom coordinate system and update corresponding voxels in the voxel grid
    def convert_trees(self, cg2voxel):
        tree_instances = []
        trees = self._cg.get_trees()
        for tree in trees:
            if self.point_in_bbox(self._sur_bbox_rw, tree[0]):
                tree_instances.append(Tree(cg2voxel.transform_point(tree[0]), tree[1], 'citygml'))
        for tree in tree_instances:
            if tree.get_height() > 1.5:
                root_voxel = tree.get_root_voxel()
                x, y, min_z = root_voxel
                max_z = int(math.ceil(min_z + tree.get_height() / self.variables["resolution"][2]))
                has_building = False
                for z in range(min_z, max_z + 2):
                    if self._voxel_grid.get_voxel((x, y, z)).has_building():
                        has_building = True
                        break
                if not has_building:
                    self._voxel_grid.get_voxel(root_voxel).set_tree(tree)

    # functions within conversion tool class #

    # check whether ifc input model lies within boundaries of citygml input model
    def ref_point_in_bbox(self):
        ifc_x, ifc_y, ifc_z = self._ifc_ref_point_in_target
        cg_lc, cg_uc, _ = self._bbox_citygml
        return cg_lc[0] < ifc_x < cg_uc[0] and cg_lc[1] < ifc_y < cg_uc[1]

    # check if point lies within bounding box (2d - xy)
    def point_in_bbox(self, bbox_rect, point):
        assert type(bbox_rect) == tuple, "expected tuple for bbox"
        assert len(bbox_rect) == 4, "expected 4 coords for bbox"
        for coord in bbox_rect:
            assert type(coord) == tuple, "expected tuple for coord in bbox"

        p1, p2, p3, p4 = bbox_rect
        pnt = Point(point[0], point[1])
        polygon = Polygon([(p1[0], p1[1]), (p2[0], p2[1]), (p3[0], p3[1]), (p4[0], p4[1])])
        return polygon.contains(pnt)

    # return collision point between line and polygon if they cross, otherwise return none
    def line_crosses_polygon(self, line, polygon):
        # define ray trough line (as point and direction)
        ray_point = np.array(line[1])
        if line[0][0] != line[1][0]:
            ray_direction = np.array((-1, 0, 0))
        elif line[0][1] != line[1][1]:
            ray_direction = np.array((0, -1, 0))
        elif line[0][2] != line[1][2]:
            ray_direction = np.array((0, 0, -1))
        else:
            raise AssertionError("invalid line")

        # define plane through polygon (as point and normal)
        x0, y0, z0 = polygon[0]
        x1, y1, z1 = polygon[1]
        x2, y2, z2 = polygon[2]

        ux, uy, uz = [x1-x0, y1-y0, z1-z0]
        vx, vy, vz = [x2-x0, y2-y0, z2-z0]

        u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]

        plane_point = np.array(polygon[0])
        plane_normal = np.array(u_cross_v)

        # where do ray on line and plane on polygon collide
        epsilon = 1e-6
        ndotu = plane_normal.dot(ray_direction)
        if abs(ndotu) < epsilon:
            return None

        w = ray_point - plane_point
        si = -plane_normal.dot(w) / ndotu
        collision_point = w + si * ray_direction + plane_point

        # check if collision point is on line and in polygon
        pnt = Point(tuple(collision_point))
        line = LineString(line)
        poly = Polygon(polygon)

        # when polygon perpendicular to z-axis: use projection on y-z plane (otherwise just uses xy-projection)
        if plane_normal[2] != 0:
            pnt_in_poly = poly.contains(pnt)
        else:
            poly_proj_list = []
            for point in polygon:
                poly_proj_list.append((point[1], point[2]))
            pnt_proj = Point((collision_point[1], collision_point[2]))
            poly_proj = Polygon(poly_proj_list)
            pnt_in_poly = poly_proj.contains(pnt_proj)
        if line.distance(pnt) < epsilon and pnt_in_poly:
            return collision_point
        else:
            return None

    # functions called from other classes #

    # return voxel grid
    def get_voxel_grid(self):
        return self._voxel_grid

    # create and return model description based on input file names
    def get_model_description(self):
        ifc_file = os.path.basename(self._ifc.filename)
        cg_file = os.path.basename(self._cg.filename)
        return "Area input file generated based on the following IFC and CityGML model: " \
               "IFC-file: {}, CityGML-file: {}".format(ifc_file, cg_file)
        # todo: multiple files and name/description,  basic metadata (lod,etc.)

    # return number of grid cells / voxels in x and y direction
    def get_grids(self):
        (res_x, res_y, res_z) = self.variables["resolution"]
        ((ifc_min_x, ifc_min_y, ifc_min_z), (ifc_max_x, ifc_max_y, ifc_max_z)) = self._ifc.bbox_site()

        # calculate number of grid cells in I and J direction
        borderwidth = self.variables["citygml_borderwidth"]
        border_grid = self.variables["border_grid"]
        length_x = (ifc_max_x - ifc_min_x)/1000.0 + 2.0 * borderwidth
        length_y = (ifc_max_y - ifc_min_y)/1000.0 + 2.0 * borderwidth
        grids_x = math.ceil(length_x / res_x) + 2 * border_grid
        grids_y = math.ceil(length_y / res_y) + 2 * border_grid

        return grids_x, grids_y

    # return model orientation
    def get_model_orientation(self):
        return self._orientation

    # return reference point from ifc, expressed in citygml coordinate system
    def get_reference_point(self):
        return self.ref_point_target

    # return longitude, latitude and elevation from ifc reference point (wgs84)
    def get_lonlatel(self):
        return self._ref_point_ifc

    # lookup location name based on longitude and latitude
    def get_location_name(self):
        try:
            coordinate = [(self._ref_point_ifc[1], self._ref_point_ifc[0])]
            result = reverse_geocode.search(coordinate)[0]
        except Exception:
            result = {'city': 'unknown', 'country': 'unknown', 'country_code': 'unknown'}
        return result

    # lookup timezone based on longitude and latitude
    def get_timezone(self):
        tf = TimezoneFinder()
        latitude, longitude = self._ref_point_ifc[1], self._ref_point_ifc[0]
        tz = tf.timezone_at(lng=longitude, lat=latitude)
        return tz

