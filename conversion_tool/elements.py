import math


# Tree definition
class Tree(object):
    def __init__(self, root_point, height, source):
        self._root_point = root_point
        self._height = height
        self._source = source

    # Return the root voxel of the tree
    def get_root_voxel(self):
        return (int(math.floor(self._root_point[0])), int(math.floor(self._root_point[1])), int(math.floor(self._root_point[2])))

    # Return the height of the tree
    def get_height(self):
        return self._height


# Building definition
class Building(object):
    def __init__(self, building, nr, building_source):
        self._building_geom = building
        self._building_nr = nr
        self._building_source = building_source
        self._bbox_building = self.bbox_building(self._building_geom)

    # Return the building geometry
    def get_geom(self):
        return self._building_geom

    # Return the bounding box of the building
    def get_bbox(self):
        return self._bbox_building

    # Return the building number
    def get_nr(self):
        return self._building_nr

    # Calculates de bounding box of a building in the voxel system, round to whole voxels
    def bbox_building(self, geom):
        min_x = float('inf')
        max_x = float('-inf')
        min_y = float('inf')
        max_y = float('-inf')
        min_z = float('inf')
        max_z = float('-inf')

        for polygon in geom:
            for point in polygon:
                assert len(point) == 3, "expected x, y and z for coordinate"
                min_x = min(min_x, point[0])
                max_x = max(max_x, point[0])
                min_y = min(min_y, point[1])
                max_y = max(max_y, point[1])
                min_z = min(min_z, point[2])
                max_z = max(max_z, point[2])

        return (math.floor(min_x), math.floor(min_y), math.floor(min_z)), \
               (math.floor(max_x), math.floor(max_y), math.floor(max_z))


# Voxel definition describing whether a voxel contains certain elements
class Voxel(object):
    def __init__(self):
        self._building = None
        self._tree = None
        self._building_coverage = 0
        self._terrain = False

    # Assign building instance to voxel (when building coverage is greater then perevious, coverage not fully implemented)
    def set_building(self, building, coverage):
        if coverage > self._building_coverage:
            self._building = building
            self._building_coverage = coverage

    # Checks whether the voxel contains building
    def has_building(self):
        return self._building is not None

    # Returns building instance associated with the voxel
    def get_building(self):
        return self._building

    # Assign tree instance to voxel, only when the voxel does not contain building
    def set_tree(self, tree):
        if self._building == None:
            self._tree = tree

    # Return tree instance associated with the voxel
    def get_tree(self):
        return self._tree

    # Marks voxel as terrain voxel
    def set_terrain(self):
        self._terrain = True

    # Return terrain
    def get_terrain(self):
        return self._terrain

    # Checks if voxel is a terrain voxel
    def has_terrain(self):
        return self._terrain


# Voxel grid definition managing the individual voxel instances in a 3D grid
class VoxelGrid(object):
    def __init__(self):
        self._voxel_dict = dict()
        self._max_z = -999
        self._building_coord_list = None
        self._building_nr_list = None

    # Return the voxel associated with the location (via dictionary). When the voxel is not yet created, a new voxel
    # is initiated, added to the dictionary with the grid location as key, and maximum height of the model is adjusted
    # accordingly when necessary.
    def get_voxel(self, voxel_coord):
        if voxel_coord not in self._voxel_dict.keys():
            self._voxel_dict[voxel_coord] = Voxel()
            self._max_z = max(self._max_z, voxel_coord[2])
        v = self._voxel_dict[voxel_coord]
        assert isinstance(v, Voxel)
        return v

    # Return the z value of the highest point of the model
    def get_max_z(self):
        if self._max_z == -999:
            raise RuntimeError("max_z not defined for empty voxel grid")
        return self._max_z

    # Return a dictionary with voxels that contain elements (key is grid location, value is voxel instance)
    def get_voxel_dict(self):
        return self._voxel_dict

    # Returns a list of coordinates of voxels that contain building and a matching list with the corresponding building
    # numbers. When the lists do not exist yet, they are created.
    def get_building_list(self):
        if self._building_coord_list is None:
            self._building_coord_list = []
            self._building_nr_list = []
            for coord, voxel in self._voxel_dict.items():
                building = voxel.get_building()
                if building is not None:
                    self._building_coord_list.append(coord)
                    self._building_nr_list.append(building.get_nr())
        return self._building_coord_list, self._building_nr_list

    # Return a list with tree coordinate and tree instance couples
    def get_tree_list(self):
        tree_list = []
        for coord, voxel in self._voxel_dict.items():
            tree = voxel.get_tree()
            if tree:
                tree_list.append((coord, tree))
        return tree_list

    # Return 2.5D representation of buildings based on output model dimensions and generated voxel grid
    # (necessary for ENVI-met)
    def get_2d_building_grids(self, gridsxy, resolution):
        ztop_grid = []
        zbottom_grid = []
        zbuilding_nr_grid = []
        zfixed_height_grid = []
        for y in range(gridsxy[1]):
            ztop_grid.append(gridsxy[0] * [0])
            zbottom_grid.append(gridsxy[0] * [999])
            zbuilding_nr_grid.append(gridsxy[0] * [0])
            zfixed_height_grid.append(gridsxy[0] * [1])

        for coord, voxel in self._voxel_dict.items():
            x, y, z = coord
            if voxel.has_building():
                ztop_grid[y][x] = max(ztop_grid[y][x], (z + 1) * resolution[2])
                zbottom_grid[y][x] = min(zbottom_grid[y][x], z * resolution[2])
                zbuilding_nr_grid[y][x] = voxel.get_building().get_nr()

        for y in range(gridsxy[1]):
            for x in range(gridsxy[0]):
                if zbottom_grid[y][x] == 999:
                    zbottom_grid[y][x] = 0

        return ztop_grid, zbottom_grid, zbuilding_nr_grid, zfixed_height_grid

    # Return a list with the coordinates of the voxels that are marked as terrain
    def get_terrain_list(self):
        terrain_coord_list = []
        for coord, voxel in self._voxel_dict.items():
            if voxel.get_terrain():
                terrain_coord_list.append(coord)
        return terrain_coord_list

    # Return 2.5D version of terrain model
    def get_2d_terrain_grid(self, gridsxy, resolution):
        grid = []
        for y in range(gridsxy[1]):
            grid.append(gridsxy[0] * [0])

        for coord, voxel in self._voxel_dict.items():
            x, y, z = coord
            if voxel.has_terrain():
                grid[y][x] = max(grid[y][x], (z + 1) * resolution[2])

        return grid

    # Return dictionary with walls. Generated by checking for each building voxel, if there's a non-building voxel at
    # any of the six sides (these would need a wall in between), and if so, binding the wall to the correct voxel.
    def get_wall_dic(self):
        building_list = self.get_building_list()[0]
        wall_dic = {}
        for building_coord in building_list:

            # left (-x) side of voxel
            neighbour_coord = (building_coord[0] - 1, building_coord[1], building_coord[2])
            if neighbour_coord not in building_list:
                if building_coord not in wall_dic:
                    wall_dic[building_coord] = ["", "", ""]
                wall_dic[building_coord][0] = '000000'

            # front (-y) side of voxel
            neighbour_coord = (building_coord[0], building_coord[1] - 1, building_coord[2])
            if neighbour_coord not in building_list:
                if building_coord not in wall_dic:
                    wall_dic[building_coord] = ["", "", ""]
                wall_dic[building_coord][1] = '000000'

            # bottom (-z) side of voxel
            neighbour_coord = (building_coord[0], building_coord[1], building_coord[2] - 1)
            if neighbour_coord not in building_list:
                if building_coord not in wall_dic:
                    wall_dic[building_coord] = ["", "", ""]
                wall_dic[building_coord][2] = '000000'

            # right (+x) side of voxel
            neighbour_coord = (building_coord[0] + 1, building_coord[1], building_coord[2])
            if neighbour_coord not in building_list:
                if neighbour_coord not in wall_dic:
                    wall_dic[neighbour_coord] = ["", "", ""]
                wall_dic[neighbour_coord][0] = '000000'

            # back (+y) side of voxel
            neighbour_coord = (building_coord[0], building_coord[1] + 1, building_coord[2])
            if neighbour_coord not in building_list:
                if neighbour_coord not in wall_dic:
                    wall_dic[neighbour_coord] = ["", "", ""]
                wall_dic[neighbour_coord][1] = '000000'

            # top (+z) side of voxel
            neighbour_coord = (building_coord[0], building_coord[1], building_coord[2] + 1)
            if neighbour_coord not in building_list:
                if neighbour_coord not in wall_dic:
                    wall_dic[neighbour_coord] = ["", "", ""]
                wall_dic[neighbour_coord][2] = '000000'

        return wall_dic
