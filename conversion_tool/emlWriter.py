import lxml.etree as et
from datetime import datetime
from conversionTool import ConversionTool


class EmlWriter:
    def __init__(self, conversion_tool):
        assert isinstance(conversion_tool, ConversionTool)
        self._conversion_tool = conversion_tool
        self._voxel_grid = self._conversion_tool.get_voxel_grid()
        self._root = et.Element('ENVI-MET_Datafile')
        self._gridsxy = self._conversion_tool.get_grids()
        self._gridsz = (self._voxel_grid.get_max_z() + 1) * 2

    def get_3dtag(self, tag, default):
        tag.set("type", "sparematrix-3D")
        tag.set("dataI", str(self._gridsxy[0]))
        tag.set("dataJ", str(self._gridsxy[1]))
        tag.set("zlayers", str(self._gridsz))
        tag.set("defaultValue", default)
        return tag

    def get_2dtag(self, tag):
        tag.set("type", "matrix-data")
        tag.set("dataI", str(self._gridsxy[0]))
        tag.set("dataJ", str(self._gridsxy[1]))
        return tag

    def create_header(self):
        header = et.SubElement(self._root, 'Header')
        et.SubElement(header, 'filetype').text = 'INPX ENVI-met Area Input File'
        et.SubElement(header, 'version').text = str(440)

        timestamp = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        et.SubElement(header, 'revisiondate').text = timestamp

        et.SubElement(header, 'remark').text = "Generated via conversion tool"
        et.SubElement(header, 'checksum').text = str(0)
        et.SubElement(header, 'encryptionlevel').text = str(0)

    def create_basedata(self):
        basedata = et.SubElement(self._root, 'baseData')

        et.SubElement(basedata, 'modelDescription').text = self._conversion_tool.get_model_description()

        et.SubElement(basedata, 'modelAuthor').text = "unknown"
        et.SubElement(basedata, 'modelcopyright').text = "unknown"

    def create_modelgeometry(self):
        modelgeometry = et.SubElement(self._root, 'modelGeometry')

        et.SubElement(modelgeometry, 'grids-I').text = str(self._gridsxy[0])
        et.SubElement(modelgeometry, 'grids-J').text = str(self._gridsxy[1])
        et.SubElement(modelgeometry, 'grids-Z').text = str(self._gridsz)

        (res_x, res_y, res_z) = self._conversion_tool.variables["resolution"]
        et.SubElement(modelgeometry, 'dx').text = "{:.5f}".format(res_x)
        et.SubElement(modelgeometry, 'dy').text = "{:.5f}".format(res_y)
        et.SubElement(modelgeometry, 'dz-base').text = "{:.5f}".format(res_z)

        # these elements are hard coded defaults now, could be added to user input if it seems necessary
        et.SubElement(modelgeometry, 'useTelescoping_grid').text = "0"  # no telescoping grid for now
        et.SubElement(modelgeometry, 'useSplitting').text = "0"  # splitting off
        et.SubElement(modelgeometry, 'verticalStretch').text = "0.00000"
        et.SubElement(modelgeometry, 'startStretch').text = "0.00000"
        et.SubElement(modelgeometry, 'has3DModel').text = "1"
        et.SubElement(modelgeometry, 'isFull3DDesign').text = "1"

    # nesting area is disabled, the software creates empty grid cells (without building/vegetation) around model, with
    # the actual terrain
    def create_nestingarea(self):
        nestingarea = et.SubElement(self._root, 'nestingArea')

        et.SubElement(nestingarea, 'numberNestinggrids').text = "0"
        et.SubElement(nestingarea, 'soilProfileA').text = "000000"
        et.SubElement(nestingarea, 'soilProfileB').text = "000000"

    def create_locationdata(self):
        locationdata = et.SubElement(self._root, 'locationData')

        # model rotation based on IFC model
        model_rotation = self._conversion_tool.get_model_orientation()
        et.SubElement(locationdata, 'modelRotation').text = "{:.5f}".format(model_rotation)

        et.SubElement(locationdata, 'projectionSystem').text = "EPSG:{}".format(self._conversion_tool.variables["reference_system_EPSG"])

        reference_point = self._conversion_tool.get_reference_point()
        et.SubElement(locationdata, 'realworldLowerLeft_X').text = str(reference_point[0])
        et.SubElement(locationdata, 'realworldLowerLeft_Y').text = str(reference_point[1])

        locationname = self._conversion_tool.get_location_name()
        et.SubElement(locationdata, 'locationName').text = "{}, {}".format(locationname["city"], locationname["country"])

        lonlatel = self._conversion_tool.get_lonlatel()
        et.SubElement(locationdata, 'location_Longitude').text = str(lonlatel[0])
        et.SubElement(locationdata, 'location_Latitude').text = str(lonlatel[1])

        timezone = self._conversion_tool.get_timezone()
        et.SubElement(locationdata, 'locationTimeZone_Name').text = str(timezone)
        et.SubElement(locationdata, 'locationTimeZone_Longitude').text = ""

    def create_defaultsetting(self):
        default = et.SubElement(self._root, 'defaultSettings')

        et.SubElement(default, 'commonWallMaterial').text = "000000"
        et.SubElement(default, 'commonRoofMaterial').text = "000000"

    def create_buildings2D(self):
        buildings2d = et.SubElement(self._root, 'buildings2D')

        ztop_grid, zbottom_grid, zbuilding_nr_grid, fixed_height_grid = \
            self._voxel_grid.get_2d_building_grids(self._gridsxy, self._conversion_tool.variables["resolution"])

        ztop_el = et.SubElement(buildings2d, 'zTop')
        self.get_2dtag(ztop_el)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = [str(x) for x in ztop_grid[y]]
            string += ",".join(row) + "\n    "
        ztop_el.text = string

        zbottom_el = et.SubElement(buildings2d, 'zBottom')
        self.get_2dtag(zbottom_el)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = [str(x) for x in zbottom_grid[y]]
            string += ",".join(row) + "\n    "
        zbottom_el.text = string

        zbuilding_nr_el = et.SubElement(buildings2d, 'buildingNr')
        self.get_2dtag(zbuilding_nr_el)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = [str(x) for x in zbuilding_nr_grid[y]]
            string += ",".join(row) + "\n    "
        zbuilding_nr_el.text = string

        zfixed_height_el = et.SubElement(buildings2d, 'fixedheight')
        self.get_2dtag(zfixed_height_el)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = [str(x) for x in fixed_height_grid[y]]
            string += ",".join(row) + "\n    "
        zfixed_height_el.text = string

    def create_simpleplants2D(self):
        simple_plants = et.SubElement(self._root, 'simpleplants2D')
        plants = et.SubElement(simple_plants, 'ID_plants1D')

        self.get_2dtag(plants)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = ["" for x in range(self._gridsxy[0])]
            string += ",".join(row) + "\n    "
        plants.text = string

    def create_3Dplants(self):
        trees = self._voxel_grid.get_tree_list()
        for coord, tree in trees:
            plants_3d = et.SubElement(self._root, 'totallynotadigit'+'3Dplants')  # hack: tags starting with digit not supported
            et.SubElement(plants_3d, 'rootcell_i').text = str(coord[0]+1)  # envi-met tree index starts at 1
            et.SubElement(plants_3d, 'rootcell_j').text = str(coord[1]+1)
            et.SubElement(plants_3d, 'rootcell_k').text = str(coord[2]+1)

            if tree.get_height() < 10:
                tree_id = '01SMDS'
                name = ".Spherical, medium trunk, dense, small (5m)"
            elif tree.get_height() < 20:
                tree_id = '01SMDM'
                name = ".Spherical, medium trunk, dense, medium (15m)"
            else:
                tree_id = '01SMDL'
                name = ".Spherical, medium trunk, dense, large (25m)"

            et.SubElement(plants_3d, 'plantID').text = tree_id
            et.SubElement(plants_3d, 'name').text = name
            et.SubElement(plants_3d, 'observe').text = "0"

    def create_soils2D(self):
        soils = et.SubElement(self._root, 'soils2D')
        profile = et.SubElement(soils, 'ID_soilprofile')

        self.get_2dtag(profile)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = ["000000" for x in range(self._gridsxy[0])]
            string += ",".join(row) + "\n    "
        profile.text = string

    def create_dem(self):
        dem = et.SubElement(self._root, 'dem')
        dem_ref = et.SubElement(dem, 'DEMReference')
        dem_ref.text = "0.00000"

        terrain_height = et.SubElement(dem, 'terrainheight')
        self.get_2dtag(terrain_height)

        grid = self._voxel_grid.get_2d_terrain_grid(self._gridsxy, self._conversion_tool.variables["resolution"])
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = [str(x) for x in grid[y]]
            string += ",".join(row) + "\n    "
        terrain_height.text = string

    def create_sources2D(self):
        sources = et.SubElement(self._root, 'sources2D')
        id_sources = et.SubElement(sources, 'ID_sources')

        self.get_2dtag(id_sources)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = ["" for x in range(self._gridsxy[0])]
            string += ",".join(row) + "\n    "
        id_sources.text = string

    def create_additionaldata(self):
        additional_data = et.SubElement(self._root, 'additionalData')
        db_link_point = et.SubElement(additional_data, 'db_link_point')
        db_link_area = et.SubElement(additional_data, 'db_link_area')

        self.get_2dtag(db_link_point)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = ["" for x in range(self._gridsxy[0])]
            string += ",".join(row) + "\n    "
        db_link_point.text = string

        self.get_2dtag(db_link_area)
        string = "\n    "
        for y in range(self._gridsxy[1] - 1, -1, -1):
            row = ["" for x in range(self._gridsxy[0])]
            string += ",".join(row) + "\n    "
        db_link_area.text = string

    def create_modelgeometry3D(self):
        modelgeom3d = et.SubElement(self._root, 'modelGeometry3D')
        et.SubElement(modelgeom3d, 'grids3D-I').text = str(self._gridsxy[0])
        et.SubElement(modelgeom3d, 'grids3D-J').text = str(self._gridsxy[1])
        et.SubElement(modelgeom3d, 'grids3D-K').text = str(self._gridsz)

    def create_buildings3D(self):
        buildings3d = et.SubElement(self._root, 'buildings3D')
        string = "\n    "
        building_flag = et.SubElement(buildings3d, 'buildingFlagAndNr')

        self.get_3dtag(building_flag, str(0))

        coords, nrs = self._voxel_grid.get_building_list()
        for i in range(len(coords)):
            t = coords[i] + (1, nrs[i])
            t = (str(v) for v in t)
            string += ",".join(t) + "\n    "

        building_flag.text = string

    def create_dem3D(self):
        dem3d = et.SubElement(self._root, 'dem3D')
        string = "\n    "
        terrain_flag = et.SubElement(dem3d, 'terrainflag')

        self.get_3dtag(terrain_flag, str(0))

        coords = self._voxel_grid.get_terrain_list()
        for coord in coords:
            t = coord + (1,)
            t = (str(v) for v in t)
            string += ",".join(t) + "\n    "

        terrain_flag.text = string

    def create_wallDB(self):
        wall = et.SubElement(self._root, 'WallDB')
        string = "\n    "
        wall_id = et.SubElement(wall, 'ID_wallDB')

        self.get_3dtag(wall_id, "")

        wall_dic = self._voxel_grid.get_wall_dic()
        for key, value in wall_dic.items():
            key = (str(v) for v in key)
            string += ",".join(key) + ',' + ",".join(value) + "\n    "

        wall_id.text = string

    def create_singlewallDB(self):
        singlewall = et.SubElement(self._root, 'SingleWallDB')
        string = "\n    "
        singlewall_id = et.SubElement(singlewall, 'ID_singlewallDB')

        self.get_3dtag(singlewall_id, "")

        singlewall_id.text = string

    def create_greeningDB(self):
        greening = et.SubElement(self._root, 'GreeningDB')
        string = "\n    "
        greening_id = et.SubElement(greening, 'ID_greeningDB')

        self.get_3dtag(greening_id, "")

        greening_id.text = string

    # create file structure and write to file
    def write(self, output_file):
        print("creating 'ENVI-met area input file'...")

        self.create_header()
        self.create_basedata()
        self.create_modelgeometry()
        self.create_nestingarea()
        self.create_locationdata()
        self.create_defaultsetting()
        self.create_buildings2D()
        self.create_simpleplants2D()
        self.create_3Dplants()
        self.create_soils2D()
        self.create_dem()
        self.create_sources2D()
        self.create_additionaldata()
        self.create_modelgeometry3D()
        self.create_buildings3D()
        self.create_dem3D()
        self.create_wallDB()
        self.create_singlewallDB()
        self.create_greeningDB()

        print('writing to file [%s]...' % output_file)
        mydata = et.tostring(self._root, pretty_print=True)
        mydata = mydata.replace(b'totallynotadigit', b'')  #hack: tags starting with digit not supported
        with open(output_file, "wb") as fh:
            fh.write(mydata)

