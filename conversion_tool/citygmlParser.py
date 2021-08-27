import xml.etree.ElementTree as ET


class CitygmlParser(object):
    def __init__(self, input_citygml):
        self.filename = input_citygml
        self.tree = ET.parse(input_citygml)
        self.root = self.tree.getroot()

        # todo: check if citygml file and correct version (+ check if srs for objects same as envelope) + dimension is 3

    # get the the lower and upper coordinate of the 3D bounding box
    def get_bbox(self):
        for envelope in self.root.iter("{http://www.opengis.net/gml}Envelope"):
            l_corner = envelope.find('{http://www.opengis.net/gml}lowerCorner').text.split(' ')
            u_corner = envelope.find('{http://www.opengis.net/gml}upperCorner').text.split(' ')
            srs = envelope.attrib.get('srsName')
            return tuple(float(i) for i in l_corner), tuple(float(i) for i in u_corner), srs

    # extracts buildings as a list of surfaces, where surfaces are a list of coordinates (as tuples)
    # (from lod2Multisurface)
    def get_buildings(self):
        buildings = []
        for building in self.root.iter("{http://www.opengis.net/citygml/building/2.0}Building"):
            building_surfaces = []
            for multisurface in building.iter("{http://www.opengis.net/citygml/building/2.0}lod2MultiSurface"):
                for surface in multisurface.iter("{http://www.opengis.net/gml}LinearRing"):
                    pos_list = surface.find('{http://www.opengis.net/gml}posList').text.split(' ')
                    i = 0
                    surface_coords = []
                    while i < len(pos_list):
                        coord = (float(pos_list[i]), float(pos_list[i+1]), float(pos_list[i+2]))
                        surface_coords.append(coord)
                        i += 3
                    building_surfaces.append(surface_coords)
            buildings.append(building_surfaces)
        return buildings

    # extracts the terrain model as a list of surfaces, where surfaces are a list of coordinates (as tuples)
    # (from tin)
    def get_tin(self):
        tin_list = []
        for tin in self.root.iter("{http://www.opengis.net/citygml/relief/2.0}tin"):
            for triangle in tin.iter("{http://www.opengis.net/gml}LinearRing"):
                pos_list = triangle.find('{http://www.opengis.net/gml}posList').text.split(' ')
                i = 0
                triangle_coords = []
                while i < len(pos_list):
                    coord = (float(pos_list[i]), float(pos_list[i + 1]), float(pos_list[i + 2]))
                    triangle_coords.append(coord)
                    i += 3
                tin_list.append(triangle_coords)
        return tin_list

    # extracts the trees as a list of trees, represented as a tuple with a root point and a height
    # (from SolitaryVegetationObject - lod1ImplicitRepresentation)
    def get_trees(self):
        tree_list = []
        for tree in self.root.iter("{http://www.opengis.net/citygml/vegetation/2.0}SolitaryVegetationObject"):
            height = tree.find("{http://www.opengis.net/citygml/vegetation/2.0}height").text
            rep = tree.find("{http://www.opengis.net/citygml/vegetation/2.0}lod1ImplicitRepresentation")
            for ref in rep.iter("{http://www.opengis.net/citygml/2.0}referencePoint"):
                for point in ref.iter("{http://www.opengis.net/gml}pos"):
                    split_point = point.text.split(' ')
                    coord = tuple(float(i) for i in split_point)
                    tree_list.append((coord, float(height)))
        return tree_list
