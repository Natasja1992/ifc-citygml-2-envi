import ifcopenshell
import numpy
import ifcopenshell.geom as openshell_geom

import OCC.Core.gp
import OCC.Core.Geom

import OCC.Core.Bnd
import OCC.Core.BRepBndLib

import OCC.Core.BRep
import OCC.Core.BRepPrimAPI
import OCC.Core.BRepAlgoAPI
import OCC.Core.BRepBuilderAPI

import OCC.Core.GProp
import OCC.Core.BRepGProp

import OCC.Core.TopoDS
import OCC.Core.TopExp
import OCC.Core.TopAbs
from OCC.Core.BRepTools import breptools_OuterWire

settings = openshell_geom.settings()
settings.set(settings.USE_PYTHON_OPENCASCADE, True)


class IfcParser(object):
    def __init__(self, input_ifc):
        self.filename = input_ifc
        self.fh = ifcopenshell.open(input_ifc)

        # check if unit is in millimetres
        assert len(self.fh.by_type('IfcProject')) == 1, "expected only one project entity"
        project = self.fh.by_type('IfcProject')[0]
        for unit in project.UnitsInContext.Units:
            if unit.UnitType == 'LENGTHUNIT':
                assert str(unit.Prefix + unit.Name) == 'MILLIMETRE', "expected SI unit millimetre"

    # finds model rotation from north in degrees
    def model_orientation(self):
        assert len(self.fh.by_type('IfcProject')) == 1, "expected only one project entity"
        project = self.fh.by_type('IfcProject')[0]
        true_north = project.RepresentationContexts[0].TrueNorth.DirectionRatios
        return -numpy.degrees(numpy.arctan2(true_north[0], true_north[1]))

    # extracts reference point
    def reference_point(self):
        # reference point in ifc site
        assert len(self.fh.by_type('IfcSite')) == 1, "expected only one site entity"
        site = self.fh.by_type('IfcSite')[0]
        lat = site.RefLatitude     # WGS84
        lon = site.RefLongitude   # WGS84
        elevation = site.RefElevation   # relative to sea level
        latitude = float(lat[0]) + float(lat[1]) / 60 + (float(lat[2]) + float(lat[3])*1e-6) / (60 * 60)
        longitude = float(lon[0]) + float(lon[1]) / 60 + (float(lon[2]) + float(lon[3])*1e-6) / (60 * 60)

        return longitude, latitude, elevation

    # bounding box based on IfcSite representation
    def bbox_site(self):
        assert len(self.fh.by_type('IfcSite')) == 1, "expected only one site entity"
        site_rep = self.fh.by_type('IfcSite')[0].Representation

        min_x = float('inf')
        max_x = float('-inf')
        min_y = float('inf')
        max_y = float('-inf')
        min_z = float('inf')
        max_z = float('-inf')

        for representation in site_rep.Representations:
            if representation.RepresentationIdentifier == "Body":
                for item in representation.Items:
                    for faceset in item.FbsmFaces:
                        for face in faceset.CfsFaces:
                            for bound in face.Bounds:
                                for point in bound.Bound.Polygon:
                                    assert len(point.Coordinates) == 3, "expected x, y and z for coordinate"
                                    min_x = min(min_x, point.Coordinates[0])
                                    max_x = max(max_x, point.Coordinates[0])
                                    min_y = min(min_y, point.Coordinates[1])
                                    max_y = max(max_y, point.Coordinates[1])
                                    min_z = min(min_z, point.Coordinates[2])
                                    max_z = max(max_z, point.Coordinates[2])

        return (min_x, min_y, min_z), (max_x, max_y, max_z)

    # find adjacent edges to form polygon (without holes)
    def edges2polygon(self, edges):
        for edge in edges:
            assert len(edge) == 2

        for i in range(0, len(edges)-1):
            edge_found = False
            for j in range(i+1, len(edges)):
                if edges[i][1] == edges[j][0]:
                    # edge found, swap edges
                    edges[i+1], edges[j] = edges[j], edges[i+1]
                    edge_found = True
                    break

                elif edges[i][1] == edges[j][1]:
                    # swap point in edge
                    edges[j] = [edges[j][1], edges[j][0]]
                    # edge found, swap edges
                    edges[i+1], edges[j] = edges[j], edges[i+1]
                    edge_found = True
                    break

            if not edge_found:
                pass
            assert edge_found  # verify connecting edge has been found

        assert edges[0][0] == edges[-1][1]  # verify last edge connects to first edge

        polygons = [e[0] for e in edges]
        return polygons

    # extracting floors and roofs of the buildings as surfaces
    def get_slabs(self):
        slabs = self.fh.by_type('IfcSlab')
        faces = []

        for slab in slabs:
            shape = ifcopenshell.geom.create_shape(settings, slab).geometry

            exp = OCC.Core.TopExp.TopExp_Explorer(shape, OCC.Core.TopAbs.TopAbs_FACE)
            while exp.More():
                face = OCC.Core.TopoDS.topods.Face(exp.Current())
                exp.Next()
                surf = OCC.Core.BRep.BRep_Tool.Surface(face)
                plane = OCC.Core.Geom.Geom_Plane.DownCast(surf)
                if surf.DynamicType().Name() == "Geom_Plane":

                    if plane.Axis().Direction().Z() != 0:
                        faces.append(face)

        polygons = []
        for face in faces:
            wire = breptools_OuterWire(face)

            edge_exp = OCC.Core.TopExp.TopExp_Explorer(wire, OCC.Core.TopAbs.TopAbs_EDGE)
            edges = []
            while edge_exp.More():
                edge = OCC.Core.TopoDS.topods.Edge(edge_exp.Current())
                edge_exp.Next()
                vertex_exp = OCC.Core.TopExp.TopExp_Explorer(edge, OCC.Core.TopAbs.TopAbs_VERTEX)
                edge = []
                while vertex_exp.More():
                    vertex = OCC.Core.TopoDS.topods.Vertex(vertex_exp.Current())
                    vertex_exp.Next()
                    coord = OCC.Core.BRep.BRep_Tool.Pnt(vertex).Coord()
                    edge.append(coord)
                edges.append(edge)
            polygon = self.edges2polygon(edges)
            polygons.append(polygon)

        slabs_mm = []
        for polygon in polygons:
            polygon_mm = []
            for point in polygon:
                polygon_mm.append(tuple(1000.0*x for x in point))
            slabs_mm.append(polygon_mm)

        return slabs_mm

########################################################################################################################
#   One of the more promising approaches for implementing a full 3D building extraction, partly based on code          #
#   by Thomas Krijnen (http://academy.ifcopenshell.org/using-ifcopenshell-and-pythonocc-to-construct-new-geometry/).   #
#   Not implemented in the current version of the conversion tool, but the code might be of use for someone trying to  #
#   solve a similar problem. It should work for simpler shaped buildings. (it might have broken, in the process of     #
#   trying to make it better, but the idea is there).
########################################################################################################################

    # extracting (multiple) point(s) that are inside the building, by taking a point just underneath the center of
    # each roof element
    def get_points_inside(self):
        slabs = self.fh.by_type('IfcSlab')
        points = []
        for slab in slabs:
            if slab.PredefinedType == 'ROOF':
                # get point below center of roof
                shape = ifcopenshell.geom.create_shape(settings, slab).geometry
                bbox = OCC.Core.Bnd.Bnd_Box()
                OCC.Core.BRepBndLib.brepbndlib_Add(shape, bbox)
                center = ifcopenshell.geom.utils.get_bounding_box_center(bbox)
                center.SetZ(center.Z()-1.0)
                points.append(center)

        return points

    # extracting many points that are inside the building by taking the center point of inside walls
    def get_points_inside2(self):
        walls = self.fh.by_type('IfcWall')
        points = []
        for wall in walls:
            for definition in wall.IsDefinedBy:
                if definition.is_a('IfcRelDefinesByProperties'):
                    property_set = definition.RelatingPropertyDefinition
                    for prop in property_set.HasProperties:
                        if prop.is_a('IfcPropertySingleValue'):
                            if prop.Name != 'IsExternal':
                                if prop.NominalValue.wrappedValue:
                                    # get point below center of roof
                                    shape = ifcopenshell.geom.create_shape(settings, wall).geometry
                                    bbox = OCC.Core.Bnd.Bnd_Box()
                                    OCC.Core.BRepBndLib.brepbndlib_Add(shape, bbox)
                                    center = ifcopenshell.geom.utils.get_bounding_box_center(bbox)
                                    points.append(center)

        return points

    # Extracts the building as a watertight surface model. For simple shapes the center of the bounding box is used
    # as inside point to base the model around. For more complex shaped building multiple points can be used to create
    # multiple building parts. Still doesn't work perfectly.
    def get_building(self, point_inside=None):
        walls = self.fh.by_type('IfcWall')
        wall_shapes = []
        bbox = OCC.Core.Bnd.Bnd_Box()

        for wall in walls:
            for definition in wall.IsDefinedBy:
                if definition.is_a('IfcRelDefinesByProperties'):
                    property_set = definition.RelatingPropertyDefinition
                    for prop in property_set.HasProperties:
                        if prop.is_a('IfcPropertySingleValue'):
                            if prop.Name == 'IsExternal':
                                if prop.NominalValue.wrappedValue:
                                    shape = ifcopenshell.geom.create_shape(settings, wall).geometry

                                    wall_shapes.append((wall, shape))
                                    OCC.Core.BRepBndLib.brepbndlib_Add(shape, bbox)

        if not point_inside:
            point_inside = ifcopenshell.geom.utils.get_bounding_box_center(bbox)

        halfspaces = []
        for wall, shape in wall_shapes:
            exp = OCC.Core.TopExp.TopExp_Explorer(shape, OCC.Core.TopAbs.TopAbs_FACE)
            while exp.More():
                face = OCC.Core.TopoDS.topods.Face(exp.Current())
                exp.Next()
                surf = OCC.Core.BRep.BRep_Tool.Surface(face)
                assert surf.DynamicType().Name() == "Geom_Plane"

                plane = OCC.Core.Geom.Geom_Plane.DownCast(surf)

                if plane.Axis().Direction().Z() == 0:
                    face_bbox = OCC.Core.Bnd.Bnd_Box()
                    OCC.Core.BRepBndLib.brepbndlib_Add(face, face_bbox)
                    face_center = ifcopenshell.geom.utils.get_bounding_box_center(face_bbox).XYZ()

                    face_normal = plane.Axis().Direction().XYZ()
                    face_towards_center = point_inside.XYZ() - face_center
                    face_towards_center.Normalize()

                    dot = face_towards_center.Dot(face_normal)

                    if dot < -0.8:
                        face_plane = plane.Pln()
                        new_face = OCC.Core.BRepBuilderAPI.BRepBuilderAPI_MakeFace(face_plane).Face()
                        halfspace = OCC.Core.BRepPrimAPI.BRepPrimAPI_MakeHalfSpace(new_face,
                                                                                   point_inside).Solid()
                        halfspaces.append(halfspace)

        slabs = self.fh.by_type("IfcSlab")
        for slab in slabs:
            for definition in slab.IsDefinedBy:
                if definition.is_a('IfcRelDefinesByProperties'):
                    property_set = definition.RelatingPropertyDefinition
                    for prop in property_set.HasProperties:
                        if prop.is_a('IfcPropertySingleValue'):
                            if prop.Name == 'IsExternal':
                                if prop.NominalValue.wrappedValue:
                                    try:
                                        shape = ifcopenshell.geom.create_shape(settings, slab).geometry

                                        exp = OCC.Core.TopExp.TopExp_Explorer(shape, OCC.Core.TopAbs.TopAbs_FACE)
                                        while exp.More():
                                            face = OCC.Core.TopoDS.topods.Face(exp.Current())
                                            exp.Next()
                                            surf = OCC.Core.BRep.BRep_Tool.Surface(face)
                                            plane = OCC.Core.Geom.Geom_Plane.DownCast(surf)

                                            assert surf.DynamicType().Name() == "Geom_Plane"
                                            if plane.Axis().Direction().Z() == 0:
                                                face_plane = plane.Pln()
                                                new_face = OCC.Core.BRepBuilderAPI.BRepBuilderAPI_MakeFace(
                                                    face_plane).Face()
                                                halfspace = OCC.Core.BRepPrimAPI.BRepPrimAPI_MakeHalfSpace(new_face,
                                                                                                           point_inside).Solid()
                                                halfspaces.append(halfspace)
                                    except:
                                        print('skipping slab...')

        # Create an initial box from which to cut the halfspaces
        common_shape = OCC.Core.BRepPrimAPI.BRepPrimAPI_MakeBox(bbox.CornerMin(), bbox.CornerMax()).Solid()
        for halfspace in halfspaces:
            common_shape = OCC.Core.BRepAlgoAPI.BRepAlgoAPI_Common(common_shape, halfspace).Shape()

        faces = []
        exp = OCC.Core.TopExp.TopExp_Explorer(common_shape, OCC.Core.TopAbs.TopAbs_FACE)
        while exp.More():
            face = OCC.Core.TopoDS.topods.Face(exp.Current())
            exp.Next()
            faces.append(face)

        polygons = []
        for face in faces:
            edge_exp = OCC.Core.TopExp.TopExp_Explorer(face, OCC.Core.TopAbs.TopAbs_EDGE)
            edges = []
            while edge_exp.More():
                edge = OCC.Core.TopoDS.topods.Edge(edge_exp.Current())
                edge_exp.Next()
                vertex_exp = OCC.Core.TopExp.TopExp_Explorer(edge, OCC.Core.TopAbs.TopAbs_VERTEX)
                edge = []
                while vertex_exp.More():
                    vertex = OCC.Core.TopoDS.topods.Vertex(vertex_exp.Current())
                    vertex_exp.Next()
                    coord = OCC.Core.BRep.BRep_Tool.Pnt(vertex).Coord()
                    edge.append(coord)
                edges.append(edge)
            polygon = self.edges2polygon(edges)
            polygons.append(polygon)

        building_mm = []
        for polygon in polygons:
            polygon_mm = []
            for point in polygon:
                polygon_mm.append(tuple(1000.0*x for x in point))
            building_mm.append(polygon_mm)

        return building_mm
