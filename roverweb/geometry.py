# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:16:22 2019
This tool creates correct poylgons from point measurements considering the footprint of the cosmic ray rover measurements
-footprint of stationary measurement is a circle
-as the device is moving the position is measured after the 10s measurement period
-->footprint is either a) geometric object looking like a cilindy with two half circle (=)
                        b) an ellipse if the circle at start time and end time intersect
@author: nixdorf
"""
import geog
import geopandas as gpd
import math
import numpy as np
import pandas as pd
import shapely
from shapely.ops import nearest_points
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from scipy.cluster.vq import kmeans, vq


def from_csv(inpt_data,
             geomtrcol=['LongDec', 'LatDec'],
             src_crs={
                 'init': 'epsg:4326'
             }):
    if inpt_data.endswith('.csv'):
        df = pd.read_csv(inpt_data)
        # delete all rows where Lat or Lon is Na
        df.dropna(subset=geomtrcol, inplace=True)
        # repair the index
        df.index = range(0, len(df))
        #  create a geodataframe
        inpt_data = gpd.GeoDataFrame(
            df,
            geometry=[
                Point(xy) for xy in zip(df[geomtrcol[0]], df[geomtrcol[1]])
            ])
        inpt_data.crs = src_crs
        # check CRS and convert to EPSG if required
        src_epsg = inpt_data.crs
        if src_epsg != 'epsg:4326':
            inpt_data['geometry'] = inpt_data['geometry'].to_crs(epsg=4326)
            print('geometry converted to EPSG:4326 for processing')
        return inpt_data
    else:
        return None


def clustering(inpt_data,
               geomtrcol=['LongDec', 'LatDec'],
               src_crs={'init': 'epsg:4326'},
               clusters=3,
               cluster_fld_name='ClusterID'):
    """
    This function takes either a csv input file or shape or geodataframe
    as input and starts the clustering,
    result is geodataframe with point geometry
    """
    if isinstance(inpt_data, str):
        if inpt_data.endswith('.csv'):
            inpt_data = from_csv(
                inpt_data(
                    geomtrcol=['LongDec', 'LatDec'],
                    src_crs={
                        'init': 'epsg:4326'
                    }))
        else:
            # try to read directly as geodataframe
            try:
                inpt_data = gpd.read_file(inpt_data)
            except Exception as e:
                print(e)
                print('only vector GIS formats supported')

    # now we start the clustering depending on the datatype
    if inpt_data.iloc[0].geometry.type == 'Polygon':
        points_xy = np.column_stack((inpt_data.geometry.centroid.x,
                                     inpt_data.geometry.centroid.y))
    if inpt_data.iloc[0].geometry.type in ['Point', 'LineString']:
        points_xy = np.column_stack((inpt_data.geometry.x,
                                     inpt_data.geometry.y))
    centroids, _ = kmeans(points_xy, clusters)
    # assign each sample to a cluster
    idx, _ = vq(points_xy, centroids)
    # add Cluster_ID to dataframe
    outpt_data = inpt_data[:]  # create a copy
    outpt_data.insert(outpt_data.shape[1], cluster_fld_name, idx)
    print(clusters, ' cluster were created by k-mean method')
    return outpt_data


#%% Define necessary subfunction (maths
def PointsInCircum(center, radius, number_of_pnts):
    """
    Simple function to create a number of points on a circle in equal dist
    """
    circle_pnts = [
        (center[0] + math.cos(2 * math.pi / number_of_pnts * x) * radius,
         center[1] + math.sin(2 * math.pi / number_of_pnts * x) * radius)
        for x in range(0, number_of_pnts + 1)
    ]

    return np.asarray(circle_pnts)


def closest_point_of_neighbor_ring(polyg0, polyg1):
    """
    Function calculates the point on ring polyg0 which is closest
    to polyg 1 centroid and provides back its ID in polyg0
    """
    index = 0
    distances = np.ones((len(polyg0.exterior.coords)))
    for vertex in polyg0.exterior.coords:
        distances[index] = ((vertex[0] - polyg1.centroid.coords[0][0])**2 +
                            (vertex[1] - polyg1.centroid.coords[0][1])**2)**0.5
        index = index + 1
    # get min_argument
    id_polyg0_closest = np.argmin(distances)
    polyg0_closest = Point(polyg0.exterior.coords[int(id_polyg0_closest)])
    return polyg0_closest, id_polyg0_closest


def points_to_circle(gdfin, crcl_radius=100, number_of_points=10):
    """
    Create circles which cover the relevant geometry for geographic coordinat
    """
    angles = np.linspace(0, 360, number_of_points)
    circle_cord = [
        geog.propagate(xy, angles, crcl_radius)
        for xy in zip(gdfin.geometry.x, gdfin.geometry.y)
    ]
    circle_geom = [
        shapely.geometry.Polygon(
            zip(circle_cord[i][:, 0], circle_cord[i][:, 1]))
        for i in range(0, len(circle_cord))
    ]
    gdfout = gdfin[:]
    gdfout.geometry = circle_geom
    return gdfout


def extrap_line_creation(point0, point1, max_extrap=100):
    """
    A function which creates an extrapolated shapely Linestring
    from two points. Line connects the points and extends line in both
    directions from points by a length defined by "max_extrap"
    Output is a shapely linestring
    """
    #determine slope
    m = (point0.coords[0][1] - point1.coords[0][1]) / (
        point0.coords[0][0] - point1.coords[0][0])
    #Use basic algebra and geometry to create a line that
    #will extend in both directions as defined
    r = np.sqrt(1 + m**2)
    polyl_extrap = LineString([(point0.coords[0][0] + max_extrap / r,
                                point0.coords[0][1] + max_extrap * m / r),
                               (point0.coords[0][0] - max_extrap / r,
                                point0.coords[0][1] - max_extrap * m / r),
                               (point1.coords[0][0] + max_extrap / r,
                                point1.coords[0][1] + max_extrap * m / r),
                               (point1.coords[0][0] - max_extrap / r,
                                point1.coords[0][1] - max_extrap * m / r)])
    return polyl_extrap


def arc_vertices(polyg, ID_vertex_start, ID_vertex_end, add_bounds=False):
    """
    A function to find all vertices on a polyline/polygon from start
    to end vertex. Output as list of points. Start/End Vertex are
    excluded as long as "add_bounds" is false
    """
    if add_bounds:
        ID_increment = 0
    else:
        ID_increment = 1

    if ID_vertex_end > ID_vertex_start:
        arc_vertices = [
            Point(polyg.exterior.coords[i])
            for i in range(ID_vertex_start + ID_increment, ID_vertex_end)
        ]
    else:  # the other case is a bit more complicated
        ID_arc_vertices = np.hstack((range(ID_vertex_start + ID_increment,
                                           len(polyg.exterior.coords) - 1),
                                     range(0, ID_vertex_end)))
        arc_vertices = [
            Point(polyg.exterior.coords[int(i)]) for i in ID_arc_vertices
        ]
    return arc_vertices

def swap_xy(geom):
    if geom.is_empty:
        return geom

    if geom.has_z:
        def swap_xy_coords(coords):
            for x, y, z in coords:
                yield (y, x, z)
    else:
        def swap_xy_coords(coords):
            for x, y in coords:
                yield (y, x)

    # Process coordinates from each supported geometry type
    if geom.type in ('Point', 'LineString', 'LinearRing'):
        return type(geom)(list(swap_xy_coords(geom.coords)))
    elif geom.type == 'Polygon':
        ring = geom.exterior
        shell = type(ring)(list(swap_xy_coords(ring.coords)))
        holes = list(geom.interiors)
        for pos, ring in enumerate(holes):
            holes[pos] = type(ring)(list(swap_xy_coords(ring.coords)))
        return type(geom)(shell, holes)
    elif geom.type.startswith('Multi') or geom.type == 'GeometryCollection':
        # Recursive call
        return type(geom)([swap_xy(part) for part in geom.geoms])
    else:
        raise ValueError('Type %r not recognized' % geom.type)


#%% the function to map real_footprint to data
def points_to_footprint(inpt_data,
                        footprint_radius=50,
                        crs_src="epsg:4326",
                        crs_dst="epsg:25833",
                        number_of_points=10,
                        inpt_geo_type='Point'):
    """
    This tool creates correct poylgons from point measurements considering the
    footprint of the cosmic ray rover measurements
    footprint of stationary measurement is a circle as the device is moving
    the position is measured after the 10s measurement period
    -->footprint is either
    a) geometric object looking like a cilindy with two half circle (=)
    b) an ellipse if the circle at start time and end time intersect
    number_of_points : Number of points to create a circle for each measurement
                """

    #convert points to stationary footprint (circles) first
    if inpt_geo_type == 'Point':
        inpt_crlcs = points_to_circle(
            inpt_data,
            crcl_radius=footprint_radius,
            number_of_points=number_of_points)
    else:
        inpt_crlcs = inpt_data[:]
    # convert coordinate system if necessary
    inpt_crlcs['geometry'] = inpt_crlcs['geometry'].to_crs(crs_dst)
    #get the geometry
    crlc_plgns = inpt_crlcs.geometry

    #%% now we start our scheme to merge the circles to
    # polygons/Ellipsoids representing the footprint
    footprnt_geomrts = [None] * (len(crlc_plgns) - 1)

    for idpnt in range(1, len(crlc_plgns)):  #the first value is ignored
        #check whether circles are intersecting with preceding one
        try:
            polyg_pair = MultiPolygon([crlc_plgns[idpnt - 1], crlc_plgns[idpnt]])
        except:
            #print('test with iloc')
            polyg_pair = MultiPolygon([crlc_plgns.iloc[idpnt - 1], crlc_plgns.iloc[idpnt]])
        polyg_pair = gpd.GeoSeries(polyg_pair)
        polyg_hull = polyg_pair.convex_hull

        # add to the geometry
        footprnt_geomrts[idpnt - 1] = polyg_hull[0]

    print('finished footprint_creation')

    # Finally we replace all geometries of the original file with the new one and save
    outpt_data = inpt_data.copy(deep=True)
    #delete the first row
    outpt_data = outpt_data.iloc[1:]
    #add new geometry
    outpt_data['geometry'] = footprnt_geomrts
    # convert coordinate system back
    outpt_data.crs=(crs_dst)
    outpt_data['geometry'] = outpt_data['geometry'].to_crs(crs_src)
    #swapping helps as pyproj 2.1.3 seems to have some errors
    #outpt_data.geometry = outpt_data.geometry.map(swap_xy)
    outpt_data.crs=(crs_src)
    print('geotransformation back to original', crs_src)

    return outpt_data
