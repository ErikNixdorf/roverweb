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
from shapely.geometry import LineString
from shapely.geometry import Point
from scipy.cluster.vq import kmeans, vq


def from_csv(inpt_data, geomtrcol=['LongDec', 'LatDec'],
             src_crs={'init': 'epsg:4326'}):
    if inpt_data.endswith('.csv'):
        df = pd.read_csv(inpt_data)
        # delete all rows where Lat or Lon is Na
        df.dropna(subset=geomtrcol, inplace=True)
        #  create a geodataframe
        inpt_data = gpd.GeoDataFrame(df,
                                     geometry=[Point(xy)for xy in
                                               zip(df[geomtrcol[0]],
                                                   df[geomtrcol[1]])])
        inpt_data.crs = src_crs
        # check CRS and convert to EPSG if required
        src_epsg = inpt_data.crs.get('init')
        if src_epsg != 'epsg:4326':
            inpt_data['geometry'] = inpt_data['geometry'].to_crs(epsg=4326)
            print('geometry converted to EPSG:4326 for processing')
        return inpt_data
    else:
        return None


def clustering(inpt_data, geomtrcol=['LongDec', 'LatDec'],
               src_crs={'init': 'epsg:4326'}, clusters=3,
               cluster_fld_name='ClusterID'):
    """
    This function takes either a csv input file or shape or geodataframe
    as input and starts the clustering,
    result is geodataframe with point geometry
    """
    if isinstance(inpt_data, str):
        if inpt_data.endswith('.csv'):
            inpt_data = from_csv(inpt_data(geomtrcol=['LongDec', 'LatDec'],
                                           src_crs={'init': 'epsg:4326'}))
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


def PointsInCircum(center, radius, number_of_pnts):
    """
    Simple function to create a number of points on a circle in equal dist
    in cartesian coordinates
    """
    circle_pnts = [(center[0] + math.cos(2*math.pi/number_of_pnts*x)*radius,
                    center[1] + math.sin(2*math.pi/number_of_pnts*x)*radius)
                   for x in range(0, number_of_pnts+1)]

    return np.asarray(circle_pnts)


def point_to_circle(gdfin, crcl_radius=100,
                               number_of_points=10):
    """
    Create circles which cover the relevant geometry for geographic coordinat
    """
    angles = np.linspace(0, 360, number_of_points)
    circle_cord = [geog.propagate(xy, angles, crcl_radius)
                   for xy in zip(gdfin.geometry.x, gdfin.geometry.y)]
    circle_geom = [shapely.geometry.Polygon(zip(circle_cord[i][:, 0],
                                                circle_cord[i][:, 1]))
                   for i in range(0, len(circle_cord))]
    gdfout = gdfin[:]
    gdfout.geometry = circle_geom
    return gdfout

def closest_point_of_neighbor_ring(polyg0, polyg1):
    """
    Function calculates the point on ring polyg0 which is closest
    to polyg 0 and provides back its ID in polyg0
    """
    polyg0_start = nearest_points(polyg1.centroid,
                                polyg0)[1]
    #extract all points from polyg_0
    polyg0_vertices = np.array(polyg0.exterior.coords)
    #find id of nearest pnt
    id_polyg0_closest_point = np.where(np.logical_and(polyg0_vertices[:,0]==
                                            polyg0_start.coords[0][0],
                                            polyg0_vertices[:,1]==
                                            polyg0_start.coords[0][1]))
    #output as integer
    return polyg0_start,int(id_polyg0_closest_point[0])

def extrap_line_creation(point0,point1,max_extrap=100):
    """
    A function which creates an extrapolated shapely Linestring
    from two points. Line connects the points and extends line in both
    directions from points by a length defined by "max_extrap"
    Output is a shapely linestring
    """
    #determine slope
    m=(point0.coords[0][1]-point1.coords[0][1])/(
            point0.coords[0][0]-point1.coords[0][0])
    #Use basic algebra and geometry to create a line that
    #will extend in both directions as defined
    r=np.sqrt(1+m**2)
    polyl_extrap=LineString([ (point0.coords[0][0]+max_extrap/r,
                           point0.coords[0][1]+max_extrap*m/r),
                          (point0.coords[0][0]-max_extrap/r,
                           point0.coords[0][1]-max_extrap*m/r),
                           (point1.coords[0][0]+max_extrap/r,
                           point1.coords[0][1]+max_extrap*m/r),
                          (point1.coords[0][0]-max_extrap/r,
                           point1.coords[0][1]-max_extrap*m/r)])
    return polyl_extrap
