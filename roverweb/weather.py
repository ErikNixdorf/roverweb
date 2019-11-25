# -*- coding: utf-8 -*-
"""
A small module which adds data from weather stations to our rovertracks,
currently it havily based on dwdweather2 library. In future also raster data
should be used
WARNING: CURRENTLY WORKS ONLY WITH DWDWEATHER 0.9 DUE TO Python 3 compat
issues
@author: nixdorf
"""
from dwdweather import DwdWeather
import numpy as np
import geopandas as gpd
from io import StringIO
import math
import pandas as pd

from shapely.ops import nearest_points


#%% some function we need to parse our datasets
def nearest_direct(row, gdf1, gdf2, src_column=None):
    """Find the nearest point and return the corresponding value from specified column.
    df 1 is the origin
    df 2 is the destination
    inspired by https://automating-gis-processes.github.io/2017/lessons/L3/nearest-neighbour.html
    """
    #create a unary union
    unary_union = gdf2.unary_union
    # Find the geometry that is closest
    nearest = gdf2['centroid'] == nearest_points(row['centroid'],
                                                 unary_union)[1]
    # Get the corresponding value from df2 (matching is based on the geometry)
    value = gdf2[nearest][src_column].get_values()[0]
    return value


def nearest_loop(row,
                 gdf2,
                 geometry_cols=['geo_lon', 'geo_lat'],
                 src_column=None):
    """
    takes longer, seems to be more precise
    """

    def haversine_distance(origin, destination):
        lon1, lat1 = origin
        lon2, lat2 = destination
        radius = 6371000  # meters

        dlat = math.radians(lat2 - lat1)
        dlon = math.radians(lon2 - lon1)
        a = math.sin(dlat / 2) * math.sin(dlat / 2) + math.cos(
            math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(
                dlon / 2) * math.sin(dlon / 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        d = radius * c
        return d

    # start the main iteration
    if row.geometry.type == 'Polygon':
        point_xy = np.array((row.geometry.centroid.x, row.geometry.centroid.y))
    if row.geometry.type in ['Point', 'LineString']:
        point_xy = np.array((row.geometry.x, row.geometry.y))
    # Select most current stations datasets.
    closest = None
    closest_distance = 99999999999
    for _, station in gdf2.iterrows():
        d = haversine_distance(
            (point_xy[0], point_xy[1]),
            (station[geometry_cols[0]], station[geometry_cols[1]]))
        if d < closest_distance:
            closest = station
            closest_distance = d
    return closest[src_column]


#%% The main function


def apnd_dwd_stationdata(inpt_data,
                         time_col='Date Time(UTC)',
                         time_format='%Y-%m-%d %H:%M:%S',
                         parameters=['airtemp_humidity'],
                         temp_resolution='hourly',
                         fixed_network=True):    
    print('Start quering data from DWD')
    # Define outpt_data
    outpt_data = gpd.GeoDataFrame()
    # next we convert our date to a datetime object
    inpt_data[time_col] = pd.to_datetime(
        inpt_data[time_col], format=time_format)
    # add a centroid column
    inpt_data['centroid'] = inpt_data.centroid
    # First we create an object from class DwDWeather
    dw = DwdWeather(resolution=temp_resolution, category_names=parameters)
    #as dwdweather 0.9 has some issues, we can get the station metadata by ourself
    if fixed_network:
        stations_all = pd.read_csv(
            '.\\roverweb\\dbase\\stations_air_temp_hourly.txt',
            usecols=[0, 1, 2, 3, 4, 5])
    else:
        #update the sql database
        dw.import_stations(
        )  # we check all available stations and create a valid list
        stations_all = pd.read_csv(StringIO(dw.stations_csv()))
    #convert to datetime
    stations_all.date_end = pd.to_datetime(
        stations_all.date_end, format='%Y%m%d')
    stations_all.date_start = pd.to_datetime(
        stations_all.date_start, format='%Y%m%d')
    # we clean to all stations which cover the campaign time
    #dt_low <= dt <= dt_high:
    dwd_stations = stations_all[
        (stations_all.date_start <= inpt_data.iloc[0][time_col])
        & (inpt_data.iloc[0][time_col] <= stations_all.date_end)]
    #make a geodataframe out of it
    dwd_stations = gpd.GeoDataFrame(
        dwd_stations,
        geometry=gpd.points_from_xy(dwd_stations.geo_lon,
                                    dwd_stations.geo_lat))
    #add a centroid
    dwd_stations['centroid'] = dwd_stations.centroid
    inpt_data['nearest_dwd_id'] = inpt_data.apply(
        nearest_direct,
        gdf1=inpt_data.copy(),
        gdf2=dwd_stations.copy(),
        src_column='station_id',
        axis=1)
    #get the same ones by using a looped approach

    #inpt_data['nearest_dwd_id'] = inpt_data.apply(
        #nearest_loop,
        #gdf2=dwd_stations,
        #geometry_cols=['geo_lon', 'geo_lat'],
        #src_column='station_id',
        #axis=1)

    #delete the centroid_column
    inpt_data.drop(columns='centroid',inplace=True)

    #%%loop trough all lines to get the parameterset
    for _, row in inpt_data.iterrows():
        # query for result
        query_result = dw.query(
            station_id=row['nearest_dwd_id'], timestamp=row[time_col])
        # extract the requested parameter from the entire data
        for parameter in parameters:
            row[parameter] = query_result[parameter]
    # append datasets
        outpt_data = outpt_data.append(row)
    print('finished quering data from DWD')
    return outpt_data
