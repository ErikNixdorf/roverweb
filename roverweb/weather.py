# -*- coding: utf-8 -*-
"""
A small module which adds data from weather stations to our rovertracks,
currently it havily based on dwdweather2 library. In future also raster data
should be used

@author: nixdorf
"""
from dwdweather import DwdWeather
from datetime import datetime
import numpy as np
import geopandas as gpd


def apnd_dwd_stationdata(inpt_data, time_col='Date Time(UTC)',
                         time_format='%Y-%m-%d %H:%M:%S',
                         parameters=['airtemp_humidity'],
                         temp_resolution='hourly'):
    print('Start quering data from DWD')
    # Define outpt_data
    outpt_data = gpd.GeoDataFrame()
    # First we create an object from class DwDWeather
    dw = DwdWeather(resolution=temp_resolution)

    # now we have to loop through each object of gdf and find nearest stat
    for _, row in inpt_data.iterrows():

        # now we extract the centroid/point xy of each entry
        if row.geometry.type == 'Polygon':
            point_xy = np.array((row.geometry.centroid.x,
                                 row.geometry.centroid.y))
        if row.geometry.type in ['Point', 'LineString']:
            point_xy = np.array((row.geometry.x, row.geometry.y))
        # next we convert our date to a datetime object
        query_hour = datetime.strptime(row[time_col], time_format)
        # get closest station
        closest = dw.nearest_station(lon=point_xy[0], lat=point_xy[1])
        # query for result
        query_result = dw.query(station_id=closest["station_id"],
                                timestamp=query_hour)
        # extract the requested parameter from the entire data
        for parameter in parameters:
            row[parameter] = query_result[parameter]
        # append datasets
        outpt_data = outpt_data.append(row)
    print('finished quering data from DWD')
    return outpt_data
