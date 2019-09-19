# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 17:42:23 2019
A small example of how to use the code of roverweb
@author: nixdorf
"""

from roverweb import *
import roverweb as rw
import geopandas as gpd
# load our dataset as geodataframe
gdfRover = rw.geometry.from_csv(
    './testdata/Mueglitz-20190708_'
    'selection.csv',
    geomtrcol=['LongDec', 'LatDec'],
    src_crs={
        'init': 'epsg:4326'
    })
# create a circular footprint (polygon) around the gps points
gdfRover_polyg = rw.geometry.point_to_circle(
    gdfRover, crcl_radius=100, number_of_points=10)

# Next is that we cluster our datasets in order to improve query speed
clusters = 3  # number of clusters
gdfRover_clust = rw.geometry.clustering(
    gdfRover_polyg, clusters=clusters, cluster_fld_name='ClusterID')
# %%having our clustered dataset we can run our osm to get desired secondary
# datasets
# group geodataframe
gdfRover_grouped = gdfRover_clust.groupby('ClusterID')
# create an empty result geodataframe
gdfRover_osm = gpd.GeoDataFrame()
# inititate the loop
for cluster_no in range(0, clusters):
    # create a subset of the original geodataframe
    gdfRover_subset = gdfRover_grouped.get_group(cluster_no)
    # create the queryboundstr
    querygeobound = rw.osm.querybound_generation(gdfRover_subset)
    # first query without special conditions
    gdfRover_subset = rw.osm.apnd_from_overpass(
        gdfRover_subset,
        querygeobound,
        queryfeatures={
            'way': ['landuse', 'water']
        })
    # another query for highway where values type service should be ignored
    gdfRover_subset = rw.osm.apnd_from_overpass(
        gdfRover_subset,
        querygeobound,
        queryfeatures={'way': ['highway']},
        value_out='service')
    # query only trees in landuse and count the no of feature
    gdfRover_subset = rw.osm.apnd_from_overpass(
        gdfRover_subset,
        querygeobound,
        queryfeatures={'node': ['natural']},
        value_in='tree',
        CountValues=True)
    # append subset back to entire dataset
    gdfRover_osm = gdfRover_osm.append(gdfRover_subset)

# %% Next we add data from soilgrid_network from wcs (faster than restapi)
gdfRover_osm_sg = rw.soilgrid.apnd_from_wcs(
    gdfRover_osm,
    soilgridlrs={
        'sg250m:ORCDRC_M_sl1_250m': 'ORCDRC_sl1',
        'sg250m:BLDFIE_M_sl1_250m': 'BLDFIE_sl1',
        'sg250m:CLYPPT_M_sl1_250m': 'CLYPPT_sl1'
    },
    raster_res=(250, 250),
    statistical_metric=['mean'],
    all_touched=True,
    output=None)

# Finally we get humidity values from dwd data, however current version of
# dwdweather library seems to have problems :-(
try:
    gdfRover_osm_sg_2 = rw.weather.apnd_dwd_stationdata(
        gdfRover_osm_sg,
        time_col='Date Time(UTC)',
        time_format='%Y-%m-%d %H:%M:%S',
        parameters=['airtemp_humidity'],
        temp_resolution='hourly')
except Exception as e:
    print(e)
    print('No DWD Station data retrieved')

#finally we write out our results to an shapefile
gdfRover_osm_sg.to_csv('./testdata/Mueglitz-20190708_'
                       'selection_roverweb.csv')
