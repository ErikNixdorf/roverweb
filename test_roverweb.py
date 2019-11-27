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
    './testdata/Mueglitz-20190708_selection.csv',
    geomtrcol=['LongDec', 'LatDec'],
    src_crs={
        'init': 'epsg:4326'
    })
# create a circular footprint (polygon) around the gps points
gdfRover_circls = rw.geometry.points_to_circle(
    gdfRover, crcl_radius=100, number_of_points=10)

# instead we can also create a true footprint which consideres that the rover is moving,
# basically a parallelogram attached on both sides with two half circles
gdfRover_polyg=rw.geometry.points_to_footprint(gdfRover,footprint_radius=50,crs_src = "epsg:4326",crs_dst = "epsg:25833",number_of_points=10,inpt_geo_type='Point')
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
        values_out=['service','footway'])
    # query only trees in landuse and count the no of feature
    gdfRover_subset = rw.osm.apnd_from_overpass(
        gdfRover_subset,
        querygeobound,
        queryfeatures={'node': ['natural']},
        values_in=['tree'],
        CountValues=True)
    # append subset back to entire dataset
    gdfRover_osm = gdfRover_osm.append(gdfRover_subset,sort=True)

# %% Next we add data from soilgrid_network from wcs (faster than restapi)
gdfRover_osm_sg = rw.soilgrid.apnd_from_wcs(
    gdfRover_osm.copy(),
    soilgridlrs={
        'sg250m:ORCDRC_M_sl1_250m': 'ORCDRC_sl1',
        'sg250m:BLDFIE_M_sl1_250m': 'BLDFIE_sl1',
        'sg250m:CLYPPT_M_sl1_250m': 'CLYPPT_sl1'
    },
    raster_res=(250, 250),
    statistical_metric=['mean'],
    all_touched=True,
    output=None)

# Finally we get humidity values from dwd data, 
gdfRover_osm_sg_2 = rw.weather.apnd_dwd_stationdata(
    gdfRover_osm_sg.copy(),
    time_col='Date Time(UTC)',
    time_format='%Y-%m-%d %H:%M:%S',
    data_category=['air_temperature'],
    parameters=['airtemp_humidity','airtemp_temperature'],
    temp_resolution='hourly',
    no_of_nearest_stations=3)    

#finally we write out our results to an shapefile
gdfRover_osm_sg_2.to_csv('./output/Mueglitz-20190708_'
                       'selection_roverweb.csv')

#change datetime to str
gdfRover_osm_sg_2['Date Time(UTC)']=gdfRover_osm_sg_2['Date Time(UTC)'].dt.strftime('%Y-%m-%d %H:%M:%S%')
#save as shp
gdfRover_osm_sg_2.to_file('./output/Mueglitz-20190708_'
                       'selection_roverweb.shp')
