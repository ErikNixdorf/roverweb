"""
A submodule which allows access to soilgrid WebmappingService
"""

# import librarys
import geog
from owslib.wcs import WebCoverageService
import numpy as np

from io import BytesIO
from rasterio.io import MemoryFile
import rasterstats
import geopandas as gpd
import requests
import time
def connect_wcs(server = 'https://data.isric.org/geoserver/ows',version='1.0.0',connected=False,trials=1):
    trial=0
    while not connected:
        if trial>=trials:
            print('Maximum connection tests of', trials, 'reached')
            return None
        try:
            wcs = WebCoverageService(server,
                             version=version)
            connected = True
    
        except:
            time.sleep(1)
            print('Reconnect to Server after 5 seconds')
            trial+=1
            pass
    return wcs


def apnd_from_wcs(inpt_data,
                          soilgridlrs={'sg250m:ORCDRC_M_sl1_250m': 'ORCDRC_sl1',
                                       'sg250m:BLDFIE_M_sl1_250m': 'BLDFIE_sl1',
                                       'sg250m:PHIHOX_M_sl1_250m': 'PHIHOX_sl1',
                                       'sg250m:SNDPPT_M_sl1_250m': 'SNDPPT_sl1',
                                       'sg250m:SLTPPT_M_sl1_250m': 'SLTPPT_sl1',
                                       'sg250m:CLYPPT_M_sl1_250m': 'CLYPPT_sl1',
                                       'sg250m:CRFVOL_M_sl1_250m': 'CRFVOL_sl1',
                                       'sg250m:CECSOL_M_sl1_250m': 'CECSOL_sl1',
                                       'sg250m:ORCDRC_M_sl2_250m': 'ORCDRC_sl2',
                                       'sg250m:BLDFIE_M_sl2_250m': 'BLDFIE_sl2',
                                       'sg250m:PHIHOX_M_sl2_250m': 'PHIHOX_sl2',
                                       'sg250m:SNDPPT_M_sl2_250m': 'SNDPPT_sl2',
                                       'sg250m:SLTPPT_M_sl2_250m': 'SLTPPT_sl2',
                                       'sg250m:CLYPPT_M_sl2_250m': 'CLYPPT_sl2',
                                       'sg250m:CRFVOL_M_sl2_250m': 'CRFVOL_sl2',
                                       'sg250m:CECSOL_M_sl2_250m': 'CECSOL_sl2',
                                       'sg250m:ORCDRC_M_sl3_250m': 'ORCDRC_sl3',
                                       'sg250m:BLDFIE_M_sl3_250m': 'BLDFIE_sl3',
                                       'sg250m:PHIHOX_M_sl3_250m': 'PHIHOX_sl3',
                                       'sg250m:SNDPPT_M_sl3_250m': 'SNDPPT_sl3',
                                       'sg250m:SLTPPT_M_sl3_250m': 'SLTPPT_sl3',
                                       'sg250m:CLYPPT_M_sl3_250m': 'CLYPPT_sl3',
                                       'sg250m:CRFVOL_M_sl3_250m': 'CRFVOL_sl3',
                                       'sg250m:CECSOL_M_sl3_250m': 'CECSOL_sl3',
                                       'sg250m:BDTICM_M_250m': 'BDTICM',
                                       },
                          raster_res=(250, 250), statistical_metric=['mean'],
                          all_touched=True, output=None):
    """
    A function to retrieve layers from soildgrid Wec Coverage Service and maps
    retrieved values on polygons
    raster_resolution= (longitudinal_res in m, latitudinal_res in m)
    gdf: Can be either a geodataframe or a shapefile, if shapefile it will be
    converted to
    #https://www.isric.org/explore/soilgrids/faq-soilgrids
    """
    # First, if dataset is shapefile, we read it to geodataframe
    if isinstance(inpt_data, str):
        inpt_data = gpd.read_file(inpt_data)
    # check source coordinate system and convert to EPSG 4326
    src_epsg = inpt_data.crs
    if src_epsg.srs != 'epsg:4326':
        inpt_data['geometry'] = inpt_data['geometry'].to_crs(epsg=4326)
        print('geometry converted to EPSG:4326 for processing')
    # contact the wcs
    wcs = connect_wcs(server = 'https://data.isric.org/geoserver/ows',version='1.0.0',connected=False)
    # check whether all layers are required
    if soilgridlrs is None:
        lrs_avlbl = list(wcs.contents)
        lsr_name = [lr_avlbl.split(sep=':')[1] for lr_avlbl in lrs_avlbl]
        lrs_dict = {lrs_avlbl[i]: lsr_name[i] for i in range(0, len(lsr_name))}
    else:
        lrs_dict = soilgridlrs
    # Calculate the required number of rows and columns per raster
    cols = int((geog.distance([inpt_data.total_bounds[0], 0],
                              [inpt_data.total_bounds[2], 0]))/raster_res[0])

    rows = int((geog.distance([inpt_data.total_bounds[1], 0],
                              [inpt_data.total_bounds[3], 0]))/raster_res[1])

    # start the loop over all entries
    for lr_name_server, lr_name_gdf in lrs_dict.items():
        response = wcs.getCoverage(identifier=lr_name_server,
                                   format='GeoTIFF', crs='EPSG:4326',
                                   bbox=tuple(inpt_data.total_bounds),
                                   width=cols, height=rows)
        print('Layer', lr_name_server, 'retrieved from WCS Server')
        # Create a BytesIO object.
        compressedFile = BytesIO()
        compressedFile.write(response.read())
        # Open the file from memory with rasterio and calculate stats
        with MemoryFile(compressedFile.getvalue()) as memfile:
            with memfile.open() as dataset:
                data_array = dataset.read()[0]
                polystats = rasterstats.zonal_stats(inpt_data, data_array,
                                                    affine=dataset.transform,
                                                    stats=statistical_metric,
                                                    all_touched=True,
                                                    nodata=dataset.nodata)
        # append on geopandas dataset
        inpt_data[lr_name_gdf] = [stat[statistical_metric[0]]
                                  for stat in polystats]
        print('Layer', lr_name_gdf, 'added to input dataset')

    # create output_data
    outpt_data = inpt_data
    # convert geometry back to original projection
    if src_epsg.srs != 'epsg:4326':
        outpt_data['geometry'] = outpt_data['geometry'].to_crs(
                src_epsg.srs)
        print('geometry converted back to EPSG:', src_epsg, ' for output')
    # create output shape if required
    if isinstance(output, str):
        inpt_data.to_file(output)

    return outpt_data


def apnd_from_restapi(inpt_data,
                              soil_attributes=['ORCDRC', 'BLDFIE', 'CLYPPT'],
                              soil_layerdepths=['sl1'], output=None):
    """
    use the Soilgrid Rest API to get point based data
    docu is here: https://www.isric.org/explore/soilgrids/faq-soilgrid
    some interesting datasets could be:
    APOBLDFIE=Bulk Density [g/m³]
    ORCDRC=Soil organic carbon mass content [g/kg]
    WWP=Available soil water capacity (volumetric fraction)to wilting point
    CLYPPT Percentage of Clay [%]
    """
    # test whether input is str and convert to database
    if isinstance(inpt_data, str):
        inpt_data = gpd.read_file(inpt_data)
    # check source coordinate system and convert to EPSG 4326
    src_epsg = inpt_data.crs
    if src_epsg.srs != 'epsg:4326':
        inpt_data['geometry'].to_crs(epsg=4326)
        print('geometry converted to EPSG:4326 for processing')
    # define attribute and depth strings for the query
    query_attributes = ','.join(soil_attributes)
    query_soildepths = ','.join(soil_layerdepths)
    # create empty columns in our geopandas dataframe which are filled later
    for soil_attribute in soil_attributes:
        for soil_layerdepth in soil_layerdepths:
            inpt_data[soil_attribute + '_' + soil_layerdepth] = np.nan

    # Loop over all points : 0.5s per request
    for index, sampl_geomtr in inpt_data.iterrows():
        # extract the coordinates from geometry depending on type
        if sampl_geomtr.geometry.geom_type == 'Point':
            # for point we take directly the coordinates
            centercoords = np.array(sampl_geomtr.geometry.xy)
        else:
            # for polygons we use the centroid
            centercoords = np.array(sampl_geomtr.geometry.centroid.xy)
        # request attributes from rest API
        while True:
            try:
                resp = requests.get("https://rest.soilgrids.org/query?lon="
                                    + str(float(centercoords[0])) + "&lat="
                                    + str(float(centercoords[1]))
                                    + "&attributes="
                                    + query_attributes+"&depths="
                                    + query_soildepths)
                soilgrd_results = resp.json()
                print('Soil Grid Query No.', index, ' of in total ',
                      inpt_data.shape[0], 'done ')
            except Exception as e:
                print(e)
                print('retry connection to soilgrids server')
                continue

            break

    # read out the dictionary and append datasets
        for soilgrdkey in soilgrd_results['properties']:
            if soilgrdkey in soil_attributes:
                soilgrd_attribute = soilgrd_results['properties'][soilgrdkey]['M']
                for soilgrd_lrkey in soilgrd_attribute:
                    soilgrd_lr_val = soilgrd_attribute[soilgrd_lrkey]
                    inpt_data.loc[index, soilgrdkey + '_'
                                  + soilgrd_lrkey] = soilgrd_lr_val
    print('Finished retrieving data from REST SoilGrids API')

    # create output_data
    outpt_data = inpt_data
    # convert geometry back to original projection
    if src_epsg.srs != 'epsg:4326':
        outpt_data['geometry'] = outpt_data['geometry'].to_crs(
                src_epsg.srs)
        print('geometry converted back to EPSG:', src_epsg, ' for output')
    if isinstance(output, str):
        inpt_data.to_file(output)

    return outpt_data