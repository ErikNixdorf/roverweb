# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:27:45 2019
Some functions which map the osm data on the input gpd
@author: nixdorf, UFZ
"""
import numpy as np
from scipy.spatial import ConvexHull
import geopandas as gpd
import overpass  # Server Communication with OSM


def getPolyCoords(gdf_in, row):
    """Returns the coordinates ('x|y') of edges/vertices of a Polygon/others
    from
    https://stackoverflow.com/questions/55659835/trying-to-
    separate-polygon-data-into-x-and-y-coordinates-but-get-error-multip


    """

    # Parse the geometries and grab the coordinate
    geometry = gdf_in.iloc[row].geometry
    # print(geometry.type)

    if geometry.type == 'Polygon':
            # Get the x coordinates of the exterior
            # Interior is more complex: xxx.interiors[0].coords.xy[0]
            return list(geometry.exterior.coords.xy)

    if geometry.type in ['Point', 'LineString']:
            return list(geometry.xy)
    if geometry.type == 'MultiLineString':
        all_xy = []
        for ea in geometry:
            all_xy.append(list(ea.xy))
        return all_xy

    if geometry.type == 'MultiPolygon':
        all_xy = []
        for ea in geometry:
                all_xy.append(list(ea.exterior.coords.xy))
        return all_xy

    else:
        # Finally, return empty list for unknown geometries
        return []




def querybound_generation(gdf_in):
    """
    This functions aims to generate the queryboundstring which is needed for
    the following query
    Input must be a geodataframe
    """
    # Create the polycoordinats using list comprehension
    polycoords = np.vstack([np.array(getPolyCoords(gdf_in, row)).T
                            for row in range(0, gdf_in.shape[0])])
    hull = ConvexHull(polycoords)
    # get the points of the polyline, latitude first
    hullpoints = np.transpose(np.vstack((np.append(polycoords[hull.vertices, 1],
                                                   polycoords[hull.vertices[0], 1]
                                                   ),
                                        np.append(polycoords[hull.vertices, 0],
                                                  polycoords[hull.vertices[0], 0]
                                                  ))))
    # We create a polygonstatement which we send to the query
    hullpoint_str = np.array2string(hullpoints.flatten(), precision=5)
    # the API is very sensitive if there is a space after the last number
    if hullpoint_str[-2] == ' ':
        hullpoint_str = hullpoint_str[:-3]+']'

    querybound = 'poly:"' + hullpoint_str[1:-1]+'"'

    return querybound


def apnd_from_overpass(gdfin, querybound,
                    queryfeatures={'way': ['landuse', 'highway']},
                    value_in=None, value_out=None, CountValues=False):
    """
    Get the queried features from the Overpass API
    # for documentation see https://wiki.openstreetmap.org/wiki/Overpass_API
    Currently supported for polygons only
    """

    # Connect to server
    api = overpass.API(timeout=2000)

    # get data from OSM and count server request try'"'
    for element, _ in queryfeatures.items():
            for key in queryfeatures[element]:

                # remove index_right if occuring
                try:
                    gdfin.drop('index_right', axis=1)
                except Exception:
                    pass
                print('start to query data for osm element:', element,
                      'and key: ', key, 'from overpass')
                while True:
                    try:
                        osm_retrieved = api.get(element + '[' + key + ']('
                                                + querybound + ');(._;>;);',
                                                verbosity='geom')
                    except Exception as e:
                        print(e)
                        print('retry connection to server')
                        continue

                    break
                # if no features were retrieved write nan and continue
                if len(osm_retrieved['features']) == 0:
                    print('no data retrieved for osm type:', element,
                          'and key ', key,)
                    gdfin[element[0] + '_' + key] = np.nan
                    continue
                # convert to geodataframe
                gdfosm = gpd.GeoDataFrame.from_features(osm_retrieved,
                                                        crs={'init':
                                                             'epsg:4326'})
                if element is 'way':
                    # delete all points
                    gdfosm = gdfosm[gdfosm.geom_type != "Point"]

                # drop all columns except of the relevant key
                retained_columns = ['geometry', key]
                gdfosm.drop(gdfosm.columns.difference(retained_columns),
                            1, inplace=True)

                # delete entries if defined in value_out, e.g the service
                if value_out is not None:
                    gdfosm = gdfosm[gdfosm[key] != value_out]
                # if value in is defined we look for this value only in the key
                if value_in is not None:
                    gdfosm = gdfosm[gdfosm[key] == value_in]
                    # in this case we have to rename the colfrom key to value
                    gdfosm = gdfosm.rename(columns={key: value_in})
                    key = value_in  # bad solution...
                # add ID column to input data if not existing
                if 'ID' not in gdfin.columns:
                    gdfin.insert(gdfin.shape[1], 'ID',
                                 range(0, gdfin.shape[0]))
                # join with inout gdf
                gdfjoined = gpd.tools.sjoin(gdfin, gdfosm, how='left')

                # if the number of entries larger than the No of TrackPoints,
                # we need to sum and find maxmimum
                if gdfjoined.shape[0] > gdfin.shape[0]:

                    # group by unique ones https://stackoverflow.com/questions/
                    # 36174624/
                    s = gdfjoined.groupby(['ID', key]).size()
                    # check whether values of the key or count is needed
                    if CountValues:
                        dfselection = s.reset_index().drop(key, axis=1)
                        dfselection = dfselection.rename(columns={0: 'No_of_'
                                                                  + key})
                        gdfmerged = gdfin.join(dfselection.set_index('ID'),
                                               on='ID')
                        gdfmerged['No_of_'+key].fillna(0, inplace=True)
                    else:
                        dfselection = s.loc[s.groupby(level=0).idxmax()
                                            ].reset_index().drop(0, axis=1)
                        # replace values in the original dataset
                        # https://stackoverflow.com/questions/53262894/
                        gdfmerged = gdfin.join(dfselection.set_index('ID'),
                                               on='ID')
                else:  # if it is equal things are more easy
                    if CountValues:
                        gdfjoined.rename(columns={key: 'No_of_'+key})
                        try:
                            gdfjoined['No_of_'+key].fillna(0, inplace=True)
                        except Exception:
                            print('No values type', key, 'queried')
                    gdfmerged = gdfjoined
                # remove index_right if occuring
                try:
                    gdfmerged.drop('index_right', axis=1, inplace=True)
                except Exception:
                    pass
                gdfin = gdfmerged.rename(columns={key: element[0] + '_' + key})
                print('data for osm type:', element, 'and key ',
                      key, 'retrieved')
    gdfout = gdfin
    return gdfout
