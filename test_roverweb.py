# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 17:42:23 2019
A small example of how to use the code of roverweb
@author: nixdorf
"""


from roverweb import *
import roverweb
# load our dataset as geodataframe
gdfRover = roverweb.geometry.from_csv('./testdata/Mueglitz-20190708.csv',
                                      geomtrcol=['LongDec', 'LatDec'],
                                      src_crs={'init': 'epsg:4326'})
# create a circular footprint (polygon) around the gps points
gdfRover_polyg = roverweb.geometry.create_circular_footprints(gdfRover,
                                                              crcl_radius=100,
                                                              number_of_points=10)
# Next is that we cluster our datasets in order to improve query speed