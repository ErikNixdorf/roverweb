# RoverWeb
roverweb is a python module which append data from three different kind of web database (**OpenStreetMap**, **Soilgrids** and **DWD_GAUGES**) to geometric features. Originally aimed to retrieve ancillary data for processing cosmic-ray neutron sense records the submodule function are written more generic to append spatial-temporal data to any geometric features provided


## Submodules
The module consists of 4 submodules
### geometry.py
the geometry tool collects different geometric operations required for (at least parts) of the other modules.
e.g. points can be extented to polygons representing circles (```points_to_circle```) and sets of those circles can be merged to trapezoids (```points_to_footprint```)

code example: 
```python
import roverweb as rw

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
gdfRover_polyg=rw.geometry.points_to_footprint(gdfRover,footprint_radius=50,crs_src = "epsg:4326",crs_dst = "epsg:4326",number_of_points=10,inpt_geo_type='Point')
# Next is that we cluster our datasets in order to improve query speed
clusters = 5  # number of clusters
gdfRover_clust = rw.geometry.clustering(
    gdfRover_polyg, clusters=clusters, cluster_fld_name='ClusterID')
```

### osm.py
osm.py module retrieves data from openstreetmap via overpass API. Input needs to be polygon-like. Different ```queryfeatures``` can be selected and specific entries can ignored as well. Usually, more than one feature can be found for each polygon. Hence, the majority entry is assigned to the output geodataframe representing the geometry. In addition, if ```CountValues=True```, the sum of counts of the majority feature is provided as output entry