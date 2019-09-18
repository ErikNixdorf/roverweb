# RoverWeb
** Information will come soon**
Radohydro is a small python tool which converts and projects time series precipitation records from DWD Radolan data on polygon shapefiles (e.g. representing subcatchments/catchments)

## Features
* Automatize download and unzip from dwd ftp server for given time interval
* Automatic reprojection to destination EPSG given by shapefile
* Equal-area approach to compute precipitation rates per basin per time
* Output 1) as csv files containing time series for each basin and <br/>2) as polygonshapes having the precipitation rates as attribute tables

## Dependencies

* Python 3
* numpy
* geopandas
* gdal
* shapely

## Limitations
<span style="color:red">**Important: Currently fiona and rasterio package cant be installed together due to a version conflict with gdal (gdal.2.x vs gdal 3.0). As soon as this is fixed radohydro 0.4 will be released which will have less lines of code and higher performance by applying rasterstats package**</span>.

* until now supports only download of hourly datasets
* for more than 500 precipation records(~21days), output as shp is not supported (dbase problem)
* large amounts of polygons per shapefile (>10000) may lead to considerable performance decrease
* required radolan files are downloaded again every time the script runs

```python
import radohydro
radohydro.processradolan('2017123020', '2018010220','einzugsgebiet.shp', output=True)
```

## Authors

* **Erik Nixdorf**
* **Marco Hannemann**



## Acknowledgments

* Thx to Nico Trauth for improving the code


