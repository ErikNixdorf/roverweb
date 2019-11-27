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
import math
import pandas as pd
from ftplib import FTP
import time
import re
from datetime import datetime,timedelta
import os


#%% some function we need to parse our datasets
def nearest_direct(row, gdf1, gdf2, src_column=None):
    from shapely.ops import nearest_points
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

def update_stationlist(time_res='hourly',dbase_dir='dbase'):
    """
    small function to get updated dwd station list
    
    """

    def connect_ftp(server = 'opendata.dwd.de',connected=False):
        while not connected:
            try:
                ftp = FTP(server)
                ftp.login()
                connected = True
        
            except:
                time.sleep(5)
                print('Reconnect to Server')
                pass
        return ftp
    
    dwd_abbr = {'air_temperature': 'TU',
                  'cloud_type': 'CS', 
                  'cloudiness': 'N',
                  'dew_point' : 'TD',
                  'extreme_temperature': 'TX',
                  'extreme_wind': 'FX',
                  'precipitation': 'RR',
                  'pressure': 'P0',
                  'soil_temperature': 'EB',
                  'solar': 'ST',
                  'sun': 'SD',
                  'visibility': 'VV',
                  'wind': 'FF',
                  'wind_synop': 'F'
                  }
    
    # lets start
    print('Updating station list')
        
    # create output directory if not existing
    
    if not os.path.exists(dbase_dir):
        os.makedirs(dbase_dir)
    #check whether we have an up-to-date-station-list-already
    stations_network_old=os.listdir(dbase_dir)[0]
    datetime_network=datetime.date(datetime.strptime(re.findall('\d+',stations_network_old)[0],'%Y%m%d'))
    #update if more than 24hours
    dt_today=datetime.date(datetime.now())
    if (dt_today-datetime_network)<timedelta(days=1):
        print('DWD network list is up-to-date, no update needed')
        filename_stations=dbase_dir+'\\'+stations_network_old
        return filename_stations
    else:
        print('DWD network list neeeds to be updated')
        os.remove(dbase_dir+'\\'+stations_network_old)
    
    
    # header
    stations_network=pd.DataFrame()
    
    # connect to ftp server and go to the folder
    
    # Connect to the Server
    server='opendata.dwd.de'
    ftp=connect_ftp(server = server,connected = False)
    #change to subfolder
    ftp.cwd('/climate_environment/CDC/observations_germany/climate/' + time_res +'/')
    #get dwd categories
    dwd_categories=ftp.nlst()
    #loop through the subfolders  to get the station lists
    for category in dwd_categories:
        print('retrieve stationlist for', category)
        #try to get historical data
        try:
            dir_path='/climate_environment/CDC/observations_germany/climate/' + time_res +'/'+category+'/historical/'
            ftp.cwd(dir_path)
        except Exception as e:
            print(e, 'try to download category', category, 'from other folder')
            try:
                dir_path='/climate_environment/CDC/observations_germany/climate/' + time_res +'/'+category+'/'
                ftp.cwd(dir_path)
            except:
                print('Category', category, 'could not have been downloaded')
                pass
        #retrieve the stationlist
        stationlist = []
        # try to retrieve file
        retrieved=False
        filename=dwd_abbr[category]+'_Stundenwerte_Beschreibung_Stationen.txt'
        while not retrieved:
            try:
                ftp.retrlines("RETR " + filename, stationlist.append)
                #ftp.retrbinary("RETR " + filestr, stationlist.write)
                retrieved = True
            except:
                ftp=connect_ftp(server = server,connected = False)
                ftp.cwd(dir_path)
        #remove first two lines
        stationlist=stationlist[2:]
        #delete uncessary blanks
        stationlist=[re.sub(' +', ' ', station.rstrip()) for station in stationlist]
        #split the list
        stationlist=[station.split(" ")[:7] for station in stationlist]
        #read as dataframe
        dfstations=pd.DataFrame(stationlist,columns=['station_id','date_start','date_end','height','geo_lat','geo_lon','name'])
        #add true information to category
        dfstations[category]=True
        
        stations_network=stations_network.append(dfstations,sort=False,ignore_index=True)
        #A=[sub.split(" ") for sub in stationlist]        
    
    #replace all Na by False
    stations_network[stations_network.isna()]=0      
    #aggregate
    stations_network=stations_network.groupby(['station_id'],as_index=False).agg('max')
    #replace zero by False in order to have pure boolean data
    stations_network.replace(0,False,inplace=True)
    
    #save to database writing the time as well
    filename_stations=dbase_dir+'\\dwd_station_network_'+datetime.now().strftime('%Y%m%d')+'.csv'
    stations_network.to_csv(filename_stations,index=False)
                  
    print('Updating station list...finished')
    
    return  filename_stations


#%% The main function
def apnd_dwd_stationdata(inpt_data,
                         time_col='Date Time(UTC)',
                         time_format='%Y-%m-%d %H:%M:%S',
                         data_category=['air_temperature'],
                         parameters=['airtemp_humidity'],
                         temp_resolution='hourly',
                         no_of_nearest_stations=3):
    """
    The Main Function of this module
    Very Important Notice, in the dwdweather2 library I experienced some bugs
    regarding the communication with the sql database
    I fixed three main things:
        1) solved problem of vanishing category properties
        - class DwdWeather in core.py:
            # a fixed categories property
            self.categories_fixed = self.resolve_categories(category_names)
            ...
            def import_measures(self, station_id, latest=True, historic=False):
                ...
                # Download and import data.
                for category in self.categories_fixed:
        2) Solved issue that sql query fails to fetch results in core.py
            def query(self, station_id, timestamp, recursion=0):
            if recursion < 2:
                sql = "SELECT * FROM %s WHERE station_id=? AND datetime=?" % self.get_measurement_table()
                c = self.db.cursor()
                c.execute("SELECT * FROM " + self.get_measurement_table() +" WHERE station_id= ? AND datetime= ? ",(int(station_id),str(timestamp.strftime(self.get_timestamp_format())),))
        3) Solved issue of geeting incomplete list of stations
            wrote an own tool to get station list uptodate
                    
    WARNING: CURRENTLY ONLY ONE CATEGORY ACCEPTED, FOR MULTIPLE CATEGORIES
    PLEASE RUN THE FUNCTION MULTIPLE TIMES  
    
    """
    if len(data_category)>1:
        print('Currently only one dwd category allowed, please run function multiple times for each category')
        return None
    print('Start quering data from DWD')
    # Define outpt_data
    outpt_data = gpd.GeoDataFrame()
    #define the database folder
    pypath=os.path.dirname(os.path.abspath(__file__))
    dbase_dir=pypath+'\\'+'dbase'
    # next we convert our date to a datetime object
    inpt_data[time_col] = pd.to_datetime(
        inpt_data[time_col], format=time_format)
    # First we create an object from class DwDWeather
    dw = DwdWeather(resolution=temp_resolution, category_names=parameters)
    #update the sql database
    dw.import_stations()  
    # we check all available stations and create a valid list
    filename_stations=update_stationlist(time_res='hourly',dbase_dir=dbase_dir)
    stations_all=pd.read_csv(filename_stations)
    # delete all stations which do not cover the category
    dwd_stations=stations_all[stations_all[data_category[0]]==True].copy()
    #correct to datetime
    dwd_stations['date_end']=pd.to_datetime(stations_all.date_end,format='%Y%m%d')
    dwd_stations['date_start']=pd.to_datetime(stations_all.date_start,format='%Y%m%d')
    # clean to stations which cover the campaign time #dt_low <= dt <= dt_high:
    dwd_stations = dwd_stations[(dwd_stations.date_start <= inpt_data.iloc[0][time_col])
        & (inpt_data.iloc[0][time_col] <= dwd_stations.date_end)]
    #make a geodataframe out of it
    dwd_stations=gpd.GeoDataFrame(dwd_stations,geometry=gpd.points_from_xy(dwd_stations.geo_lon, dwd_stations.geo_lat))
    
    #loop through all rows to get the n closest points
    distances=pd.DataFrame()
    for _, station in dwd_stations.iterrows():
        distances[station.station_id]=inpt_data.distance(station.geometry)
     
    # get the n stations with smallest distance
    id_nearest_stations=distances.apply(lambda s: s.nsmallest(no_of_nearest_stations).index.tolist(), axis=1).values.tolist() #station ids
    dist_nearest_stations=pd.DataFrame(np.sort(distances.values)[:,:no_of_nearest_stations]).values.tolist() #distances themself
    #write it together to a dictionary
    inpt_data['nearest_stations']=[dict(zip(id_nearest_stations[i], dist_nearest_stations[i])) for i in range(0,len(id_nearest_stations))]
    
    print('Start quering data from DWD and using IDW algorithm for parameter interpolation')
    
    #add additional columns to the inpt data
    for parameter in parameters:
        inpt_data[parameter]=-9999
    #%% some older ideas to get the neighboring stations
    #inpt_data['centroid']=inpt_data.centroid
    #dwd_stations['centroid']=dwd_stations.centroid
    #inpt_data['nearest_dwd_id'] = inpt_data.apply(nearest_direct, gdf1=inpt_data.copy(), gdf2=dwd_stations.copy(), src_column='station_id', axis=1)
    
    #get the same ones by using a looped approach
    #a=time.time()
    #inpt_data['nearest_dwd_id'] = inpt_data.apply(nearest_loop,gdf2=dwd_stations,geometry_cols=['geo_lon','geo_lat'],src_column='station_id', axis=1)
    #print('Loop took', time.time()-a,'seconds')
    #http://www.gitta.info/ContiSpatVar/en/html/Interpolatio_learningObject2.xhtml
    #delete the centroid_column
    #%%
    # Define outpt_data
    outpt_data = gpd.GeoDataFrame()
    
    #loop trough all lines to get the parameterset and apply idw #http://www.gitta.info/ContiSpatVar/en/html/Interpolatio_learningObject2.xhtml
    for _,row in inpt_data.iterrows():
        # get the information of the nearest stations
        # query for result
        ii=0
        #creating result array  and inverse distance
        inverse_dist=0
        result_matrix=np.zeros((no_of_nearest_stations,len(parameters)))
        for station_id, station_dist in dict(row['nearest_stations']).items():        
            query_result = dw.query(station_id=station_id,
                                timestamp=row[time_col])
            # extract the requested parameter from the entire data              
            for i in range(0,len(parameters)):
                result_matrix[ii,i]=query_result[parameters[i]]*(1/station_dist**2)
            ii+=1
            inverse_dist+=(1/station_dist**2)
        row[parameters]=np.sum(result_matrix,axis=0)/inverse_dist
    # append datasets
        outpt_data = outpt_data.append(row)
    #delete the nearest station information
    outpt_data.drop(columns='nearest_stations',inplace=True)
    print('finished quering data from DWD')
    return outpt_data
