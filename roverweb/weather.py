
"""
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
from ftplib import FTP
from datetime import datetime,timedelta
import numpy as np
import geopandas as gpd
import pandas as pd
import math
import time
from shapely.ops import nearest_points
import re
import os
import json
import xarray as xr
from zipfile import ZipFile
import io
import sys

def findkeys(node, kv):
    """
    from https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-dictionaries-and-lists
    """
    if isinstance(node, list):
        for i in node:
            for x in findkeys(i, kv):
               yield x
    elif isinstance(node, dict):
        if kv in node:
            yield node[kv]
        for j in node.values():
            for x in findkeys(j, kv):
                yield x

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
def update_stationlist(time_res='hourly',dbase_dir='dbase',data_category='air_tremperature'):
    """
    small function to get updated dwd station list
    
    """

    #dwd categories
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
                  'wind_synop': 'F',
                  'kl':'KL',
                  'more_precip':'RR',
                  'soil_temperature':'EB',
                  'water_equiv':'Wa',
                  'weather_phenomena':'wetter'                  
                  }

    # lets start
    print('Updating station list')
    
    #load metadata table
    dwd_datasets_meta=dwd_datasets_meta=json.load(open(dbase_dir+"\\dwd_station_meta.txt"))
        
    # create output directory if not existing
    
    if not os.path.exists(dbase_dir):
        os.makedirs(dbase_dir)
        
    #check whether we have an up-to-date-station-list-already
    try:
        stations_network_old=[s for s in os.listdir(dbase_dir) if 'dwd_station_network_'+time_res in s][0]
        datetime_network=datetime.date(datetime.strptime(re.findall('\d+',stations_network_old)[0],'%Y%m%d'))
        #update if more than 24hours
        dt_today=datetime.date(datetime.now())
        if (dt_today-datetime_network)<timedelta(days=1):
            print('DWD network list is up-to-date, no update needed')
            filename_stations=dbase_dir+'\\'+stations_network_old
            time_res_dbase=time_res
            return (filename_stations, time_res_dbase)
        else:
            print('DWD network list neeeds to be updated')
            os.remove(dbase_dir+'\\'+stations_network_old)
    except:
        print('DWD network list neeeds to be updated')
        pass
    
    
    # header
    stations_network=pd.DataFrame()
    
    # connect to ftp server and go to the folder
    
    # Connect to the Server
    server='opendata.dwd.de'
    ftp=connect_ftp(server = server,connected = False)
    #change to subfolder
    ftp.cwd('/climate_environment/CDC/observations_germany/climate/' + time_res +'/')
    
    #if data is not daily, we have to check categorywise the availability of information
    #get dwd categories
    dwd_categories=ftp.nlst()
    
    if data_category in dwd_categories:
        time_res_dbase=time_res
    elif time_res=='daily':
            #check hourly database
            time_res_dbase='hourly'
            ftp.cwd('/climate_environment/CDC/observations_germany/climate/' + time_res_dbase +'/')
            dwd_categories=ftp.nlst()
            if data_category in dwd_categories:
                print(data_category,' is not provided at the required resolution, daily_mean of hourly data used instead')
            else:
                time_res_dbase='10_minutes'
                if data_category in dwd_categories:
                    ftp.cwd('/climate_environment/CDC/observations_germany/climate/' + time_res_dbase +'/')
                    dwd_categories=ftp.nlst()  
                if data_category in dwd_categories:
                    print(data_category,' is not provided at the required resolution, daily_mean of 10_minutes data used instead')
                else:
                    print(data_category, 'not available')
                    sys.exit(1)       
    elif time_res=='hourly':
        time_res_dbase='10_minutes'
        if data_category in dwd_categories:
            ftp.cwd('/climate_environment/CDC/observations_germany/climate/' + time_res_dbase +'/')
            dwd_categories=ftp.nlst()  
        if data_category in dwd_categories:
            print(data_category,' is not provided at the required resolution, daily_mean of 10_minutes data used instead')
        else:
            print(data_category, 'not available')        
            sys.exit(1)      
    
    #loop through the subfolders  to get the station lists

    for category in dwd_categories:
        #retreive station_list   
        print('retrieve stationlist for', category)
        #try to get historical data
        try:
            dir_path='/climate_environment/CDC/observations_germany/climate/' + time_res_dbase +'/'+category+'/historical/'
            ftp.cwd(dir_path)
        except Exception as e:
            print(e, 'try to download category', category, 'from other folder')
            try:
                dir_path='/climate_environment/CDC/observations_germany/climate/' + time_res_dbase +'/'+category+'/'
                ftp.cwd(dir_path)
            except:
                print('Category', category, 'could not have been downloaded')
                pass
        #get the filename of the station list
        if time_res_dbase=='hourly':
            filename=dwd_abbr[category]+'_Stundenwerte_Beschreibung_Stationen.txt'
        if time_res_dbase=='daily':
            if dwd_abbr[category] =='wetter':
                filename=dwd_abbr[category]+'_tageswerte_Beschreibung_Stationen.txt'
            else:                    
                filename=dwd_abbr[category]+'_Tageswerte_Beschreibung_Stationen.txt'
        if time_res_dbase=='10_minutes':
            filename='zehn_min_'+dwd_abbr[category].lower()+'_Beschreibung_Stationen.txt'
        #retrieve the stationlist
        stationlist = []
        # try to retrieve file
        retrieved=False
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
        dfstations=pd.DataFrame(stationlist,columns=['STATIONS_ID','date_start','date_end','height','geo_lat','geo_lon','name'])
        #add true information to category
        dfstations[category]=True
        
        stations_network=stations_network.append(dfstations,sort=False,ignore_index=True)
        #A=[sub.split(" ") for sub in stationlist]        
    
    #replace all Na by False
    stations_network[stations_network.isna()]=0      
    #aggregate
    stations_network=stations_network.groupby(['STATIONS_ID'],as_index=False).agg('max')
    #replace zero by False in order to have pure boolean data
    stations_network.replace(0,False,inplace=True)
    #fix the error with station 14138 and 05614 and 07325, which does not have pressure cord
    stations_network.loc[stations_network.STATIONS_ID=='14138','pressure']=False
    stations_network.loc[stations_network.STATIONS_ID=='05614','pressure']=False
    stations_network.loc[stations_network.STATIONS_ID=='07325','pressure']=False
    stations_network.loc[stations_network.STATIONS_ID=='01572','pressure']=False
    #for temperature the same
    stations_network.loc[stations_network.STATIONS_ID=='14138','air_temperature']=False
    #save to database writing the time as well
    filename_stations=dbase_dir+'\\dwd_station_network_' + time_res_dbase +'_' + datetime.now().strftime('%Y%m%d')+'.csv'
    stations_network.to_csv(filename_stations,index=False)
                  
    print('Updating station list...finished')
    
    return  (filename_stations, time_res_dbase)


def nearest_direct(row, gdf1, gdf2, src_column=None):
    """Find the nearest point and return the corresponding value from specified column.
    df 1 is the origin
    df 2 is the destination
    inspired by https://automating-gis-processes.github.io/2017/lessons/L3/nearest-neighbour.html
    """
    #create a unary union
    unary_union = gdf2.unary_union    
    # Find the geometry that is closest
    nearest = gdf2['centroid'] == nearest_points(row['centroid'], unary_union)[1]
    # Get the corresponding value from df2 (matching is based on the geometry)
    value = gdf2[nearest][src_column].get_values()[0]
    return value

#inpt_df=pd.read_csv('Mueglitz-20190708_selection.csv')
def nearest_loop(row, gdf2,geometry_cols=['geo_lon','geo_lat'],src_column=None,surrounding=False):
    """
    takes longer, seems to be more precise and allows multiple locations
    surrounding: int: m around the measurement point, default is false
    """
    def haversine_distance(origin, destination):
        lon1, lat1 = origin
        lon2, lat2 = destination
        radius = 6371000 # meters
    
        dlat = math.radians(lat2-lat1)
        dlon = math.radians(lon2-lon1)
        a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
            * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        d = radius * c
        return d

    # start the main iteration
    if row.geometry.type == 'Polygon':
        point_xy = np.array((row.geometry.centroid.x,
                     row.geometry.centroid.y))
    if row.geometry.type in ['Point', 'LineString']:
        point_xy = np.array((row.geometry.x, row.geometry.y))    
    # Select most current stations datasets.
    closest = None
    closest_distance = 99999999999
    for _, station in gdf2.iterrows():
        d = haversine_distance((point_xy[0], point_xy[1]),
            (station[geometry_cols[0]], station[geometry_cols[1]]))
        if d < closest_distance:
            closest = station
            closest_distance = d
    # if surroung 
    if surrounding:
        closest1 = []
        closest_distance = closest_distance+surrounding
        i = 0
        for _, station in gdf2.iterrows():
            d = haversine_distance((point_xy[0], point_xy[1]),
                                   (station[geometry_cols[0]], station[geometry_cols[1]]))
            if d < closest_distance:
                closest1.append(station)
                i += 1
            closest = closest1
    return closest[src_column]

def dwd_age_test(date):
    # Compute timerange labels / subfolder names.

    age = (datetime.utcnow() - date).total_seconds() / 86400
    if age < 360:
        latest=True
        historic=False
    elif age >= 360 and age <= 370:
        latest=True
        historic=True
    else:
        historic=True
        latest=False
    timeranges = []
    if latest:
        timeranges.append("recent")
    if historic:
        timeranges.append("historical")
    return timeranges


#%% We start to connect to dwd server and download
def import_stations(time_res='hourly',time_format='%Y%m%d%H',
                    campaign_time=[datetime(2018,12,9), datetime(2018,12,12)],
                    data_category='air_temperature', station_ids=['00044','00091'],
                    dbase_dir='dbase', table_dir='tables',Output=True,
                    memory_save=True):
    """
    Imports stations to the existing netcdf database
    Warning: Currently, the only way to update the database for existing
    stations is to delete the entire database  
    WARNING: Import does not seems to be perfect, sometimes double import         
    """
    timeranges=['recent','historical']
    #%%load the datasets available at each timestep
    dwd_datasets_meta=dwd_datasets_meta=json.load(open(table_dir+"\\dwd_station_meta.txt"))
    #try to get a variable from the category, otherwise use interpolation of higher frequency data
    resample_frequency=None
    time_res_dbase=time_res
    try:
        dwd_datasets_meta[time_res][data_category]
    except Exception:
        if time_res=='daily':
            try:
              dwd_datasets_meta['hourly'][data_category]
              print(data_category,' is not provided at the required resolution, daily_mean of hourly data used instead')
              resample_frequency='D'
              time_res_dbase='hourly'
            except Exception:
                try: 
                    dwd_datasets_meta['10_minutes'][data_category]
                    print(data_category,' is not provided at the required resolution, daily_mean of 10_minutes data used instead')
                    resample_frequency='D'
                    time_res_dbase='10_minutes'
                except Exception:
                    print(data_category, 'not available')
                    sys.exit(1)
        if time_res=='hourly':
            try: 
                dwd_datasets_meta['10_minutes'][data_category]
                print(data_category,' is not provided at the required resolution, hourly_mean of 10_minutes data used instead')
                resample_frequency='H'
                time_res_dbase='10_minutes'
            except Exception:
                print(data_category, 'not available')
                sys.exit(1)
    
    
    #%% download from dwd if necessary
    #connect to server
    server='opendata.dwd.de'
    ftp=connect_ftp(server = server,connected = False)
    #get the mean time of the campaign
    date_mean=campaign_time[0]+(campaign_time[1]-campaign_time[0])/2   
    # load the inititial ds
    dbase_path=dbase_dir+'\\db_stations_'+time_res+'_'+data_category+'.nc'
    if os.path.exists(dbase_path):
        with xr.open_dataset(dbase_path) as dwd_dbase:
            dwd_dbase.load()
            print('Existing database imported')
         #get the non_nans stations
        current_stations=np.array(dwd_dbase[list(dwd_dbase.keys())[0]].sel(time=date_mean,method='nearest').dropna('STATIONS_ID').coords['STATIONS_ID'])
    else:
        print(dbase_path, 'does not exist, we create a new netcdf_file')
        dwd_dbase=xr.Dataset()
        current_stations=np.array((-9999)).reshape(1)
    #change directory on server
    for timerange in timeranges:
        archive_url='/climate_environment/CDC/observations_germany/climate/'+time_res_dbase+'/'+data_category+'/'+timerange       
        ftp.cwd(archive_url)
        #get the archive
        for station_id in station_ids:
            #we check whether the station is in the database with this parameter already
            if int(station_id) in current_stations:
                print('Station', station_id, 'with category', data_category,'in ',timerange,'dbase already')
                continue
            try:
                archive_name=[s for s in ftp.nlst() if station_id in s][0]
            except:
                print('No ',timerange,'data for station',station_id)
                continue
            print('Retrieving {}...'.format(archive_name))
            retrieved = False
            archive = io.BytesIO()
            # try to retrieve file
            while not retrieved:
                try:
                    ftp.retrbinary("RETR " + archive_name, archive.write)
                    retrieved = True
                except:
                    ftp=connect_ftp(server = server,connected = False)
                    ftp.cwd(archive_url)
            archive.seek(0)
            with ZipFile(archive) as myzip:
                for f in myzip.infolist():
                    # This is the data file
                    #print('zip content:', f.filename)
                    if f.filename.startswith('produkt_'):
                        product = io.StringIO(str(myzip.read(f.filename),'utf-8'))
            #get dataframe from product            
            dwd_product=pd.read_csv(product,sep=';',skipinitialspace=True)
            #get datetime
            dwd_product['time']=pd.to_datetime(dwd_product['MESS_DATUM'],format=time_format)    
            dwd_product=dwd_product.rename(columns=dwd_datasets_meta[time_res_dbase][data_category])
            dwd_product=dwd_product.reset_index()
            dwd_product=dwd_product.set_index(['time','STATIONS_ID'])
            dwd_product=dwd_product.drop(columns=['MESS_DATUM','quality_level_of_next_columns','end_of_record','index'])
            #append to database
            dwd_xr=dwd_product.to_xarray()
            #replace all values equal to -999 to nan
            for data_var in dwd_xr.data_vars:
                dwd_xr[data_var]=dwd_xr[data_var].where(dwd_xr[data_var]>-999)
            if station_id=='05009':
                print('ok')            
            #only add relevant dates if available memoryis rather small
            
            if memory_save and timerange=='historical':
                dwd_xr=dwd_xr.sel(time=slice(campaign_time[0]-timedelta(days=1),campaign_time[1]+timedelta(days=1)))
                #dwd_xr=dwd_xr.squeeze()
            
            try:
                dwd_dbase=xr.merge([dwd_dbase,dwd_xr])
            except Exception as e:
                print(e)
                print('try merging with compat=override')
                dwd_dbase=xr.merge([dwd_dbase,dwd_xr],compat='override')
            print(archive_name,' added to database')
    #upscale to required temporal resolution
    if resample_frequency is not None:
        dwd_dbase=dwd_dbase.resample(time=resample_frequency).mean(skipna=True)
        print('DWD data upscaled to',time_res,'averages')
    if Output==True:
        dwd_dbase.to_netcdf(dbase_path)
        print('Updated database' ,dbase_path)
    return dwd_dbase
#%% Start the main function
def Find_nearest_dwd_stations(inpt_data,
                         date_start='20051201',
                         date_end='20201231',
                         data_category='air_temperature',
                         temp_resolution='hourly',
                         no_of_nearest_stations=4,
                         memory_save=True,
                         Output='True'):
    """
    The Main function which is written from skretch and uses an netcdf database
    Run for each category individually
    """
    if isinstance(data_category,list):
        if len(list(data_category)) > 1:
            print(
                'Currently only one dwd category allowed, please run function multiple times for each category'
            )
            return None
    #define dwd time format for each temporal_resolution
    dwd_time_formats={'hourly':'%Y%m%d%H',
                      'daily':'%Y%m%d',
                      '10_minutes':'%Y%m%d%H%M'}
    
    dwd_time_format=dwd_time_formats[temp_resolution]
        
    #convert time to datetime
    dt_start=datetime.strptime(date_start,'%Y%m%d')
    dt_end=datetime.strptime(date_end,'%Y%m%d')
    print('Start quering data from DWD')
    #define the database folder
    pypath = os.path.dirname(os.path.abspath(__file__))
    table_dir = pypath + '\\' + 'tables'
    dbase_dir = pypath + '\\' + 'dbase'    
    #%% we check all available stations and create a valid list
    filename_stations,time_res_dbase=update_stationlist(time_res=temp_resolution,dbase_dir=table_dir,data_category=data_category)
    
    #update time format
    dwd_time_format=dwd_time_formats[time_res_dbase]
    
    #read station data
    stations_all=pd.read_csv(filename_stations, dtype={'STATIONS_ID': object})
    # delete all stations which do not cover the category
    dwd_stations=stations_all[stations_all[data_category]==True].copy()
    #correct to datetime
    dwd_stations['date_end']=pd.to_datetime(stations_all.date_end,format='%Y%m%d')
    dwd_stations['date_start']=pd.to_datetime(stations_all.date_start,format='%Y%m%d')
    # clean to stations which cover the campaign time #dt_low <= dt <= dt_high:
    dwd_stations=dwd_stations[(dwd_stations.date_start<=dt_start) & (dwd_stations.date_end>=dt_end)]
    #make a geodataframe out of it
    dwd_stations=gpd.GeoDataFrame(dwd_stations,geometry=gpd.points_from_xy(dwd_stations.geo_lon, dwd_stations.geo_lat))
    
    #loop through all rows to get the n closest points
    distances=pd.DataFrame()
    for _, station in dwd_stations.iterrows():
        distances[station.STATIONS_ID]=inpt_data.distance(station.geometry)
     
    #%% get the n stations with smallest distance and update database
    id_nearest_stations=distances.apply(lambda s: s.nsmallest(no_of_nearest_stations).index.tolist(), axis=1).values.tolist() #station ids
    #get them as unique values by sum a list of lists https://bit.ly/353iZQB
    id_dwd_stations=list(set(sum(id_nearest_stations,[])))
    
    #update the database
    db_dwd_stations=import_stations(time_res=temp_resolution,time_format=dwd_time_format,campaign_time=[dt_start,dt_end],data_category=data_category,station_ids=id_dwd_stations,dbase_dir=dbase_dir,Output=Output,table_dir=table_dir,memory_save=memory_save)
    
    #distance of nearest stattions
    dist_nearest_stations=pd.DataFrame(np.sort(distances.values)[:,:no_of_nearest_stations]).values.tolist() #distances themself
    #create new columns in the input data
    station_col_nm=list()
    for i in range(0,no_of_nearest_stations):
        station_col_nm.append(data_category+'_station_'+str(i))
    for i in range(0,no_of_nearest_stations):
        station_col_nm.append(data_category+'_distance_'+str(i))
    #create new dataframe
    distance_data=pd.concat([pd.DataFrame(id_nearest_stations).astype(int),pd.DataFrame(dist_nearest_stations)],axis=1)
    distance_data.columns=station_col_nm
    #add to main dataset
    inpt_data=pd.concat([inpt_data, distance_data.set_index(inpt_data.index)],axis=1)   
    
    return inpt_data,db_dwd_stations

def Apnd_dwd_data(inpt_data,dwd_dbase,
                         time_col='Date Time(UTC)',
                         data_time_format='%Y-%m-%d %H:%M:%S',
                         data_category='air_temperature',
                         parameters=['2m_air_temperature','2m_relative_humidity'],
                         no_of_nearest_stations=3,
                         idw_exponent=1): 
    print('Start quering data from DWD and using IDW algorithm for parameter interpolation')    
    # convert input time 
    df_times = pd.to_datetime(inpt_data[time_col], format=data_time_format)
    #correct if input parameter is str but not list
    if isinstance(parameters,str):
        parameters=[parameters]    
    #add additional columns to the inpt data
    inverse_dist=np.zeros((len(inpt_data),no_of_nearest_stations))
    for i in range (0,no_of_nearest_stations):
        inverse_dist[:,i]=(1/inpt_data[data_category+'_distance_'+str(i)]**2) 
    for parameter in parameters:
        inpt_data[parameter]=-9999   
        #calculate the result for each matrix
        result_matrix=np.zeros((no_of_nearest_stations,len(inpt_data)))            
        #loop over all the available stations
        for i in range (0,no_of_nearest_stations):
            #add station data colummn-wise
            station_data=dwd_dbase[parameter].sel(STATIONS_ID=xr.DataArray(inpt_data[data_category+'_station_'+str(i)]), time=xr.DataArray(df_times),method='nearest').values
            #http://xarray.pydata.org/en/stable/indexing.html                
            #check for nan values
            station_data[np.where(station_data==-999)]=np.nan
            #replace distance to 100000 if value is nan
            inpt_data.loc[np.isnan(station_data),data_category+'_distance_'+str(i)]=10000000
            #replace if distance is zero (can happen if device is at weather station the value with a value one digit smaller then GPS precision
            inpt_data.loc[inpt_data[data_category+'_distance_'+str(i)]==0,data_category+'_distance_'+str(i)]=0.00000001
            #add inverse distances in dependence on the data
            inverse_dist[:,i]=(1/inpt_data[data_category+'_distance_'+str(i)]**idw_exponent)
            #depending whether there is data or not, 
            result_matrix[i,:]=np.array(station_data*(1/inpt_data[data_category+'_distance_'+str(i)]**idw_exponent))
        #add to new columns
        inpt_data[parameter]=np.nansum(result_matrix,axis=0)/inverse_dist.sum(axis=1)
    print('Finished querying from DWD')
    return inpt_data

