import pandas as pd
import geopandas as gpd

import numpy as np
from shapely.geometry import Point, LineString, Polygon
from zipfile import ZipFile, Path
import datetime
import matplotlib as mpl
from shapely.wkt import loads
import cityImage as ci
import math

def GTFS_to_gpd(companies, crs, input_path, filter_service_strings = [], route_shape_geometries = False, 
                route_shapes_from_stops = False, times = []):
    """
    Reads in GTFS data for multiple companies and returns two GeoDataFrames:
    one containing information on routes/trips, and another containing information
    on stops along those routes.
    
    Args:
        companies (list): A list of strings, where each string is the name of a GTFS zip file.
        filter_service_strings (list): A list of strings to filter trips by, one for each GTFS zip file in companies.
        crs (str): The coordinate reference system to project the output GeoDataFrames to.
        shape_geometries (bool): A flag indicating whether to use shapes.txt to form geometries.
    
    Returns:
        routes_gdf (GeoDataFrame): A GeoDataFrame containing information on routes/trips.
        stops_gdf (GeoDataFrame): A GeoDataFrame containing information on stops along those routes.
    """
    # Loop over each GTFS file in the input list
    for n, company in enumerate(companies):
        print("importing", company)
        
         # Open the GTFS file as a ZipFile object
        with ZipFile(input_path+company+".zip") as gtfs_zip:

            # Load each relevant table from the GTFS data as a pandas DataFrame
            routes = pd.read_csv(gtfs_zip.open("routes.txt"), low_memory=False)
            trips = pd.read_csv(gtfs_zip.open("trips.txt"), low_memory=False)
            try:
                calendar = pd.read_csv(gtfs_zip.open("calendar.txt"), low_memory=False)
            except:
                calendar = pd.DataFrame()
                
            stops = pd.read_csv(gtfs_zip.open("stops.txt"), low_memory=False)
            stop_times = pd.read_csv(gtfs_zip.open("stop_times.txt"), low_memory=False)
            
            # Add a 'geometry' column to the stops table, which contains shapely Point objects
            stops['geometry'] = gpd.points_from_xy(stops['stop_lon'], stops['stop_lat'])
            
            # Sort the stop_times table and select only the relevant columns from the routes table
            stop_times = stop_times.sort_values(by = ['trip_id','stop_sequence', 'stop_id'])
            columns = ['route_id','route_short_name', 'route_long_name', 'route_type']
            routes = routes[columns]
            routes['company'] = company
            
            # Merge the trips and routes tables together and select only the relevant columns from the stop_times and stops tables
            trips = trips.merge(routes, on = 'route_id')
            stop_times = stop_times.merge(trips[['trip_id', 'route_type']], on = 'trip_id')
            stop_times = stop_times.merge(stops[['geometry', 'stop_id','stop_name']], on = 'stop_id')
            
            # peak time
            if len(times) == 2:
                stop_times = stop_times[(stop_times['arrival_time'] >= times[0]) & (stop_times['arrival_time'] < times[1])]
                trips = trips[trips.trip_id.isin(stop_times.trip_id.unique())]
            
            trips, stop_times = filter_week_days_trips(trips, calendar, stop_times = stop_times, filter_strings = filter_service_strings[n])            
            trips = assign_duration(trips, stop_times)                              
            
            # If shape geometries, form geometries from the shapes.txt file.
            if (route_shape_geometries):
                shapes = pd.read_csv(gtfs_zip.open("shapes.txt"), low_memory=False)
                trips = form_geometries_from_shapes(trips, shapes)
            
            # Otherwise, form geometries from stops
            else:         
                trips = form_geometries_from_stops(trips, stop_times)
           
            trips['wkt'] = trips['geometry'].apply(lambda x: x.wkt)
            trips = trips[['route_id', 'service_id', 'trip_id', 'shape_id', 'wkt', 'duration']].copy()
               
            # frequency trips:  
            if len(trips.trip_id.unique()) == len(trips.service_id.unique()):
                trips = trips.groupby(['route_id', 'wkt', 'service_id'], as_index = False).agg({'trip_id': 'count','duration' : 'mean'})
            else:
                trips = trips.groupby(['route_id', 'wkt'], as_index = False).agg({'trip_id': 'count', 'duration' : 'mean'})
            
            trips = trips.merge(routes, on = 'route_id') # reattaching route_type info
            trips['geometry'] = trips.apply(lambda x: loads(x['wkt']), axis = 1)
            trips.drop('wkt', axis = 1, inplace = True)
            trips.rename(columns = {'trip_id' : 'freq'}, inplace = True)
            
            stops_tmp = stop_times[['stop_id', 'stop_name', 'geometry', 'route_type']].drop_duplicates(subset=['stop_id', 'route_type'], keep = 'first')
           
            if (n == 0):
                trips_summary = trips.copy()
                stops_df = stops_tmp.copy()
            else:
                trips_summary = pd.concat([trips_summary, trips], ignore_index=True)
                stops_df = pd.concat([stops_df, stops_tmp], ignore_index=True)
                
    routes_gdf = gpd.GeoDataFrame(trips_summary, geometry= trips_summary['geometry'], crs = 'EPSG:4326')
    routes_gdf = routes_gdf.to_crs(crs)
    routes_gdf['length'] = routes_gdf.geometry.length
    
    routes_gdf.drop('service_id', inplace = True, axis = 1, errors = 'ignore')
    stops_gdf = gpd.GeoDataFrame(stops_df, geometry = stops_df['geometry'], crs = 'EPSG:4326')
    stops_gdf = stops_gdf.to_crs(crs)
    
    return routes_gdf, stops_gdf
   
def GTFS_to_pd(companies, input_path, filter_service_strings = [], times = [], polygon = None, crs = None):
    
    for n, company in enumerate(companies):
        print(company)
        
        with ZipFile(input_path+company+".zip") as gtfs_zip:
            routes_tmp = pd.read_csv(gtfs_zip.open("routes.txt"), low_memory=False)
            trips_tmp = pd.read_csv(gtfs_zip.open("trips.txt"), low_memory=False)
            
            try:
                calendar = pd.read_csv(gtfs_zip.open("calendar.txt"), low_memory=False)
            except:
                calendar = pd.DataFrame()

            stops_tmp = pd.read_csv(gtfs_zip.open("stops.txt"), low_memory=False)
            stops_tmp['geometry'] = gpd.points_from_xy(stops_tmp['stop_lon'], stops_tmp['stop_lat'])
            
            stop_times_tmp = pd.read_csv(gtfs_zip.open("stop_times.txt"), dtype={'stop_id': 'str'}, low_memory=False)
            stop_times_tmp = stop_times_tmp.sort_values(by = ['trip_id', 'stop_sequence', 'stop_id'])
            
            # peak time: select services
            if len(times) == 2:
                stop_times_tmp = stop_times_tmp[(stop_times_tmp['arrival_time'] >= times[0]) & (stop_times_tmp['arrival_time'] < times[1])]
                trips_tmp = trips_tmp[trips_tmp.trip_id.isin(stop_times_tmp.trip_id.unique())]            
            
            trips_tmp = trips_tmp.merge(routes_tmp, on = 'route_id')
            stop_times_tmp = stop_times_tmp.merge(trips_tmp[['trip_id','route_id', 'route_type']], on = 'trip_id')
           
            trips_tmp, stop_times_tmp = filter_week_days_trips(trips_tmp, calendar, stop_times = stop_times_tmp, 
                                                               filter_strings = filter_service_strings[n])                       

            
            trips_tmp = trips_tmp[['route_id', 'service_id', 'trip_id', 'direction_id','route_type']]
            stops_tmp = stops_tmp[['stop_id', 'stop_name', 'geometry']]   
            stops_tmp['stop_id'] = stops_tmp['stop_id'].astype(str)
            stop_times_tmp['stop_id'] = stop_times_tmp['stop_id'].astype(str)
            
            if polygon is not None:
                trips_tmp, stop_times_tmp = trips_within_area(trips_tmp, stops_tmp, stop_times_tmp, polygon, crs)
                trips_tmp, stop_times_tmp = urban_services(trips_tmp, stop_times_tmp) 
            
            stop_times_tmp = frequencies_at_stops(stop_times_tmp)
            stop_times_tmp = stop_times_tmp.merge(stops_tmp, on = 'stop_id') 
            stop_times_tmp = stop_times_tmp[['freq', 'route_id', 'route_type', 'stop_id', 'stop_name', 'geometry']]  
            trips_tmp['agency'] = company
            stop_times_tmp['agency'] = company
            
        if (n == 0):
            trips = trips_tmp.copy()
            stop_times = stop_times_tmp.copy()
        else:
            trips = pd.concat([trips, trips_tmp], ignore_index=True)
            stop_times = pd.concat([stop_times, stop_times_tmp], ignore_index=True)    
    
    return trips, stop_times


def assign_duration(trips, stop_times):
    
    stop_times = stop_times.copy()
    
    def trip_duration(departure_time, arrival_time):
        dh, dm,_ = map(int, departure_time.split(':'))
        ah, am,_ = map(int, arrival_time.split(':'))

        if ah < dh:
            ah += 24

        return (ah*60+am)-(dh*60+dm)
    
    grouped = stop_times.groupby(['trip_id', 'stop_sequence'], as_index=False).agg({'arrival_time': 'first', 'departure_time': 'first'})
    grouped = grouped.sort_values(['trip_id', 'stop_sequence'])
    times = pd.concat([grouped.groupby('trip_id', as_index=False).first(), grouped.groupby('trip_id', as_index=False).last()]
                            ).sort_values(['trip_id', 'stop_sequence'])
    times_sy = times.groupby('trip_id', as_index=False).agg({'departure_time': 'first', 'arrival_time': lambda x: x.iloc[1]})
    times_sy['duration'] = times_sy.apply(lambda row: trip_duration(row.departure_time, row.arrival_time), axis=1)
    trips = trips.merge(times_sy[['duration', 'trip_id']], on='trip_id')
    
    return trips

def form_geometries_from_stops(trips, stop_times):
    
    stop_times_shapes = stop_times[['trip_id', 'stop_sequence', 'geometry']].copy()
    stop_times_shapes = stop_times_shapes.sort_values(by = ['trip_id', 'stop_sequence'])
    grouped = stop_times_shapes.groupby('trip_id', as_index = False).agg({'geometry': list})
    grouped = grouped[grouped.geometry.str.len() > 1]
    grouped['geometry'] = [LineString(point) for point in grouped['geometry']]
    trips = trips.merge(grouped, on = 'trip_id')
    return trips
    
def form_geometries_from_shapes(trips, shapes):   
    
    shapes['geometry'] = gpd.points_from_xy(shapes['shape_pt_lon'], shapes['shape_pt_lat'])
    shapes = shapes.sort_values(by = ['shape_id', 'shape_pt_sequence'])
    grouped = shapes.groupby('shape_id', as_index = False).agg({'geometry': list})
    grouped['geometry'] = [LineString(point) for point in grouped['geometry']]
    trips = trips.merge(grouped, on = 'shape_id')
    return trips
       
def filter_week_days_trips(trips, calendar = pd.DataFrame(), stop_times = None, filter_strings = []):
    
    if len(calendar) > 0:   
        week_services = calendar[calendar.wednesday == 1].service_id.unique()
        return trips[trips.service_id.isin(week_services)], stop_times
    
    elif len(filter_strings)>0:
        
        trips = trips.sort_values(by = 'route_type', ascending = True)
        route_types = trips.route_type.unique()
        ids = []
        for n, route_type in enumerate(route_types):
            stop_times_filter = stop_times[stop_times.route_type == route_type].copy()
            mask = stop_times_filter['trip_id'].astype(str).str.contains('|'.join(filter_strings[n]), case=False)
            stop_times_filter = stop_times_filter[mask]
            ids.extend(list(stop_times_filter.trip_id.unique()))
        
        return trips[trips.trip_id.isin(ids)], stop_times[stop_times.trip_id.isin(ids)]
    
    return trips, stop_times

def frequencies_at_stops(stop_times):

    groups = stop_times.groupby(['stop_id', 'route_id', 'route_type'], as_index = False)[['trip_id']].count()
    groups.rename(columns = {'trip_id' : "freq"}, inplace = True)   
    return groups

def urban_services(trips, stop_times):

    # filter out service that only stop once within the study_area
    single_trips_id = stop_times.groupby('trip_id').size().eq(1).reset_index(name='count')
    single_trips_id = single_trips_id[single_trips_id['count']].trip_id.unique()
    stop_times = stop_times[~stop_times.trip_id.isin(single_trips_id)]
    trips = trips[trips.trip_id.isin(stop_times.trip_id.unique())]
    
    return trips, stop_times

# keep stops/stations in study area
def trips_within_area(trips, stops, stop_times, polygon, crs):
    
    stops_gdf = gpd.GeoDataFrame(stops, geometry = stops['geometry'], crs = 'EPSG:4326')
    stops_gdf = stops_gdf.to_crs(crs)
    stops_gdf = stops_gdf[stops_gdf.geometry.within(polygon)]
    
    stop_times = stop_times[stop_times.stop_id.isin(stops_gdf.stop_id.unique())]
    trips = trips[trips.trip_id.isin(stop_times.trip_id.unique())].copy()
    
    return trips, stop_times

def split_routes_at_stops(routes_gdf):

    routes_gdf['route_length'] = routes_gdf.geometry.length
    
    def split_route_at_stops(row):
        coords = list(row['geometry'].coords)
        segments = [LineString([coords[i], coords[i+1]]) for i in range(len(coords)-1)]
        to_populate = [row] * len(segments)
        df = pd.DataFrame(to_populate, columns=routes_gdf.columns)
        df['geometry'] = segments
        return df
    
    routes_gdfs = routes_gdf.apply(lambda row: split_route_at_stops(row), axis=1)
    routes_df_split = pd.concat(list(routes_gdfs), ignore_index=True)
    routes_gdf_split = gpd.GeoDataFrame(data = routes_df_split, geometry = 'geometry', crs = routes_gdf.crs)
    routes_gdf_split['length'] = routes_gdf_split['geometry'].length
    
    return routes_gdf_split

    
    