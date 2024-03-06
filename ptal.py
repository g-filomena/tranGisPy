import geopandas as gpd, networkx as nx

from tqdm import tqdm, tqdm_pandas
tqdm.pandas()

from .routing import walking_time_to_stops, from_nx_to_igraph
from .utilities import min_distance_geometry_gdf
import cityImage as ci

def prepare_network_for_PTAL(nodes_gdf, edges_gdf, stops, grid):
    
    nodes_gdf, edges_gdf = ci.reset_index_graph_gdfs(nodes_gdf, edges_gdf, nodeID = "nodeID")
    edges_gdf['length'] = edges_gdf.geometry.length
    graph = ci.multiGraph_fromGDF(nodes_gdf, edges_gdf, 'nodeID')
    stops['node'] = stops.geometry.apply(lambda row: min_distance_geometry_gdf(row, nodes_gdf)[1])
    
    grid['centroid'] = grid['geometry'].centroid
    grid['node'] = grid.centroid.apply(lambda row: min_distance_geometry_gdf(row, nodes_gdf)[1])
    ig_graph = from_nx_to_igraph(graph, weight = 'length')
    
    return ig_graph, stops, grid

def compute_ptal(polygons_gdf, stops, ig_graph, ptal_params = {'max_dist_bus': 500, 'max_dist_tram': 600, 'max_dist_metro' : 1000, 'max_dist_rail': 1500, 
                'walking_speed': 4000}, route_types_to_disregard = [], routeIDs_to_disregard = []):
    
    polygons_gdf = polygons_gdf.copy()
    stops = stops.copy()
    stops['wTime'] = 0.0

    max_dist_bus = ptal_params['max_dist_bus']
    max_dist_tram = ptal_params['max_dist_tram']
    max_dist_metro = ptal_params['max_dist_metro']
    max_dist_rail = ptal_params['max_dist_rail']
    
    walking_speed = ptal_params['walking_speed']
    threshold_bus = max_dist_bus*60/walking_speed
    threshold_tram = max_dist_tram*60/walking_speed
    threshold_metro = max_dist_metro*60/walking_speed
    threshold_rail = max_dist_rail*60/walking_speed
    cutoff = max_dist_rail*2
    
    # Computing Awaiting Time
    stops['k'] = [2 if route_type == 3 else 0.7 for route_type in stops['route_type']]
    stops['swt'] = 0.5*(60/stops['freq'])
    stops['awt'] = stops['swt']+stops['k']
        
    if route_types_to_disregard:
        stops = stops[~stops.route_type.isin(route_types_to_disregard)].copy()
    if routeIDs_to_disregard:
        stops = stops[~stops.route_id.isin(routeIDs_to_disregard)].copy()
    
    pt_stops_sindex = stops.sindex  
    
    def compute_polygon_ptal(centroid, centroid_nodeID):
        buffer = centroid.buffer(cutoff)
        possible_pt_stops_index = list(pt_stops_sindex.intersection(buffer.bounds))
        matches = stops.iloc[possible_pt_stops_index].copy()    
        
        if len(matches) < 1:
            return 0.0
        
        matches = matches.drop_duplicates(subset=['route_id', 'stop_id'], keep='first').reset_index(drop=True)
        targets = matches['node'].unique()
        walking_times = walking_time_to_stops(ig_graph, source = centroid_nodeID, stop_nodes = targets)     
        matches['wTime'] = [walking_times[matches.loc[ix].node] for ix in list(matches.index.values)]
        matches['tat'] = matches['awt'] + matches['wTime']
        matches = matches.sort_values(by=[ 'stop_id', 'route_id', 'wTime'])
        matches['edf'] = 30/matches['tat']   
        
        # considering trams as buses
        bus_condition = ((matches.route_type == 3) & (matches.wTime <= threshold_bus))
        tram_condition = ((matches.route_type == 0) & (matches.wTime <= threshold_tram))
        metro_condition = ((matches.route_type == 1) & (matches.wTime <= threshold_metro))
        train_condition = ((matches.route_type == 2) & (matches.wTime <= threshold_rail))
        
        matches = matches[bus_condition | tram_condition | metro_condition | train_condition]
        accessibility_index = 0.0
        groups = matches.groupby('route_type')
        
        for _, group in groups:
            group = group.sort_values(by='edf', ascending=False).reset_index(drop=True)
            if len(group) > 1:
                highest_edf = group.iloc[0]['edf']
                accessibility_index += highest_edf + 0.5 * group.iloc[1:]['edf'].sum()
            elif len(group) == 1:
                accessibility_index += group.iloc[0]['edf']
        
        return accessibility_index
    
    
    ptal_series = polygons_gdf.progress_apply(lambda row: compute_polygon_ptal(row.centroid, row.node), axis = 1)
    
    ptal_series.index = polygons_gdf.index    
    return ptal_series
    
   