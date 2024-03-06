import pandas as pd
import geopandas as gpd

import numpy as np
from shapely.geometry import Point, LineString

import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 

from .utilities import min_distance_geometry_gdf

def create_multiplex_gdfs(nodes_gdfs, edges_gdfs, index_main_network, Zs, layer_labels, transfer_distances, inter_layer_speed = 4000, keep_same_vertexes_edges = False):
    
    """
    This function creates a multiplex network from a list of node and edge GeoDataFrames, where each GeoDataFrame represents a layer of the multiplex network. The function assigns z-coordinates to
    the nodes based on the given layer labels and z-values, and creates inter-layer links between the nodes based on the given transfer distances. The function also adds a new column 'm_nodeID' 
    to each node GeoDataFrame to represent the unique node ID within the multiplex network.
    
    
    """
    # assigning z coordinates and creating inter layer links
    nodes_gdfs, edges_gdfs = assign_z_coordinates(nodes_gdfs, edges_gdfs, layer_labels, Zs, False) 
    
    # adding MultiPlex nodeIDs
    for n, nodes_gdf in enumerate(nodes_gdfs):
        nodes_gdf.index = nodes_gdf['nodeID']
        if n == index_main_network: 
            nodes_gdf['m_nodeID'] = nodes_gdf.index
        else: 
            nodes_gdf['m_nodeID'] = nodes_gdf.index + nodes_gdfs[n-1].m_nodeID.max()+1
        nodes_gdfs[n] = nodes_gdf
    
    # create the multiplex dataframes
    multiplex_nodes = pd.DataFrame(columns = nodes_gdfs[0][['x', 'y', 'z', 'geometry', 'nodeID','m_nodeID', 'layer']].columns)
    multiplex_edges = pd.DataFrame(columns = edges_gdfs[0][['u', 'v', 'geometry', 'length','layer']].columns)
    
    multiplex_edges = create_inter_layer_links(nodes_gdfs, index_main_network = 0, transfer_distances = transfer_distances)  
    multiplex_edges['time'] = multiplex_edges.length*60/inter_layer_speed

    # recoding edges based on the new multiplex nodeIDs
    multiplex_nodes, multiplex_edges = assign_multiplex_nodesID_to_edges(nodes_gdfs, edges_gdfs, multiplex_edges,
                                                                         index_main_network)
    
    multiplex_nodes.drop(['nodeID', 'coordinates', 'height'], axis = 1, inplace = True)
    multiplex_nodes, multiplex_edges = reset_index_graph_gdfs(multiplex_nodes, multiplex_edges, nodeID = "m_nodeID")
    return multiplex_nodes, multiplex_edges

def assign_z_coordinates(nodes_gdfs, edges_gdfs, layer_labels, Zs, keep_same_vertexes_edges = False):
    
    for n, nodes_gdf in enumerate(nodes_gdfs):
                
        nodes_gdf['z'] = Zs[n]
        edges_gdf = edges_gdfs[n].copy()
        nodes_gdf['layer'], edges_gdf['layer'] = layer_labels[n], layer_labels[n]
        nodes_gdf['geometry'] = nodes_gdf.geometry.apply(lambda row: Point(row.coords[0][0], row.coords[0][1], Zs[n]))
        edges_gdf['geometry'] = edges_gdf.geometry.apply(lambda row: LineString([(coor[0], coor[1], Zs[n]) for coor in row.coords]))
        

        if not keep_same_vertexes_edges:
            edges_tmp = edges_gdf.copy()
            
            if 'freq' not in edges_gdf.columns:
                # sort values in 'u' and 'v' columns so that "5"-"6" and "6"-"5" are considered the same
                edges_tmp[['u', 'v']] = pd.DataFrame(np.sort(edges_tmp[['u', 'v']].values, axis=1), index = edges_tmp.index)
            else:
                # create a new column 'key' based on the count of each unique pair of values in 'u' and 'v' columns
                edges_gdf['key'] = edges_gdf.groupby(['u', 'v']).cumcount()
                # sum 'freq' column for each unique pair of values in 'u' and 'v' columns
                freq_mapping = edges_gdf.groupby(['u', 'v'])['freq'].sum()
                # update 'freq' column for edges with 'key' = 0 to be the sum of 'freq' for the corresponding unique pair
                edges_gdf.loc[edges_gdf.key == 0, 'freq'] = edges_gdf[edges_gdf.key == 0][['u', 'v']].apply(lambda row: freq_mapping.loc[(row['u'], row['v'])], axis=1).values    

            # two services operating the same route between two stops or stations (not interested in frequencies here) 
            edges_gdf = edges_gdf[edges_gdf.key == 0] 
        
        nodes_gdfs[n] = nodes_gdf
        edges_gdfs[n] = edges_gdf
        
    return nodes_gdfs, edges_gdfs

def create_inter_layer_links(nodes_gdfs, index_main_network, transfer_distances):
    
    edges_list = []
    nodes_main_network = nodes_gdfs[index_main_network]
    
    for n, nodes_gdf in enumerate(nodes_gdfs): 
        if n == index_main_network: 
            continue
        
        nodes_gdf = nodes_gdf.copy()
        nodes_main_network = nodes_gdfs[index_main_network].copy()
        dist = [min_distance_geometry_gdf(node_geometry, nodes_main_network) for node_geometry in nodes_gdf.geometry]
        closest = [(d[0], d[1]) for d in dist]
        nodes_gdf[['dist', 'id_main']] = pd.DataFrame(closest, index = nodes_gdf.index)
        
        uID = [nodes_main_network.loc[id_main]['m_nodeID'] for id_main in nodes_gdf['id_main']]
        vID = list(nodes_gdf['m_nodeID'])
        dist = list(nodes_gdf['dist']+[transfer_distances[n]]*len(nodes_gdf))
        
        fictious_lines = [LineString([nodes_main_network.loc[other_id].geometry, geo]) for other_id, geo in 
                          zip(nodes_gdf.id_main, nodes_gdf.geometry)]
        
        layer = ['inter_layer']*len(fictious_lines)
        key = [0]*len(fictious_lines)
        multiplex_data = {
            'u': uID,
            'v': vID,
            'geometry': fictious_lines,
            'length': dist,
            'layer':  layer,
            'key': key,
        }
        edges_list.append(pd.DataFrame.from_records(multiplex_data))
        
    return pd.concat(edges_list, ignore_index=True)

def assign_multiplex_nodesID_to_edges(nodes_gdfs, edges_gdfs, multiplex_edges, index_main_network):
    edges_gdfs_recoded = []
    
    for n, edges_gdf in enumerate(edges_gdfs):
        if n == index_main_network:
            edges_gdfs_recoded.append(edges_gdf)
            continue
        
        edges_gdf = edges_gdf.copy()
        edges_gdf.rename(columns={'u': 'old_u', 'v': 'old_v'}, inplace=True)
        
        for col in ['u', 'v']:
            edges_gdf = edges_gdf.merge(
                nodes_gdfs[n][['nodeID', 'm_nodeID']].rename(columns={'nodeID': f'old_{col}', 'm_nodeID': col}),
                how='left', on=f'old_{col}'
            )
            edges_gdf.drop(f'old_{col}', axis=1, inplace=True)
        
        edges_gdfs_recoded.append(edges_gdf)
    
    multiplex_edges = pd.concat(edges_gdfs_recoded+[multiplex_edges], ignore_index=True)
    multiplex_edges['edgeID'] = multiplex_edges.index
    multiplex_nodes = pd.concat(nodes_gdfs, ignore_index=True)
    multiplex_nodes.index = multiplex_nodes['m_nodeID']
    multiplex_nodes.index.name = None
    
    multiplex_nodes = gpd.GeoDataFrame(data = multiplex_nodes, geometry = multiplex_nodes['geometry'], crs = nodes_gdfs[0].crs)
    multiplex_edges = gpd.GeoDataFrame(data = multiplex_edges, geometry = multiplex_edges['geometry'], crs = nodes_gdfs[0].crs)

    return multiplex_nodes, multiplex_edges


def waiting_time(multiplex_nodes, multiplex_edges):
    multiplex_nodes = multiplex_nodes.copy()
    multiplex_edges = multiplex_edges.copy()

    condition_nodes = ~multiplex_nodes.layer.isin(['ped', 'inter_layer'])
    condition_edges = ~multiplex_edges.layer.isin(['ped', 'inter_layer'])
    transport_edges = multiplex_edges[condition_edges]
    
    multiplex_nodes.loc[condition_nodes,'freq'] = [transport_edges[(transport_edges.u == node)|
                                                  (transport_edges.v == node)]['freq'].sum() 
                                                  for node in multiplex_nodes[condition_nodes].m_nodeID]
    # one direction only
    multiplex_edges.loc[condition_edges,'freq_uH'] = [multiplex_nodes.loc[u].freq/2 for u in 
                                                           multiplex_edges[condition_edges].u]

    multiplex_edges.loc[condition_edges,'waitTime'] = multiplex_edges[condition_edges]['freq_uH']/60 
    multiplex_edges['waitTime'] = multiplex_edges['waitTime'].fillna(0.0)
    multiplex_edges['time_wt'] = multiplex_edges['time']+multiplex_edges['waitTime']
    
    return multiplex_edges
    
def reset_index_graph_gdfs(nodes_gdf, edges_gdf, nodeID = "nodeID"):
    """
    The function simply resets the indexes of the two dataframes.
     
    Parameters
    ----------
    nodes_gdf: Point GeoDataFrame
        Nodes (junctions) GeoDataFrame.
    edges_gdf: LineString GeoDataFrame
        Street segments GeoDataFrame.
   
    Returns
    -------
    nodes_gdf, edges_gdf: tuple
        The junction and street segment GeoDataFrames.
    """

    edges_gdf = edges_gdf.rename(columns = {"u":"old_u", "v":"old_v"})
    nodes_gdf["old_nodeID"] = nodes_gdf[nodeID].values.astype("int64")
    nodes_gdf = nodes_gdf.reset_index(drop = True)
    nodes_gdf[nodeID] = nodes_gdf.index.values.astype("int64")
    
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[["old_nodeID", nodeID]], how="left", left_on="old_u", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {nodeID:"u"})
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[["old_nodeID", nodeID]], how="left", left_on="old_v", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {nodeID:"v"})

    edges_gdf.drop(["old_u", "old_nodeID_x", "old_nodeID_y", "old_v"], axis = 1, inplace = True)
    nodes_gdf.drop(["old_nodeID", "index"], axis = 1, inplace = True, errors = "ignore")
    edges_gdf = edges_gdf.reset_index(drop=True)
    edges_gdf["edgeID"] = edges_gdf.index.values.astype(int)
        
    return nodes_gdf, edges_gdf