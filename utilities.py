import pandas as pd
import geopandas as gpd

import numpy as np
from shapely.geometry import Point, Polygon
import math

def create_grid_hexagons(gdf, side_length = 150):
    """
    Create a grid of hexagons, for a given GeoDataFrame's extent.
    
    Parameters
    ----------
    gdf: GeoDataFrame
    side_length: float
        length of the hexagon's edgeline_geometries
    
    Return
    ----------
    grid: Polygon GeoDataFrame
        the resulting grid of hexagons
    """
    xmin, ymin,xmax,ymax = gdf.total_bounds # lat-long of 2 corners
    EW = Point(xmin,ymin).distance(Point(xmax,ymin))
    NS = Point(xmin,ymin).distance(Point(xmin,ymax))

    height = int(side_length*1.73)
    width = side_length*2

    cols = list(range(int(np.floor(xmin)), int(np.ceil(xmax)), width))
    rows = list(range(int(np.floor(ymin)), int(np.ceil(ymax)), height))
    rows.reverse()
    polygons = []
    odd = False

    to_reach = cols[-1]
    x = cols[0]
    
    while (x < to_reach):
        if odd: 
            x = x-side_length/2
        for y in rows:
            if odd: 
                y = y-height/2
            centroid = Polygon([(x,y), (x+side_length, y), (x+side_length, y-side_length), (x, y-side_length)]).centroid       
            polygons.append(create_hexagon(side_length, centroid.coords[0][0], centroid.coords[0][1] ))
        if odd: 
            x = x + width-side_length/2
        else: x = x+width
        odd = not odd


    grid = gpd.GeoDataFrame({'geometry':polygons}, crs = gdf.crs)
    return grid 
 
def create_grid_squares(gdf, side_length = 200):

    xmin, ymin, xmax, ymax =  study_area.total_bounds
    rows = int(np.ceil((ymax-ymin) /  side_length))
    cols = int(np.ceil((xmax-xmin) / side_length))
    
    x_leftOrigin = xmin
    x_rightOrigin = xmin + side_length
    y_topOrigin = ymax
    y_bottomOrigin = ymax - side_length
    polygons = []
    
    for i in range(cols):
        y_top = y_topOrigin
        y_bottom =y_bottomOrigin
        for j in range(rows):
            polygons.append(Polygon([(x_leftOrigin, y_top), (x_rightOrigin, y_top), (x_rightOrigin, y_bottom), 
                                     (x_leftOrigin, y_bottom)])) 
            y_top = y_top - side_length
            y_bottom = y_bottom - side_length

        x_leftOrigin = x_leftOrigin + side_length
        x_rightOrigin = x_rightOrigin + side_length

    grid = gpd.GeoDataFrame({'geometry':polygons}, crs = gdf.crs)
    return grid 
    
def create_hexagon(side_length, x, y):
    """
    Create a hexagon centered on (x, y)
    
    Parameters
    ----------
    side_length: float
        length of the hexagon's edgeline_geometries
    x: float
        x-coordinate of the hexagon's center
    y: float
        y-coordinate of the hexagon's center
    
    Return
    ----------
    polygon: Polygon
        the resulting hexagon
    """
       
    c = [[x + math.cos(math.radians(angle)) * side_length, y + math.sin(math.radians(angle)) * side_length] for angle in range(0, 360, 60)]
    polygon = Polygon(c)
    return polygon
    
def rescale_geometry(geometry, factor):
    """
    It rescales a geometry of a certain factor. The same factor is applied to each dimension (x, y,z), from the center of the geometry.
    See shapely.affinity.scale for details
    
    Parameters
    ----------
    geometry: Polygon
        a polygon
    factor: float
        the factor/distance to use for rescaling
    
    Return
    ----------
    rescaled_geometry: Polygon
        the resulting rescaled geometry        
    """
    rescaled_geometry = scale(geometry, xfact= factor, yfact= factor, zfact=factor, origin='center') 
    return rescaled_geometry
    
def average_value_nodes_in_hexagons(hexagons_gdf, nodes_gdf, value_column_name):
    # Create a new column in the hexagons GeoDataFrame to store the average value
    hexagons_gdf[value_column_name] = 0.0
    
    # Create a spatial index for the nodes GeoDataFrame
    nodes_sindex = nodes_gdf.sindex
    
    def compute_average_value(hexagon):
        possible_matches_indexes = list(nodes_sindex.intersection(hexagon.geometry.bounds))
        possible_matches = nodes_gdf.iloc[possible_matches_indexes]
        matches = possible_matches[possible_matches.geometry.within(hexagon.geometry)]
                       
        # If there are no containing nodes, return NaN
        if len(matches) == 0:
            return float("nan")
        
        # Compute the average value for the containing nodes and return it
        return matches[value_column_name].mean()
    
    # Compute the average value for each hexagon using apply
    hexagons_gdf[value_column_name] = hexagons_gdf.apply(compute_average_value, axis=1)
        
    return hexagons_gdf
    
    
 # replace the values in the 'from' and 'to' columns with the corresponding identifier
def assign_inbound_outbound(routes_gdf):
    
    # routes_gdf['from'] = [geo.coords[0] for geo in routes_gdf.geometry]
    routes_gdf['to'] = [geo.coords[-1] for geo in routes_gdf.geometry]
    routes_gdf['nr_stops'] = [len(list(geo.coords)) for geo in routes_gdf.geometry]
    routes_gdf = assign_inbound_outbound(routes_gdf)
    
    
    routes_gdf = routes_gdf.copy()
    tuples = list(routes_gdf[['from', 'to']].drop_duplicates().values)
    tuples_dict = {tuple_[0]: i for i, tuple_ in enumerate(tuples)}
    routes_gdf['from'] = routes_gdf['from'].map(tuples_dict)
    routes_gdf['to'] = routes_gdf['to'].map(tuples_dict)
    routes_gdf['coords'] = [(list(geo.coords)) for geo in routes_gdf.geometry]
    routes_gdf["coords"] = [str(coords[::-1]) if to > fr else str(coords) for coords, to, fr in 
                            zip(routes_gdf['coords'], routes_gdf['from'], routes_gdf['to'])]
    
    routes_gdf['duplicate'] = routes_gdf.duplicated(subset=['route_id', 'coords'], keep=False)

    # assign 'a' to the row with higher value in col3
    routes_gdf.loc[routes_gdf.groupby(['route_id', 'coords'])['freq'].idxmax(), 'direction'] = 'outbound'

    # assign 'r' to the other duplicates
    routes_gdf.loc[routes_gdf['duplicate'] & ~routes_gdf['direction'].notna(), 'direction'] = 'inbound'

    # fill the NaN values with 'a'
    routes_gdf['direction'].fillna('outbound', inplace=True)
    routes_gdf.drop('duplicate', inplace = True, axis = 1)
    return(routes_gdf)
    
    
def fill_grid(grid, column):
    grid = grid.copy()
    not_empty = grid[~grid[column].isna()]
    empty = grid[grid[column].isna()]

    empty[column] = empty.geometry.apply(lambda row: not_empty.loc[min_distance_geometry_gdf(row, not_empty)[1]][column])
    grid.loc[grid[column].isna(), column] = empty[column]
    return grid
    
def min_distance_geometry_gdf(geometry, gdf):
    """
    Given a geometry and a GeoDataFrame, it returns the minimum distance between the geometry and the GeoDataFrame. 
    It provides also the index of the closest geometry in the GeoDataFrame.
    
    Parameters
    ----------
    geometry: Point, LineString or Polygon
    
    gdf: GeoDataFrame
    
    Returns:
    ----------
    distance, index: tuple
        The closest distance from the geometry, and the index of the closest geometry in the gdf.
    """
    sindex = gdf.sindex
    closest = sindex.nearest(geometry, return_distance = True)
    iloc = closest[0][1][0]
    distance = closest[1][0]
    index = gdf.iloc[iloc].name 
    return distance, index
    
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