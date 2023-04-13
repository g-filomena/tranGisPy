import igraph as ig
import networkx as nx

def paths_sources_target(graph, sources, targets, weight='length', mode = 'out'):
    
    """
    Finds the shortest paths between all pairs of sources and targets in a graph.

    Parameters:
    graph (NetworkX Graph): 
        The graph on which to perform the operation.
    sources (list): 
        A list of source nodes for which to find shortest paths.
    targets (list): 
        A list of target nodes for which to find shortest paths.
    weight (str): 
        The name of the edge attribute to be used as weight. Defaults to 'length'.

    Returns:
    list: A list of tuples representing the shortest path between each source and each target. 
          Each tuple contains three elements: source node, target node and cost of the path.

    """
    
    ig_graph = from_nx_to_igraph(graph, weight)
    
    # Calculate the shortest paths using igraph
    shortest_paths = ig_graph.distances(source=sources, target=targets, weights='weight')
    zipped =  zip(sources, shortest_paths)
    paths_list = [[(source, target, cost) for target, cost in zip(targets, costs)] for source, costs in zipped]
    paths = [item for sublist in paths_list for item in sublist]
    

    return paths
 
def from_nx_to_igraph(graph, weight):
   
   # Convert the NetworkX graph to an igraph object
    graph_int = nx.convert_node_labels_to_integers(graph)
    ig_graph = ig.Graph()
    ig_graph.add_vertices(list(graph_int.nodes()))
    ig_graph.add_edges(list(graph_int.edges()))
    ig_graph.es['weight'] = list(nx.get_edge_attributes(graph_int, weight).values())
    
    return ig_graph
    
def walking_time_to_stops(ig_graph, source, stop_nodes, walking_speed = 4000):
    walking_distances = ig_graph.distances(source = [source], target = stop_nodes, weights = 'weight', mode = 'all')
    paths_list = [[(target_node, cost) for target_node, cost in zip(stop_nodes, costs)] for costs in walking_distances]
    paths = [item for sublist in paths_list for item in sublist]
    walking_time = {path[0]: path[1]*60/walking_speed for path in paths}
    return walking_time


    
    