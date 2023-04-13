import numpy as np
from shapely.geometry import Point, LineString

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def plot_multiplex_network(multiplex_nodes, multiplex_edges, layers, nodes_sizes, edges_widths, colors, 
                           stretch = 2000, standard_theme = False):
    """
    Plots a multiplex network graph with 3D visualization.
    
    Parameters:
    - multiplex_graph (networkx.MultiGraph): A multiplex network graph object.
    - multiplex_edges (GeoDataFrame): A GeoDataFrame containing the edges of the graph.
    
    Returns:
    - matplotlib.figure.Figure: A 3D figure object containing the plotted multiplex network.
    """
    
    # Extract node coordinates and attributes
    fig_height = 40
    west, south, east, north = multiplex_edges.total_bounds
    bottom, top,  = multiplex_nodes['z'].min()*stretch, multiplex_nodes['z'].max()*stretch
    
    bbox_aspect_ratio = (north - south) / (east - west)*1.5
    fig_width = fig_height+90 / bbox_aspect_ratio/1.5
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(projection='3d')
    
    multiplex_nodes = multiplex_nodes.copy()
    nodes_x = [geometry.coords[0][0] for geometry in multiplex_nodes['geometry']]
    nodes_y = [geometry.coords[0][1] for geometry in multiplex_nodes['geometry']]
    nodes_z = [geometry.coords[0][2] for geometry in multiplex_nodes['geometry']]
       
    marker_sizes = [nodes_sizes[layers.index(layer)] for layer in multiplex_nodes['layer']]
    marker_colors = [colors[layers.index(layer)] for layer in multiplex_nodes['layer']]       
    ax.scatter(nodes_x, nodes_y, nodes_z, s=marker_sizes, c=marker_colors, zorder=0)
 
    for layer in multiplex_edges.layer.unique():
        tmp = multiplex_edges[multiplex_edges.layer == layer].copy()
        if layer != "inter_layer": 
            line_width = edges_widths[layers.index(layer)]
            line_color = colors[layers.index(layer)]
        else:
            line_width = 0.10
            line_color = 'gray'
           
        lines_tmp = [[(list(coord)) for coord in list(geo.coords)] for geo in tmp.geometry for coord in 
                     list(geo.coords)]
        lines = [[[coords[0], coords[1], coords[2]*stretch] for coords in line] for line in lines_tmp]
        lines_widths = [line_width]*len(tmp)
        ax.add_collection3d(Line3DCollection(lines, linewidths=lines_widths, alpha=1.0, color=line_color))
    
    ax.set_ylim(south, north)
    ax.set_xlim(west, east)
    ax.set_zlim(bottom, top)
    
    ax.margins(0)
    ax.tick_params(which="both", direction="in")
    ax.axis("off")
    
    if not standard_theme:

        ax.set_facecolor("black")
        ax.set_aspect("equal")

    fig.canvas.draw()

    return fig