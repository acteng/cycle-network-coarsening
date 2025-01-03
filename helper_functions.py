import pyproj
import geopandas as gpd
import pandas as pd
from shapely.geometry import box, Point
from shapely import wkt
from functools import partial
from shapely.ops import transform
from shapely.ops import split
import networkx as nx
import heapq
from sortedcontainers import SortedDict
import math
import numpy as np
from scipy.spatial import ConvexHull
import ast
from shapely import to_geojson
from geojson_length import calculate_distance, Unit

def expand_bbox(original_bbox, expansion_distance_km=5):
    # Create a Shapely geometry object for the original bounding box
    original_geometry = box(*original_bbox)
    # Define a function to project the geometry to a new coordinate reference system
    project = partial(
        pyproj.transform,
        pyproj.Proj(init='epsg:4326'),  # WGS 84 coordinate reference system
        pyproj.Proj(proj='utm', zone=33, ellps='WGS84')  # Example: UTM Zone 33
    )
    # Project the original geometry to the new coordinate reference system
    projected_geometry = transform(project, original_geometry)
    # Calculate the expansion distance in the projected coordinate system
    expansion_distance_meters = expansion_distance_km * 1000
    # Expand the geometry by the specified distance
    expanded_geometry = projected_geometry.buffer(expansion_distance_meters)
    # Project the expanded geometry back to the original coordinate reference system
    expanded_geometry = transform(partial(pyproj.transform, pyproj.Proj(proj='utm', zone=33, ellps='WGS84'), pyproj.Proj(init='epsg:4326')), expanded_geometry)
    # Get the coordinates of the expanded bounding box
    expanded_bbox = expanded_geometry.bounds
    return expanded_bbox, expanded_geometry

def create_bounding_box(geometry1, geometry2):
    # Calculate the union of all polygons in each multipolygon
    union_geometry1 = geometry1.convex_hull
    union_geometry2 = geometry2.convex_hull
    # Calculate the union of the convex hulls of the two multipolygons
    union_geometry = union_geometry1.union(union_geometry2)
    # Get the bounding box of the union geometry
    bounding_box = union_geometry.bounds
    return bounding_box

def get_dist_metres(line):
    new_dict = {}
    new_dict['geometry'] = ast.literal_eval(to_geojson(line))
    return calculate_distance(new_dict, Unit.meters)

def create_graph(osm_edges):

    #Node Tracking
    seen_coords = set()
    junctions = set()
    coords_to_edges = {}
    coords_ids = {}
    c_id = 0

    #Edge Tracking
    edges_dict = {}
    edge_id = 0

    count_its = 0


    for i,r in osm_edges.iterrows():
        #Get geom points
        #Get length of geom
        count_its += 1
        
        edge_coords = list(r['geometry'].coords)
        num_geoms = len(edge_coords)
        
        u = edge_coords[0]
        remaining_length = r['length']
        remaining_geom = r['geometry']
        edge_index = 1
        
        for coord in edge_coords:    
            dist_to_start = edge_index - 1
            dist_to_end = num_geoms - edge_index
            num_coords_in_geom = edge_coords.count(coord)
            
            #If coord is an existing junction in system
            if coord in junctions:
                
                #Check if junction along current edge
                if (dist_to_start > 0) & (dist_to_end > 0) & (num_coords_in_geom == 1):
                    #Split edge
                    split_geom = split(remaining_geom,Point(coord))
                    geom0_prop = get_dist_metres(split_geom.geoms[0]) / (get_dist_metres(split_geom.geoms[0]) + get_dist_metres(split_geom.geoms[1]))
                    geom1_prop = get_dist_metres(split_geom.geoms[1]) / (get_dist_metres(split_geom.geoms[0]) + get_dist_metres(split_geom.geoms[1]))
                    length_geom0 = geom0_prop * remaining_length
                    length_geom1 = geom1_prop * remaining_length
                
                    #Add first part of edge to edges
                    v = coord
                    edge_features = {}
                    edge_features['source'] = coords_ids[u]
                    edge_features['destination'] = coords_ids[v]
                    edge_features['length'] = length_geom0
                    edge_features['oneway'] = r['oneway']
                    edge_features['lts'] = r['lts']
                    edge_features['geometry'] = split_geom.geoms[0]
                    edges_dict[edge_id] = edge_features
                    edge_id = max(edges_dict.keys()) + 1
                    
                    #Update tracking matrices for remaining part of edge
                    u = v
                    remaining_length = length_geom1
                    remaining_geom = split_geom.geoms[1]
                
            #First time seeing a repeated coord which is not in self same edge - class this as a junction  
            elif (coord in seen_coords) & (num_coords_in_geom == 1):
                
                #Add node to junctions
                junctions.add(coord)
                
                #Check if intersection
                if (dist_to_start > 0) & (dist_to_end > 0):
                    #Split edge
                    split_geom = split(remaining_geom,Point(coord))
                    geom0_prop = get_dist_metres(split_geom.geoms[0]) / (get_dist_metres(split_geom.geoms[0]) + get_dist_metres(split_geom.geoms[1]))
                    geom1_prop = get_dist_metres(split_geom.geoms[1]) / (get_dist_metres(split_geom.geoms[0]) + get_dist_metres(split_geom.geoms[1]))
                    length_geom0 = geom0_prop * remaining_length
                    length_geom1 = geom1_prop * remaining_length
                        
                    v = coord
                    edge_features = {}
                    edge_features['source'] = coords_ids[u]
                    edge_features['destination'] = coords_ids[v]
                    edge_features['length'] = length_geom0
                    edge_features['oneway'] = r['oneway']
                    edge_features['lts'] = r['lts']
                    edge_features['geometry'] = split_geom.geoms[0]
                    edges_dict[edge_id] = edge_features
                    edge_id = max(edges_dict.keys()) + 1
                    
                    #Update tracking matrices for remaining part of edge
                    u = v
                    remaining_length = length_geom1
                    remaining_geom = split_geom.geoms[1]
                
                
                #Retrieve original edge and split if junction
                original_edge_id = coords_to_edges[coord]
                original_edge_dict = edges_dict[original_edge_id]
                geom_to_split_coords = list(original_edge_dict['geometry'].coords)
                coord_from_start = geom_to_split_coords.index(coord)
                coord_from_end = len(geom_to_split_coords) - geom_to_split_coords.index(coord) - 1
                num_coords_in_geom_temp = geom_to_split_coords.count(coord)
            
                if (coord_from_start > 0) & (coord_from_end > 0) & (num_coords_in_geom_temp == 1):
                    
                    split_geom = split(original_edge_dict['geometry'],Point(coord))
                    geom0_prop = get_dist_metres(split_geom.geoms[0]) / (get_dist_metres(split_geom.geoms[0]) + get_dist_metres(split_geom.geoms[1]))
                    geom1_prop = get_dist_metres(split_geom.geoms[1]) / (get_dist_metres(split_geom.geoms[0]) + get_dist_metres(split_geom.geoms[1]))
                    length_geom0 = geom0_prop * original_edge_dict['length']
                    length_geom1 = geom1_prop * original_edge_dict['length']
                    
                    temp_id1 = max(edges_dict.keys()) + 2
                    temp_id2 = max(edges_dict.keys()) + 3
                    
                    #Add edge 1
                    u_temp = geom_to_split_coords[0]
                    v_temp = coord                    
                    edge_features = {}
                    edge_features['source'] = coords_ids[u_temp]
                    edge_features['destination'] = coords_ids[v_temp]
                    edge_features['length'] = length_geom0
                    edge_features['oneway'] = r['oneway']
                    edge_features['lts'] = r['lts']
                    edge_features['geometry'] = split_geom.geoms[0]
                    edges_dict[temp_id1] = edge_features
                    
                    #Update coord references in geom
                    for c in list(split_geom.geoms[0].coords):
                        coords_to_edges[c] = temp_id1
                    
                    
                    
                    #Add edge 2
                    u_temp = coord
                    v_temp = geom_to_split_coords[-1]
                    edge_features = {}
                    edge_features['source'] = coords_ids[u_temp]
                    edge_features['destination'] = coords_ids[v_temp]
                    edge_features['length'] = length_geom1
                    edge_features['oneway'] = r['oneway']
                    edge_features['lts'] = r['lts']
                    edge_features['geometry'] = split_geom.geoms[1]
                    edges_dict[temp_id2] = edge_features

                    #Update coord references in geom
                    
                    for c in list(split_geom.geoms[1].coords):
                        coords_to_edges[c] = temp_id2
                    
                    #Delete original edge
                    del edges_dict[original_edge_id]

            else:
                seen_coords.add(coord)
                coords_to_edges[coord] = edge_id
                coords_ids[coord] = c_id
                c_id += 1
            
            edge_index += 1
        
        #Add remaining edge   
        v = coord
        edge_features = {}
        edge_features['source'] = coords_ids[u]
        edge_features['destination'] = coords_ids[v]
        edge_features['length'] = remaining_length
        edge_features['oneway'] = r['oneway']
        edge_features['lts'] = r['lts']
        edge_features['geometry'] = remaining_geom
        edges_dict[edge_id] = edge_features
        edge_id = max(edges_dict.keys()) + 1
        
    # Construct network

    all_edges = pd.DataFrame(edges_dict).T

    # Initialize an empty directed graph
    G = nx.DiGraph()

    # Iterate over the rows of the DataFrame and add edges to the graph
    for index, row in all_edges.iterrows():
        G.add_node(row['source'], pos=(row['geometry'].coords[0]))
        G.add_node(row['destination'], pos=(row['geometry'].coords[-1]))
        G.add_edge(row['source'], row['destination'], length=row['length'], lts = row['lts'])
        # If the edge is not one-way, add the reverse edge too
        if row['oneway'] == 0:
            G.add_edge(row['destination'], row['source'], length=row['length'], lts = row['lts'])
    

    # Extract nodes with their coordinates
    nodes_data = {
        'node': [],
        'latitude': [],
        'longitude': []
    }

    for node, data in G.nodes(data=True):
        nodes_data['node'].append(node)
        nodes_data['latitude'].append(data['pos'][0])
        nodes_data['longitude'].append(data['pos'][1])

    
    return G,pd.DataFrame(nodes_data),all_edges


def get_edges_from_geom(geom, con, project, table):
    
    # Define your query with a placeholder for the polygon WKT
    query = """
    SELECT r.*, ST_AsText(r.geometry) as geom_str
    FROM import.{} as r
    WHERE ST_Intersects(r.geometry, ST_GeomFromText(%s, 3857));
    """.format(table)
    
    # Execute the query
    with con.cursor() as cursor:
        cursor.execute(query, (transform(project.transform, geom).wkt,))
        results = cursor.fetchall()
        # Get column names
        colnames = [desc[0] for desc in cursor.description]
        
    # Convert results to a DataFrame
    expanded_edges = pd.DataFrame(results, columns=colnames)
    
    # Convert the geometry column to Shapely geometries
    expanded_edges['geometry'] = expanded_edges['geom_str'].apply(wkt.loads)
    
    # Return as a GeoDataFrame
    return gpd.GeoDataFrame(expanded_edges, geometry='geometry').set_crs(3857).to_crs(4326)

def sample_dataframe(df, min_samples ,rate=0.1):
    """
    Samples rows from the input DataFrame based on a variable sampling rate.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame.
    rate (float): The sampling rate for DataFrames with more than 10 rows.
    
    Returns:
    pd.DataFrame: The sampled DataFrame.
    """
    num_rows = len(df)
    
    if num_rows <= min_samples:
        # If the number of rows is 10 or less, return the entire DataFrame
        return df
    else:
        # Calculate the sample size
        sample_size = max(min_samples, int(num_rows * rate))
        
        # Sample the DataFrame
        sampled_df = df.sample(n=sample_size, random_state=1)
        return sampled_df

def mean_of_list(lst):
    return sum(lst)/len(lst)

def coarsen_network(od_bbox,base_graph,base_graph_nodes,oas,oa_centroids):
    od_polygon = gpd.GeoSeries([od_bbox])
    oas_in_area = oa_centroids[oa_centroids.intersects(od_polygon.unary_union)]
    
    # # Network nodes as points
    network_node_points_list = []
    for i,r in base_graph_nodes.iterrows():
        network_node_points_list.append(Point([r['latitude'],r['longitude']]))

    base_graph_nodes['geometry'] = network_node_points_list
    base_graph_nodes = gpd.GeoDataFrame(base_graph_nodes, geometry = 'geometry')
    
    area_to_nodes = {}
    num_nodes_per_oa = []
    for oa in list(oas_in_area.index):
        oa_geom = gpd.GeoSeries([oas.loc[oa]['geometry']])
        oa_nodes = base_graph_nodes[base_graph_nodes.intersects(oa_geom.unary_union)]
        area_to_nodes[oa] = oa_nodes
        num_nodes_per_oa.append(len(oa_nodes))
        
    #LTS1
    lts_1_edges = [(u, v, d) for u, v, d in base_graph.edges(data=True) if d.get('lts', 0) <= 1]
    G_lts1 = nx.Graph()
    G_lts1.add_nodes_from(list(base_graph.nodes))
    G_lts1.add_edges_from(lts_1_edges)

    #LTS2
    lts_2_edges = [(u, v, d) for u, v, d in base_graph.edges(data=True) if d.get('lts', 0) <= 2]
    G_lts2 = nx.Graph()
    G_lts2.add_nodes_from(list(base_graph.nodes))
    G_lts2.add_edges_from(lts_2_edges)

    #LTS3
    lts_3_edges = [(u, v, d) for u, v, d in base_graph.edges(data=True) if d.get('lts', 0) <= 3]
    G_lts3 = nx.Graph()
    G_lts3.add_nodes_from(list(base_graph.nodes))
    G_lts3.add_edges_from(lts_3_edges)

    #LTS4
    lts_4_edges = [(u, v, d) for u, v, d in base_graph.edges(data=True) if d.get('lts', 0) <= 4]
    G_lts4 = nx.Graph()
    G_lts4.add_nodes_from(list(base_graph.nodes))
    G_lts4.add_edges_from(lts_4_edges)
    
    oas_in_area_list = list(oas_in_area.index)
    
    paths_lts_1 = {}
    paths_lts_2 = {}
    paths_lts_3 = {}
    paths_lts_4 = {}

    for o in oas_in_area_list:
        paths_lts_1[o] = {}
        paths_lts_2[o] = {}
        paths_lts_3[o] = {}
        paths_lts_4[o] = {}
        for d in oas_in_area_list:
            paths_lts_1[o][d] = []
            paths_lts_2[o][d] = []
            paths_lts_3[o][d] = []
            paths_lts_4[o][d] = []

    # Define the original layers as MultiGraphs

    for o in oas_in_area_list:
        origin_nodes = sample_dataframe(area_to_nodes[o],4,0.1)
        for next_o_node in list(origin_nodes['node']):

            all_sp_lts1 = nx.single_source_dijkstra_path_length(G_lts1,next_o_node, weight = 'length')
            all_sp_lts2 = nx.single_source_dijkstra_path_length(G_lts2,next_o_node, weight = 'length')
            all_sp_lts3 = nx.single_source_dijkstra_path_length(G_lts3,next_o_node, weight = 'length')
            all_sp_lts4 = nx.single_source_dijkstra_path_length(G_lts4,next_o_node, weight = 'length')

            for d in oas_in_area_list:
                destination_nodes = list(area_to_nodes[d]['node'])

                #Add LTS 1 Routes to Dict
                paths_to_dest = {key: all_sp_lts1[key] for key in destination_nodes if key in all_sp_lts1}
                #Remove 0 elements
                paths_to_dest = [element for element in list(paths_to_dest.values()) if element != 0]
                #Add found routes to dictionary
                paths_lts_1[o][d].extend(paths_to_dest)

                #Add LTS 2 Routes to Dict
                paths_to_dest = {key: all_sp_lts2[key] for key in destination_nodes if key in all_sp_lts2}
                #Remove 0 elements
                paths_to_dest = [element for element in list(paths_to_dest.values()) if element != 0]
                #Add found routes to dictionary
                paths_lts_2[o][d].extend(paths_to_dest)
                
                #Add LTS 3 Routes to Dict
                paths_to_dest = {key: all_sp_lts3[key] for key in destination_nodes if key in all_sp_lts3}
                #Remove 0 elements
                paths_to_dest = [element for element in list(paths_to_dest.values()) if element != 0]
                #Add found routes to dictionary
                paths_lts_3[o][d].extend(paths_to_dest)
                
                #Add LTS 4 Routes to Dict
                paths_to_dest = {key: all_sp_lts4[key] for key in destination_nodes if key in all_sp_lts4}
                #Remove 0 elements
                paths_to_dest = [element for element in list(paths_to_dest.values()) if element != 0]
                #Add found routes to dictionary
                paths_lts_4[o][d].extend(paths_to_dest)
                
            
    G_c = nx.DiGraph()
    G_c.add_nodes_from(oas_in_area_list)

    for o in oas_in_area_list:
        for d in oas_in_area_list:
            distances = []
            if len(paths_lts_1[o][d]) > 0:
                distances.append(mean_of_list(paths_lts_1[o][d]))
            else:
                distances.append(0)

            if len(paths_lts_2[o][d]) > 0:
                distances.append(mean_of_list(paths_lts_2[o][d]))
            else:
                distances.append(0)

            if len(paths_lts_3[o][d]) > 0:
                distances.append(mean_of_list(paths_lts_3[o][d]))
            else:
                distances.append(0)
                
            if len(paths_lts_4[o][d]) > 0:
                distances.append(mean_of_list(paths_lts_4[o][d]))
            else:
                distances.append(0)
            
            if sum(distances) > 0:
                G_c.add_edge(o, d, distance=distances)

    return G_c, oas_in_area_list


def assign_weights(J,types):
    #------------------------------------------#
    #input: network with 'distance' attribute and type of risk function to apply
    # output: netowrk with 'distance' and 'risk' attribute
    # assigns risk (or discomfort) attribute to each edge of the network
    # risk is a quadruple of values, just like distance, one relative to each 'street-layer' considered 
    #------------------------------------------#
    if types=='Linear':
        for edge in J.edges:
            beta = [None]*4
            beta[0] = J[edge[0]][edge[1]]['distance'][0]
            beta[1] = J[edge[0]][edge[1]]['distance'][1]*0.66
            beta[2] = J[edge[0]][edge[1]]['distance'][2]*0.33
            beta[3] = 0
            J[edge[0]][edge[1]]['risk'] = beta
    return J

def updateParetoFront(new_length, new_risk, best_length_risk_states):
    #------------------------------------------#
    # Input:  dist [int], risk[int], paretoFront [sorted dictionary]
    #output:[bolean] if the new_length and new_dist have been added to the ParetoFront returns True, False otherwise 
    # adds new dist, and new risk if the path is pareto optimal (non dominated) 
    #------------------------------------------#
    
    # Find before and after points on length-axis of relation
    after_idx = best_length_risk_states.bisect_right( new_length );
    before_idx = after_idx - 1;
    after_risk = -1
    before_risk = -1 
    
    # Find associated risks
    keys = best_length_risk_states.keys()
    if after_idx < len(best_length_risk_states):
        after_risk = best_length_risk_states.get(keys[after_idx])
    if before_idx >= 0:
        before_risk = best_length_risk_states.get(keys[before_idx])
       
    # Exclude new path if the shorter path has same or less risk
    if before_risk <= new_risk and not before_risk == -1:
        return False

    # Exclude longer or equal paths with higher risks:
    final_exclude_idx = after_idx
    while(final_exclude_idx < len(best_length_risk_states)):
        risk = best_length_risk_states.get(keys[final_exclude_idx])
        if risk < new_risk:
            break
        final_exclude_idx = final_exclude_idx + 1
    del best_length_risk_states.keys()[after_idx:final_exclude_idx]
    best_length_risk_states[new_length] = new_risk
    
    return True

def compute_pareto_fronts(J, o):
    #------------------------------------------#
    #input:  aggregated network [networkx graph], node [int]
    #output: dictionary of all pareto fronts from node 'o' to all other nodes in the network [dictionary] 
    #this implements the multi-objective optimization criteria on all possible OD pairs from node 'o' to all 
    #other nodes of the network (thus it computes ParetoFronts based on distance and risk criteria)
    #------------------------------------------#

    states_to_explore = []
    best_length_risk_states_for_node = {o : SortedDict({0: 0})}
    
    states_to_explore = []
    cnt = 0
    heapq.heappush(states_to_explore, (0, 0, cnt, { 'distance': 0, 'risk': 0, 'prevstate': None, 'node': o }))
    
    while len(states_to_explore) > 0:
        current_distance, current_risk, dummy, current_state = heapq.heappop(states_to_explore)
        current_node = current_state.get('node')
        
        # Check if the explored path is still part of the optimal length risk curve.
        # Otherwise, discard.        
        best_length_risk_states = best_length_risk_states_for_node.get(current_node)
        if best_length_risk_states.get(current_distance, -1) != current_risk:
            continue
        
        for tt in J.out_edges(current_node):
            next_node = tt[1]
            edge = J[current_state.get('node')][next_node]
            
            best_length_risk_states = best_length_risk_states_for_node.get(next_node, None)
            if best_length_risk_states is None:
                best_length_risk_states = SortedDict()
                best_length_risk_states_for_node[next_node] = best_length_risk_states
            
            # Loop over different network choices for link
            for network_i in range(len(edge['distance'])):
                
                # Compute new length
                edge_length = edge['distance'][network_i]
                if edge_length == 0:
                    continue
                new_length = current_distance + edge_length
                
                # Compute new risk
                edge_risk = edge['risk'][network_i]
                new_risk = current_risk + edge_risk
                
                inserted = updateParetoFront(new_length, new_risk, best_length_risk_states)
                if not inserted:
                    continue
                
                new_state = { 'distance': new_length, 'risk': new_risk, 'prevstate': current_state, 'node': next_node}
                # Remove removed states from states_to_explore
                
                heapq.heappush(states_to_explore, (new_length, new_risk, cnt, new_state))
                cnt = cnt + 1
                
    return best_length_risk_states_for_node

def dict_to_array(dict):
    # Convert dictionary to list of lists (each key-value pair becomes a list)
    data_list = [[key, value] for key, value in dict.items()]

    # Convert list of lists to a 2D NumPy array
    return np.array(data_list)

def generate_reference_point(pareto_front, margin=1.0):
    max_values = np.max(pareto_front, axis=0)
    reference_point = max_values + margin
    return reference_point

def hypervolume(pareto_front, reference_point):
    pareto_sorted = pareto_front[np.argsort(pareto_front[:, 0])]
    hv = 0.0
    for i in range(len(pareto_sorted)):
        if i == 0:
            width = reference_point[0] - pareto_sorted[i, 0]
        else:
            width = pareto_sorted[i-1, 0] - pareto_sorted[i, 0]
        height = reference_point[1] - pareto_sorted[i, 1]
        hv += width * height
    return hv


# Spread/diversity calculation
def spread_diversity(pareto_front):
    distances = np.sqrt(np.sum(np.diff(pareto_front, axis=0)**2, axis=1))
    diversity = np.std(distances)
    return diversity


# Convexity (signed area method)
def compute_signed_area(pareto_front):
    hull = ConvexHull(pareto_front)
    hull_area = hull.volume
    pareto_area = np.trapz(pareto_front[:, 1], pareto_front[:, 0])
    signed_area = hull_area - pareto_area
    return signed_area, hull


# Function to compute the trade-off rate
def compute_trade_off_rate(pareto_front):
    trade_off_rates = []
    for i in range(len(pareto_front) - 1):
        delta_x = pareto_front[i+1, 0] - pareto_front[i, 0]
        delta_y = pareto_front[i+1, 1] - pareto_front[i, 1]
        trade_off_rate = delta_y / delta_x
        trade_off_rates.append(trade_off_rate)
    return np.mean(trade_off_rates), trade_off_rates


def haversine_distance(point1, point2):
    """
    Calculate the Haversine distance between two shapely Point objects.
    
    :param point1: shapely Point object for the first point
    :param point2: shapely Point object for the second point
    :return: distance in kilometers
    """
    
    # Extract latitude and longitude from Point objects
    latlon1 = (point1.y, point1.x)  # Point.y is latitude, Point.x is longitude
    latlon2 = (point2.y, point2.x)  # Point.y is latitude, Point.x is longitude
    
    # Haversine formula
    lat1, lon1 = map(math.radians, latlon1)
    lat2, lon2 = map(math.radians, latlon2)
    
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.asin(math.sqrt(a))
    
    # Radius of earth in kilometers
    r = 6371
    
    distance = c * r
    
    return distance