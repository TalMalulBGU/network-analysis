import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import json
import ast
import os
import glob
import itertools
from decimal import *
from collections import Counter
from functools import partial

def get_vertex_cluster_connections(graph, algorithm, weighted = False, weights = 'alignment_score'):

	# Get the number of clusters
	num_clusters = len(algorithm)

	# Initialize an empty table
	table = []
	memory = {}
	
	vertex_cluster_lookup = {v: c for c, cluster in enumerate(algorithm) for v in cluster}

	for vertex in graph.vs:

		# Get the connections to other clusters
		
		connections = []
		connected_clusters = None
		if weighted:
			connected_clusters = [(vertex_cluster_lookup[n], graph.es[graph.get_eid(vertex, n, directed=False)][weights]) for n in graph.neighbors(vertex) if n != vertex]
		else:
			connected_clusters = [(vertex_cluster_lookup[n], 1) for n in graph.neighbors(vertex) if n != vertex]
		connections.extend((cluster[0], cluster[1]) for cluster in connected_clusters)
			
		consolidated_connections = {}
		for cluster, count in connections:
			consolidated_connections[cluster] = consolidated_connections.get(cluster, 0) + count

		# Append cluster and connections to the table
		table.append((vertex.index, vertex, consolidated_connections))

	return table



def vertex_z_score(graph, algorithm):

	z = np.zeros(graph.vcount())

	clusters_mean_degree = []
	clusters_std_degree = []

	for cluster_id, cluster in enumerate(algorithm):
		cluster_node_degrees = graph.subgraph(cluster).vs.degree()
		# cluster_node_degrees = graph.vs[cluster].degree()
		cluster_mean_degree = np.mean(cluster_node_degrees)
		cluster_std_degree = np.std(cluster_node_degrees)
		clusters_mean_degree.append(cluster_mean_degree)
		clusters_std_degree.append(cluster_std_degree)


		z[cluster] = (np.array(cluster_node_degrees) - cluster_mean_degree) / cluster_std_degree
		
	return z

def vertex_participation_coefficient(graph, algorithm, weighted=False, weights='alignment_score'):

	p = np.zeros(graph.vcount())
	vertex_cluster = get_vertex_cluster_connections(graph, algorithm, weighted, weights)

	for vertex_index, cluster, connections in vertex_cluster:
		if weighted:
			node_weights = sum(graph.es[graph.incident(vertex_index)][weights])
			p[vertex_index] = 1 - np.sum(np.power(np.array(list(connections.values())) / node_weights, 2))
		else:		
			node_degree = graph.vs[vertex_index].degree()
			p[vertex_index] = 1 - np.sum(np.power(np.array(list(connections.values())) / node_degree, 2))
	
	return p
	
	
def create_functional_role_array(z_array, c_array, thresholds_z, thresholds_c, roles):
	num_roles = len(roles)
	
	conditions = []
	
	thresholds_z = [-np.inf] + thresholds_z + [np.inf]
	thresholds_c = [-np.inf] + thresholds_c + [np.inf]

	for i in range(len(thresholds_z) - 1):
		lower_threshold_z = thresholds_z[i]
		upper_threshold_z = thresholds_z[i+1]
		
		for j in range(len(thresholds_c) - 1):
			lower_threshold_c = thresholds_c[j]
			upper_threshold_c = thresholds_c[j+1]
			
			condition = (z_array > lower_threshold_z) & (z_array <= upper_threshold_z) & (c_array > lower_threshold_c) & (c_array <= upper_threshold_c)
			conditions.append(condition)
			
	roles_array = np.select(conditions, roles, default='???')
	
	return roles_array

def array_converter(value, sep=';', cast_type = str):
	result = []
	
	if pd.notnull(value) and value != '':
		if sep in value:
			result = [ v.lstrip().strip() for v in str(value).split(sep)[:-1] ]
		else:
			result.append(cast_type(value))
			
	return list(result)

def load_db(location):

	uniprot_db = pd.read_excel(location, engine="openpyxl")

	if len(uniprot_db) > 0:
		uniprot_db['Entry'] = uniprot_db['Entry'].astype(str)
		uniprot_db['Status'] = uniprot_db['Status'].map(lambda v:  v == 'reviewed')
		uniprot_db['Protein names'] = uniprot_db['Protein names'].astype(str)
		uniprot_db['Organism'] = uniprot_db['Organism'].astype(str)
		uniprot_db['Organism ID'] = uniprot_db['Organism ID'].astype(int)
		uniprot_db['Sequence'] = uniprot_db['Sequence'].astype(str)
		uniprot_db['Gene names'] = uniprot_db['Gene names'].map(partial(array_converter, sep=' ', cast_type=str))
		uniprot_db['Cross-reference (Pfam)'] = uniprot_db['Cross-reference (Pfam)'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Cross-reference (InterPro)'] = uniprot_db['Cross-reference (InterPro)'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['EC number'] = uniprot_db['EC number'].map(str)
		uniprot_db['Function [CC]'] = uniprot_db['Function [CC]'].map(str)
		uniprot_db['Rhea ID'] = uniprot_db['Rhea ID'].map(str)
		uniprot_db['Keywords'] = uniprot_db['Keywords'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Keyword ID'] = uniprot_db['Keyword ID'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Gene ontology IDs'] = uniprot_db['Gene ontology IDs'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Gene ontology (molecular function)'] = uniprot_db['Gene ontology (molecular function)'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Gene ontology (GO)'] = uniprot_db['Gene ontology (GO)'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Gene ontology (cellular component)'] = uniprot_db['Gene ontology (cellular component)'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Gene ontology (biological process)'] = uniprot_db['Gene ontology (biological process)'].map(partial(array_converter, sep=';', cast_type=str))
		uniprot_db['Fragment'] = uniprot_db['Fragment'].map(lambda v:  v == 'fragment' or v == 'fragments')

	uniprot_db = uniprot_db.rename(columns={'Status': 'Reviewed'})
	uniprot_db = uniprot_db.set_index('Entry')
	
	return uniprot_db


def load_network(location, db = None):
	f = open(location, 'rb')
	nx_g = networkxgmml.XGMMLReader(f)
	f.close()
	graph = ig.Graph.from_networkx(nx_g)
	nx_g.clear()
	
	graph.to_undirected(mode='each')
	graph.vs['name'] = graph.vs['_nx_name']
	
	
	if db is not None:
		db_missing_entries = set(graph.vs['name']) - set(db.index)
		graph.delete_vertices(graph.vs.select(lambda v: v['name'] in db_missing_entries))
		graph_fragments_entries = set(graph.vs['name']).intersection(set(db[(db['Fragment'])].index))
		graph.delete_vertices(graph.vs.select(lambda v: v['name'] in graph_fragments_entries))

	return graph

def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--image', help='image file location', type=str, required=True)
	parser.add_argument('--algorithm_file', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--c_thresholds', help="list of participation coefficient thresholds", type=float, nargs='+', required=True)
	parser.add_argument('--z_thresholds', help="list of whitin module degree thresholds", type=float, nargs='+', required=True)
	parser.add_argument('--roles', help="list of roles shuold be equivalent to (c_thresholds + 1) * (z_thresholds + 1)", type=str, nargs='+', required=True)
	options = parser.parse_args()


	print("protein = " + options.protein)
	print("image = " + options.image)
	print("roles = {}".format(options.roles))
	print("c_thresholds = {}".format(options.c_thresholds))
	print("z_thresholds = {}".format(options.z_thresholds))
	print("algorithm_file = " + str(options.algorithm_file))
	
	

	with open(options.algorithm_file) as json_file:
		data = json.load(json_file)
		
		xgmml_file = data['xgmml_file']
		algorithm_type = data['type']
		algorithm_clusters = data['clusters']
		
		database = data.get('database', None)
		names = data.get('names', None)
		alignment_score = data.get('alignment_score', None)
		alignment_weight = data.get('alignment_weight', None)
		alignment_score_rule = data.get('alignment_score_rule', None)
		
		uniprot_db = None
		if database:
			print("database = {}".format(database))
			uniprot_db = load_db(database)
			
		graph = load_network(xgmml_file, uniprot_db)
		
		if names:
			print("names = {}".format(names))
			db_cluster_vertices = graph.vs.select(name_in = names)
			if len(db_cluster_vertices) == 0:
				print("names {}, was not found in graph".format(names))
				return
			graph.subgraph(db_cluster_vertices)
			
		if all([item is not None for item in [alignment_score, alignment_weight, alignment_score_rule]]):
			print("alignment_score = {}".format(alignment_score))
			print("alignment_weight = {}".format(alignment_weight))
			print("alignment_score_rule = {}".format(alignment_score_rule))
			es_selection = edges_function(alignment_score_rule)
			graph = graph.subgraph_edges(es_selection(graph, alignment_weight, alignment_score))
		
		membership = np.zeros(graph.vcount(), dtype=int)
		for cluster_id in algorithm_clusters.keys():
			ids = [v.index for v in graph.vs.select(name_in = algorithm_clusters[cluster_id])]
			membership[ids] = cluster_id
		algorithm = ig.clustering.VertexClustering(graph, membership, modularity_params={'weights': alignment_weight})

	z = vertex_z_score(graph, algorithm)
	c = vertex_participation_coefficient(graph, algorithm)
	vertices_roles = create_functional_role_array(z, np.log(c), options.z_thresholds, options.c_thresholds, options.roles)
	
	fig, ax = plt.subplots(1, 1, figsize=(10,10))
	title = 'Protein roles scatter'

	for role in options.roles:
		z_role = z[(vertices_roles == role) & (~np.isnan(z))]
		c_role = c[(vertices_roles == role) & (~np.isnan(z))]
		
		ax.scatter(np.log(c_role), z_role, label=role, alpha=0.3)
		
	ax.set_title('{}'.format(title))
	ax.set_xlabel('Participant coefficient, c')
	ax.set_ylabel('Within-module degree, z')

	ax.legend()

	plt.savefig(options.image, bbox_inches='tight', transparent=True)

	plt.show()
	
if __name__ == '__main__':
	_main()


