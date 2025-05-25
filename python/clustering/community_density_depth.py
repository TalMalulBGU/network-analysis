import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import json
import ast
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import ast
import os
import glob
import itertools
from collections import Counter, deque
from decimal import *
from functools import partial

def load_vertexclustering(algorithm_file, graph, alignment_weight):
	
	algorithm = None
	with open(algorithm_file) as json_file:
		data = json.load(json_file)
		graph = graph.subgraph(graph.vs.select(name_in = [name for lst_names in data.values() for name in lst_names]))
		membership = np.zeros(graph.vcount(), dtype=int)
		for cluster_id in data.keys():
			ids = [v.index for v in graph.vs.select(name_in = data[cluster_id])]
			membership[ids] = cluster_id
		algorithm = ig.clustering.VertexClustering(graph, membership, modularity_params={'weights': alignment_weight})
		
	return algorithm

def read_relevancy(path):
	df = pd.read_excel(path)
	df['relevance'] = df['relevance'].map(lambda r: str(r) == 'yes' or str(r).lower() == 'true')
	df['Color'] = df['Color'].map(lambda c: tuple(map(float, [y for y in [x.replace('(','').replace(')','') for x in str(c).split(',')] if len(y) > 0])))
	return df
	
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
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--data_folder', help='data folder location', type=str, required=True)
	parser.add_argument('--alignment_weight', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--algorithm_file', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--min_in_cluster', help='name of edge attribute', type=int, required=False, default=0)
	parser.add_argument('--out_file', help='output file location', type=str, required=True)
	parser.add_argument('--graph_name', help='output file location', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)


	options = parser.parse_args()


	   
	print("protein = " + options.protein)
	print("xgmml_file = " + options.xgmml_file)
	print("alignment_weight = " + str(options.alignment_weight))
	print("algorithm_file = " + str(options.algorithm_file))
	print("data_folder = " + str(options.data_folder))
	print("min_in_cluster = " + str(options.min_in_cluster))
	print("out_file = {}".format(options.out_file))
	print("graph_name = {}".format(options.graph_name))
	

	algorithm = None
	uniprot_db = None
	
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)

	graph = load_network(options.xgmml_file, uniprot_db)
		
	
	clustering_dictionary_queue = deque()
	algorithm = None

	algorithm = load_vertexclustering(options.algorithm_file, graph, options.alignment_weight)
	algorithms_densities = []
	clustering_dictionary_queue.append((algorithm, options.graph_name, 0, algorithms_densities))

	while(len(clustering_dictionary_queue) > 0):
		algorithm, file_name, cluster_id, parent = clustering_dictionary_queue.pop()
		
		me = {'cluster_id': cluster_id, 'file_name': file_name, 'density': algorithm.graph.density(),
			  'n_vertices': algorithm.graph.vcount(), 'n_edges': algorithm.graph.ecount(), 'childs': [] }
		parent.append(me)
		
		clusters = {idx: cluster for idx, cluster in enumerate(algorithm)}

		for cluster_id, cluster in clusters.items():
			cluster_graph = algorithm.subgraph(cluster_id)
			if cluster_graph.vcount() < options.min_in_cluster:
				continue

			algorithm_file_path = '{}/{}_{}.json'.format(options.data_folder, file_name, cluster_id)
			if os.path.exists(algorithm_file_path):
				cluster_community_vc = load_vertexclustering(algorithm_file_path, cluster_graph, options.alignment_weight)
				clustering_dictionary_queue.append((cluster_community_vc, '{}_{}'.format(file_name, cluster_id), cluster_id, me['childs']))
			
	json_object = json.dumps(algorithms_densities, indent = 4)

	with open(options.out_file, "w") as outfile:
		outfile.write(json_object)

if __name__ == '__main__':
	_main()


