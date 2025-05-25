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

def create_clustering_size_distribution(ax, clusters, title):
	communities_sizes_counter = Counter(clusters.sizes())
	cuc = dict(sorted(communities_sizes_counter.items()))
	
	ax.bar(list(map(str,cuc.keys())), cuc.values())
	ax.set_title(title)
	ax.set_xlabel('Size')
	ax.set_ylabel('Ocurrences')
	ax.text(0.9, 0.9, "Mean Size = {:.2f}".format(np.mean(list(cuc.values()))), horizontalalignment='center',
			verticalalignment='center', transform=ax.transAxes)

	for tick in ax.get_xticklabels():
		tick.set_rotation(70)
		

def _main():

	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--out_file', help='output folder location', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--alignment_weight', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--algorithm_file', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--image_title', help='the name of the protien', type=str, required=False, default='Distribution of Communities Sizes')

	options = parser.parse_args()


	print("protein = " + options.protein)
	print("xgmml_file = " + options.xgmml_file)
	print("out_file = " + options.out_file)
	print("alignment_weight = " + str(options.alignment_weight))
	print("algorithm_file = " + str(options.algorithm_file))
	print("image_title = {}".format(options.image_title))
	print("database = " + options.database)


	graph = load_network(options.xgmml_file)
	graph.to_undirected(mode='each')
	
	algorithm = None
	uniprot_db = None
	
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)

	graph = load_network(options.xgmml_file, uniprot_db)
	
	
	
	with open(options.algorithm_file) as json_file:
		data = json.load(json_file)
		graph = graph.subgraph(graph.vs.select(name_in = [name for lst_names in data.values() for name in lst_names]))
		membership = np.zeros(graph.vcount(), dtype=int)
		for cluster_id in data.keys():
			ids = [v.index for v in graph.vs.select(name_in = data[cluster_id])]
			membership[ids] = cluster_id
		algorithm = ig.clustering.VertexClustering(graph, membership, modularity_params={'weights': options.alignment_weight})
	
	
	
	fig, ax = plt.subplots(1, 1, figsize=(20,6))

	create_clustering_size_distribution(ax, algorithm, options.image_title)

	plt.savefig(options.out_file, bbox_inches='tight', transparent=True)
	plt.show()

			

if __name__ == '__main__':
	_main()


