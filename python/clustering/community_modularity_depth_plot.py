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
from distinctipy import distinctipy


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
	parser.add_argument('--out_file', help='output file location', type=str, required=True)
	parser.add_argument('--image_title', help='image title', type=str, required=True)
	parser.add_argument('--data_file', help='name of edge attribute', type=str, required=True)

	options = parser.parse_args()
	

	print("protein = " + options.protein)
	print("out_file = " + options.out_file)
	print("data_file = {}".format(options.data_file))
	print("image_title = {}".format(options.image_title))



	clustering_dictionary_queue = deque()
	branches = []

	with open(options.data_file) as json_file:
		data = json.load(json_file)
		algorithms_modularities = data


	clustering_dictionary_queue.append((algorithms_modularities[0], []))


	while(len(clustering_dictionary_queue) > 0):
		me, path = clustering_dictionary_queue.pop()
		if len(me['childs']) == 0:
			branches.append(path + [(me['cluster_id'], me['modularity'], me['n_vertices'], me['n_edges'])])
		else:
			for child in me['childs']:
				clustering_dictionary_queue.append((child, path + [(me['cluster_id'], me['modularity'],  me['n_vertices'], me['n_edges'])]))



	
	fig, ax = plt.subplots(1, 1, figsize=(20, 8))

	colors = distinctipy.get_colors(len(branches),pastel_factor=0.3)

	for index, branch in enumerate(branches):
		cluster_ids = []
		modularities = []
		n_vertices = []
		n_edges = []
		for cluster_id, modularity, vertices, edges in branch:
			cluster_ids.append(cluster_id)
			modularities.append(modularity)
			n_vertices.append(vertices)
			n_edges.append(edges)
		depths = range(len(branch))


		ax.plot(depths, modularities, label='_'.join(map(str,cluster_ids)), color=colors[index], alpha=0.5)
		
		for i, txt in enumerate(n_vertices):
			ax.annotate(txt, (depths[i], modularities[i]), alpha=0.5)


	ax.legend(ncol=3, bbox_to_anchor=(1.1, 1.0), fancybox=True)	
	ax.set_title(options.image_title)
	ax.set_xlabel('Depth')
	ax.set_ylabel('Modularity')
	
	plt.show()

	plt.savefig(options.out_file, transparent=True, dpi=300)




if __name__ == '__main__':
	_main()


