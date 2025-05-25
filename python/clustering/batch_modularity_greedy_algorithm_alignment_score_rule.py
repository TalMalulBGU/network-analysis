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


def read_relevancy(path):
	df = pd.read_excel(path)
	df['relevance'] = df['relevance'].map(lambda r: r == 'yes')
	df['Color'] = df['Color'].map(lambda c: tuple(map(float, str(c).split(','))))
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
	
def get_edges_le(graph, alignment_weight, alignment_score):
	return graph.es.select(lambda e: e[alignment_weight] <= alignment_score)

def get_edges_ge(graph, alignment_weight, alignment_score):
	return graph.es.select(lambda e: e[alignment_weight] >= alignment_score)
	
def edges_function(rule):
	if rule == 'GE':
		return get_edges_ge
	elif rule == 'LE':
		return get_edges_le

def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--alignment_score', help='output file location', type=int, required=True)
	parser.add_argument('--alignment_weight', help='output file location', type=str, required=True)
	parser.add_argument('--alignment_score_rule', help='name of edge attribute', type=str, required=False, default='GE')
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)

	options = parser.parse_args()


	print("protein = " + options.protein)
	print("output = " + options.output)
	print("xgmml_file = {}".format(options.xgmml_file))
	print("alignment_score = {}".format(options.alignment_score))
	print("alignment_weight = {}".format(options.alignment_weight))
	print("alignment_score_rule = {}".format(options.alignment_score_rule))
	
	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)
		

	graph = load_network(options.xgmml_file, uniprot_db)
	
	es_selection = edges_function(options.alignment_score_rule)

	modularities = {}
	n_vertices = {}
	n_edges = {}
	n_clusters = {}
	
	ig_graph = graph.subgraph_edges(es_selection(graph, options.alignment_weight, options.alignment_score))
	cluster_community_dendogram = ig_graph.community_fastgreedy(weights=options.alignment_weight)
	cluster_community_vc = cluster_community_dendogram.as_clustering()
	modularity = cluster_community_vc.modularity
	modularities[options.alignment_score] = modularity
	n_vertices[options.alignment_score] = ig_graph.vcount()
	n_edges[options.alignment_score] = ig_graph.ecount()
	n_clusters[options.alignment_score] = len(cluster_community_vc)


	result_dictionary = {}

	result_dictionary[options.alignment_score] = {'modularity': modularities[options.alignment_score],
							'n_vertices': n_vertices[options.alignment_score],
							'n_edges': n_edges[options.alignment_score], 'n_clusters': n_clusters[options.alignment_score], 'rule': options.alignment_score_rule
						   }

	json_object = json.dumps(result_dictionary, indent = 4)

	with open(options.output, "w") as outfile:
		outfile.write(json_object)
	
	
	
	
if __name__ == '__main__':
	_main()


