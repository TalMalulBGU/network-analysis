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
	df['relevance'] = df['relevance'].map(lambda r: str(r) == 'yes' or str(r).lower() == 'true')
	df['Color'] = df['Color'].map(lambda c: tuple(map(float, [y for y in [x.replace('(','').replace(')','') for x in str(c).split(',')] if len(y) > 0])))
	return df
	
def get_graph_function_modularity(graph, weights, db, relevancy = None):
	
	f_col_name = 'Gene ontology (molecular function)'
	entry_db = db.loc[graph.vs['name']]
	entry_relevant_function_db = entry_db.loc[entry_db[f_col_name].map(lambda lst: any(f in lst for f in list(relevancy.index)))]

	modularity = {}
	lst_vertices_functions = {entry: lst_function[0] for entry, lst_function in entry_relevant_function_db[f_col_name].items()}
	set_vertices_functions = set(lst_vertices_functions.values())
	cluster_function = dict(zip(range(len(set_vertices_functions)),set_vertices_functions))
	function_cluster = dict(zip(set_vertices_functions, range(len(set_vertices_functions))))
	init_vertex_labels = [function_cluster[f] for f in lst_vertices_functions.values()]
	label_propagation_vertex_clustering = graph.community_label_propagation(weights=weights, initial=init_vertex_labels,fixed = [True] * len(init_vertex_labels))
	modularity = label_propagation_vertex_clustering.modularity

	return modularity

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
	parser.add_argument('--weight', help='output file location', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--relevance_function_xlsx', help='the taxonomy relevancy group', type=str, required=True)
	parser.add_argument('--reverse', required=False, action='store_true')
	parser.add_argument('--no-reverse', required=False, dest='reverse', action='store_false')
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)

	options = parser.parse_args()

	print("protein = {}".format(options.protein))
	print("graph = {}".format(options.xgmml_file))
	print("weight = {}".format(options.weight))
	print("reverse = {}".format(options.reverse))
	print("output = {}".format(options.output))
	print("relevance_function_xlsx = {}".format(options.relevance_function_xlsx))


	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)

	graph = load_network(options.xgmml_file, uniprot_db)


	alignment_modularity = {}
	n_alignment_scores_edges = {}


	min_alignment_score = int(min(graph.es[options.weight]))
	max_alignment_score = int(max(graph.es[options.weight]))
	alignment_step = 0
	edge_filter_function = None
	alignment_score = min_alignment_score
	n_max_edges = graph.ecount()

	
	if options.reverse:
		alignment_start = max_alignment_score
		alignment_end = min_alignment_score - 1
		alignment_step = -1
		edge_filter_function = lambda e, alignment_score: e[options.weight] <= alignment_score
		ax.invert_xaxis()

	else:
		alignment_start = min_alignment_score
		alignment_end = max_alignment_score + 1
		alignment_step = 1
		edge_filter_function = lambda e, alignment_score: e[options.weight] >= alignment_score


	f_col_name = 'Gene ontology (molecular function)'
	function_relevancy = read_relevancy(options.relevance_function_xlsx)
	function_relevancy = function_relevancy[function_relevancy['relevance']]
	function_relevancy = function_relevancy.set_index(f_col_name)

	entry_db = uniprot_db.loc[graph.vs['name']]
	entry_relevant_function_db = entry_db.loc[entry_db[f_col_name].map(lambda lst: any(f in lst for f in list(function_relevancy.index)))]
	graph = graph.subgraph(entry_relevant_function_db.index)
	
	for alignment_score in range(alignment_start, alignment_end, alignment_step):
		edges = graph.es.select(lambda e: edge_filter_function(e, alignment_score))
		graph = graph.subgraph_edges(edges, delete_vertices=False)
		alignment_modularity[alignment_score] = get_graph_function_modularity(graph, options.weight, uniprot_db, function_relevancy)
		n_alignment_scores_edges[alignment_score] = len(edges)


	result_dictionary = {}
	for i in range(alignment_start, alignment_end, alignment_step):
		result_dictionary[i] = {'modularity': alignment_modularity[i], 'n_edges': n_alignment_scores_edges[i]}

	json_object = json.dumps(result_dictionary, indent = 4)

	with open(options.output, "w") as outfile:
		outfile.write(json_object)



if __name__ == '__main__':
	_main()


