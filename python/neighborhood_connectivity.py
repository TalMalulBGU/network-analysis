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
	parser.add_argument('--xgmml_file', help='xgmml_file file location', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--alignment_weight', help='output file location', type=str, required=False, default=None)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)

	options, unknown_args = parser.parse_known_args()

	
	print("protein = {}".format(options.protein))
	print("output = {}".format(options.output))
	print("xgmml_file = {}".format(options.xgmml_file))
	print("alignment_weight = {}".format(options.alignment_weight))
	
	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)
	
	graph = load_network(options.xgmml_file, uniprot_db)
	
	graph_degrees = None
	weighted_graph_degrees = None
	
	if options.alignment_weight:
		weighted_graph_degrees = []
		for vertex in graph.vs:
			vertex_edges = vertex.incident(mode='all')
			alignment_scores = [edge[options.alignment_weight] for edge in vertex_edges]
			weighted_graph_degrees.append(sum(alignment_scores))


	graph_degrees = graph.vs.degree()
	
	
	
	result_dictionary = {}
	for vertex in graph.vs:
		vertex_neighbors = graph.neighbors(vertex)
		vertex_degree = len(vertex_neighbors)
		result_dictionary[vertex['name']] = {}
		result_dictionary[vertex['name']]['degree'] = vertex_degree
		result_dictionary[vertex['name']]['neighbors_degree'] = []
		
		if weighted_graph_degrees:
			result_dictionary[vertex['name']]['weighted_degree'] = weighted_graph_degrees[vertex.index]
			result_dictionary[vertex['name']]['neighbors_weighted_degree'] = []
			
		for neighbor in vertex_neighbors:
			result_dictionary[vertex['name']]['neighbors_degree'].append(graph_degrees[neighbor])
			if weighted_graph_degrees:
				result_dictionary[vertex['name']]['neighbors_weighted_degree'].append(weighted_graph_degrees[neighbor])


	result_dictionary['database'] = options.database
	result_dictionary['xgmml_file'] = options.xgmml_file
	result_dictionary['alignment_weight'] = options.alignment_weight
	
	json_object = json.dumps(result_dictionary, indent = 4)



	with open(options.output, "w") as outfile:
		outfile.write(json_object)

	
	
if __name__ == '__main__':
	_main()


