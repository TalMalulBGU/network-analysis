import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import json
import ast
import os
import glob
import itertools
from decimal import *

def array_converter(value):
	result = []

	if pd.notnull(value) and value != '':
		if ';' in value:
			result = str(value).strip().split(';')[:-1]
		else:
			result = ast.literal_eval(value)

	return list(result)

def load_db(location):

	uniprot_db = pd.read_excel(location, engine="openpyxl")

	if len(uniprot_db) > 0:
		uniprot_db['Entry'] = uniprot_db['Entry'].astype(str)
#		uniprot_db['Status'] = uniprot_db['Status'].map(lambda v:  v == 'reviewed')
#		uniprot_db['Protein names'] = uniprot_db['Protein names'].astype(str)
#		uniprot_db['Organism'] = uniprot_db['Organism'].astype(str)
#		uniprot_db['Organism ID'] = uniprot_db['Organism ID'].astype(int)
		uniprot_db['Sequence'] = uniprot_db['Sequence'].astype(str)
#		uniprot_db['Gene names'] = uniprot_db['Gene names'].map(partial(array_converter, sep=' ', cast_type=str))
#		uniprot_db['Cross-reference (Pfam)'] = uniprot_db['Cross-reference (Pfam)'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['Cross-reference (InterPro)'] = uniprot_db['Cross-reference (InterPro)'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['EC number'] = uniprot_db['EC number'].map(str)
#		uniprot_db['Function [CC]'] = uniprot_db['Function [CC]'].map(str)
#		uniprot_db['Rhea ID'] = uniprot_db['Rhea ID'].map(str)
#		uniprot_db['Keywords'] = uniprot_db['Keywords'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['Keyword ID'] = uniprot_db['Keyword ID'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['Gene ontology IDs'] = uniprot_db['Gene ontology IDs'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['Gene ontology (molecular function)'] = uniprot_db['Gene ontology (molecular function)'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['Gene ontology (GO)'] = uniprot_db['Gene ontology (GO)'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['Gene ontology (cellular component)'] = uniprot_db['Gene ontology (cellular component)'].map(partial(array_converter, sep=';', cast_type=str))
#		uniprot_db['Gene ontology (biological process)'] = uniprot_db['Gene ontology (biological process)'].map(partial(array_converter, sep=';', cast_type=str))
		if 'Fragment' in uniprot_db.columns:
			uniprot_db['Fragment'] = uniprot_db['Fragment'].map(lambda v:  v == 'fragment' or v == 'fragments')
		uniprot_db['Fragment'] = uniprot_db['Fragment'].map(lambda v:  v == 'fragment' or v == 'fragments')

	uniprot_db = uniprot_db.rename(columns={'Status': 'Reviewed'})
	uniprot_db = uniprot_db.set_index('Entry')
	
	return uniprot_db


def load_network(location, db = None, filter_db = False):
	f = open(location, 'rb')
	nx_g = networkxgmml.XGMMLReader(f)
	f.close()
	graph = ig.Graph.from_networkx(nx_g)
	nx_g.clear()
	
	graph.to_undirected(mode='each')
	graph.vs['name'] = graph.vs['_nx_name']
	
	
	if db is not None and filter_db:
		db_missing_entries = set(graph.vs['name']) - set(db.index)
		graph.delete_vertices(graph.vs.select(lambda v: v['name'] in db_missing_entries))
   		if 'Fragment' in db.columns:
			graph_fragments_entries = set(graph.vs['name']).intersection(set(db[(db['Fragment'])].index))
			graph.delete_vertices(graph.vs.select(lambda v: v['name'] in graph_fragments_entries))

	return graph


def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)
	parser.add_argument('--alignment_weight', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--alignment_length', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--json_files', help='name of edge attribute', type=str, required=True, nargs='*', dest='files')


	options = parser.parse_args()
	
	uniprot_db = None


	print("database = " + options.database)
	print("protein = " + options.protein)
	print("output = " + options.output)
	print("alignment_weight = " + str(options.alignment_weight))
	print("alignment_length = " + str(options.alignment_length))
	print("json_files = {}".format(options.files))
	print("n_json_files = {}".format(len(options.files)))



	if options.database:
		uniprot_db = load_db(options.database)
	
	
	alignment_scores_dictionaries = []
	for file in sorted(options.files):
		with open(file) as json_file:
			data = json.load(json_file)
			alignment_scores_dictionaries = itertools.chain(alignment_scores_dictionaries, data)

	df_alignment_scores = pd.DataFrame(alignment_scores_dictionaries)
	df_alignment_scores.rename(columns = {'soruce':'source'}, inplace = True)
	df_alignment_scores.set_index(['source', 'target'], inplace = True)

	
	
	
	alignment_scores = []
	alignment_lengths = []

	df_edges_weights = df_alignment_scores.reset_index()[[options.alignment_weight, options.alignment_length]]
	df_edges_weights = df_edges_weights.set_index(options.alignment_weight)
	
	
	for alignment_score in sorted(set(df_edges_weights.index)):
		selected_lengths = df_edges_weights.loc[alignment_score, options.alignment_length].tolist()
		listed_lengths = selected_lengths
		if type(selected_lengths) is not list:
			listed_lengths = [selected_lengths]
		alignment_lengths.append(list(map(lambda f: float(round(Decimal(f), 2)), listed_lengths)))
		alignment_scores.append(alignment_score)
	
	result_dictionary = {}
	for index, alignment_score in enumerate(alignment_scores):
		result_dictionary[int(alignment_score)] = alignment_lengths[index]

	json_object = json.dumps(result_dictionary, indent = 4)

	with open(options.output, "w") as outfile:
		outfile.write(json_object)

if __name__ == '__main__':
	_main()


