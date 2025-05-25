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
		uniprot_db['Sequence'] = uniprot_db['Sequence'].astype(str)
		uniprot_db['Cross-reference (Pfam)'] = uniprot_db['Cross-reference (Pfam)'].map(array_converter)
		uniprot_db['Cross-reference (InterPro)'] = uniprot_db['Cross-reference (InterPro)'].map(array_converter)
		uniprot_db['Fragment'] = uniprot_db['Fragment'].map(lambda v:  v == 'fragment')

	uniprot_db = uniprot_db.set_index('Entry')
	return uniprot_db


def load_network(location):
	f = open(location, 'rb')
	nx_g = networkxgmml.XGMMLReader(f)
	f.close()

	return nx_g


def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--alignment_weight', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--alignment_threshold', help='name of edge attribute', type=int, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--remove_attributes', help='output file location', required=False, nargs='*', default=[])
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--json_files', help='name of edge attribute', type=str, required=True, nargs='*', dest='files')


	options, unknown_args = parser.parse_known_args()
	

	print("protein = " + options.protein)
	print("output = " + options.output)
	print("alignment_weight = " + str(options.alignment_weight))
	print("alignment_threshold = " + str(options.alignment_threshold))
	print("xgmml_file = {}".format(options.xgmml_file))
	print("files = {}".format(options.files))

	graph = load_network(options.xgmml_file)
	graph_empty_copy = nx.create_empty_copy(graph).to_undirected()
	
	alignment_scores_dictionaries = []

	
	for file in sorted(options.files):
		with open(file) as json_file:
			data = json.load(json_file)
			for edge_data in data:			
				edge_source = edge_data['source']
				edge_target = edge_data['target']
				del edge_data['source']
				del edge_data['target']			   
				
				for att in options.remove_attributes:
					if att in edge_data and att != options.alignment_weight:
						del edge_data[att]

				if edge_data[options.alignment_weight] >= options.alignment_threshold:
					existing_edge = graph_empty_copy.get_edge_data(edge_source, edge_target, default=None)
					if existing_edge and existing_edge[options.alignment_weight] >= edge_data[options.alignment_weight]:
						graph_empty_copy.remove_edge(edge_source, edge_target)

					graph_empty_copy.add_edge(edge_source,
											 edge_target,
											 **edge_data)
	
	f = open(options.output, 'w')
	networkxgmml.XGMMLWriter(f ,graph_empty_copy, '',directed=False)
	f.close()

if __name__ == '__main__':
	_main()


