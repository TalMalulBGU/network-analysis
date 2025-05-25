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
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--json_files', help='name of edge attribute', type=str, required=True, nargs='*', dest='files')


	options, unknown_args = parser.parse_known_args()
	

	print("output = {}".format(options.output))
	print("files = {}".format(options.files))

	
	alignment_scores_dictionaries = []

	
	result_dictionary = {}
	for file in sorted(options.files):
		with open(file) as json_file:
			result_dictionary.update(json.load(json_file))
	
	json_object = json.dumps(result_dictionary, indent = 4)

	with open(options.output, "w") as outfile:
			outfile.write(json_object)


if __name__ == '__main__':
	_main()


