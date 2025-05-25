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
import importlib.util
import sys


def create_module(file_path, module_name):
	spec = importlib.util.spec_from_file_location(module_name, file_path)
	new_module = importlib.util.module_from_spec(spec)
	sys.modules[module_name] = new_module
	spec.loader.exec_module(new_module)

def create_modules_from_directory(directory):
	stack = []
	for foldername, subfolders, filenames in os.walk(directory):
			for filename in filenames:
				if filename.endswith(".py"):
					file_path = os.path.join(foldername, filename)
					module_name = f"{foldername.replace(f'{os.path.dirname(directory)}/', '')}.{os.path.splitext(filename)[0]}"
					module_name = module_name.replace(os.path.sep, '.')
					stack.append((file_path, module_name))

	while (len(stack) > 0):
		file_path, module_name = stack.pop(0)
		try:
			create_module(file_path, module_name)
		except (ModuleNotFoundError, ImportError):
			print(f"Module {module_name} not found. Adding to retry stack...")
			stack.append((file_path, module_name))


directory_path = "/home/talmalu/thesis/projects/python/tutils"
create_modules_from_directory(directory_path)

from tutils.databases import *


def array_converter(value, sep=';', cast_type = str):
	result = []
	
	if pd.notnull(value) and value != '':
		if sep in value:
			result = [ v.lstrip().strip() for v in str(value).split(sep)[:-1] ]
		else:
			result.append(cast_type(value))
			
	return list(result)

def load_db(location):

	uniprot_db = UniprotDB.read(location)

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
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--alignment_weight', help='output file location', type=str, required=False)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)


	options, unknown_args = parser.parse_known_args()

	print("protein = {}".format(options.protein))
	print("xgmml_file = {}".format(options.xgmml_file))
	print("output = {}".format(options.output))

	result_dictionary = {}
	result_dictionary['protein'] = options.protein
	result_dictionary['xgmml_file'] = options.xgmml_file
	
	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)
		result_dictionary['database'] = options.database
		
	graph = load_network(options.xgmml_file, uniprot_db)
	
	
	result_dictionary['n_vertices'] = graph.vcount()
	result_dictionary['n_edges'] = graph.ecount()
	
	degrees = graph.vs.degree()
	result_dictionary['min_degree'] = min(degrees)
	result_dictionary['max_degree'] = max(degrees)
	
	if options.alignment_weight:
		print("alignment_weight = {}".format(options.alignment_weight))
		result_dictionary['alignment_weight'] = options.alignment_weight
		result_dictionary['min_alignment_weight'] = min(graph.es[options.alignment_weight])
		result_dictionary['max_alignment_weight'] = max(graph.es[options.alignment_weight])

	
	
	json_object = json.dumps(result_dictionary, indent = 4)

	with open(options.output, "w") as outfile:
		outfile.write(json_object)

if __name__ == '__main__':
	_main()


