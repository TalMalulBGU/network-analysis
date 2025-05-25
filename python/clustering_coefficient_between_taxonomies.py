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


def get_cluster_coefficient_between_taxonomy(graph, find, t_level, min_nodes = 0):
   
	results = dict()
	
	taxomony_in_level_counter = Counter(taxonomy for lst_taxonomy in graph.vs[t_level] for taxonomy in lst_taxonomy)
	taxonomies = list(taxomony_in_level_counter.keys())
	for s_index in range(len(taxonomies)):
		s_taxonomy = taxonomies[s_index]
		for t_index in range(s_index, len(taxonomies)):
			t_taxonomy = taxonomies[t_index]
			if s_taxonomy in find and t_taxonomy in find and taxomony_in_level_counter[s_taxonomy] >= min_nodes and taxomony_in_level_counter[t_taxonomy] >= min_nodes:
				if s_taxonomy not in results:
					results[s_taxonomy] = dict()
				if t_taxonomy not in results[s_taxonomy]:
					taxonomy_subgraph = graph.subgraph(graph.vs.select(lambda v: s_taxonomy in v[t_level] or t_taxonomy in v[t_level]))
					results[s_taxonomy][t_taxonomy] = taxonomy_subgraph.transitivity_undirected()
					 
	
	reformed_dict = {}
	for outerKey, innerDict in results.items():
		for innerKey, values in innerDict.items():
			multi_key = tuple(sorted([outerKey, innerKey]))
			if multi_key not in reformed_dict:
				reformed_dict[multi_key] = 0
			reformed_dict[multi_key] = values

	multiIndex_df = pd.DataFrame(reformed_dict.values(), index=reformed_dict.keys(), columns=['Count'])
	pivoted_df = multiIndex_df.reset_index().rename(columns={'level_0': 'Node_1', 'level_1': 'Node_2'}).pivot(
		columns='Node_1',index='Node_2',values='Count')			
	return pivoted_df
	
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
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--taxonomy_level', help='the taxonomy level', type=str, required=True)
	parser.add_argument('--min_size_threshold', help='number of taxonomy to calculate', type=int, default=0)
	parser.add_argument('--relevance_taxonomy_xlsx', help='the taxonomy relevancy group', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)

	options = parser.parse_args()


	print("protein = " + options.protein)
	print("output = " + options.output)
	print("xgmml_file = {}".format(options.xgmml_file))
	print("taxonomy_level = {}".format(options.taxonomy_level))
	print("min_size_threshold = {}".format(options.min_size_threshold))
	print("relevance_taxonomy_xlsx = {}".format(options.relevance_taxonomy_xlsx))
	
	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)
		

	graph = load_network(options.xgmml_file, uniprot_db)
	df_relevancy = read_relevancy(options.relevance_taxonomy_xlsx)
	df_relevancy = df_relevancy.set_index(options.taxonomy_level)
	
	graph_taxonomy_links = get_cluster_coefficient_between_taxonomy(graph, df_relevancy.index, options.taxonomy_level, options.min_size_threshold)
	graph_taxonomy_links.to_excel(options.output, index=True)
	
	
	
if __name__ == '__main__':
	_main()


