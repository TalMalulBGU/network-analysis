import sys
import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import json
import ast
import os
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

	uniprot_db = None
	if location.endswith('.csv'):
		uniprot_db = pd.read_csv(location)
	else:
		uniprot_db = pd.read_excel(location, engine="openpyxl")

	if len(uniprot_db) > 0:
		uniprot_db['Entry'] = uniprot_db['Entry'].astype(str)
		uniprot_db['Sequence'] = uniprot_db['Sequence'].astype(str)
		

	uniprot_db = uniprot_db.set_index('Entry')
	return uniprot_db


def load_network(location):
	f = open(location, 'rb')
	nx_g = networkxgmml.XGMMLReader(f)
	f.close()
	graph = ig.Graph.from_networkx(nx_g)
	graph.vs['name'] = graph.vs['_nx_name']
	nx_g.clear()

	return graph


def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('protein', help='the name of the protien', type=str)
	parser.add_argument('jsonfile', help='input json file location', type=str)
	parser.add_argument('database', help='database file location for filtering', type=str)
	parser.add_argument('weights', help='name of edge attribute', type=str, default=None)


	options, unknown_args = parser.parse_known_args()
	

	my_root_protein_location = '/home/talmalu/thesis/data/' + options.protein + '/'
	my_root_output_location = my_root_protein_location + 'output/'
	normalized_weight_column_name = 'normalized_{}'.format(options.weights)
	
	print("database = " + options.database)
	print("protein = " + options.protein)
	print("jsonfile = " + options.jsonfile)
	print("weights = " + str(options.weights))

	
	uniprot_db = None
	if options.database:
		uniprot_db = load_db(options.database)
		

	df_alignment_scores = None
	with open(options.jsonfile) as json_file:
		data = json.load(json_file)
		df_alignment_scores = pd.DataFrame(data)
		df_alignment_scores.rename(columns = {'soruce':'source'}, inplace = True)
		df_alignment_scores.set_index(['source', 'target'], inplace = True)
		
		
	if df_alignment_scores is None:
		sys.exit()
		
	custom_db = uniprot_db.loc[list(set([entry for entries in df_alignment_scores.index for entry in entries]))]
	if 'Length' not in custom_db.columns:
		custom_db['Length'] = custom_db['Sequence'].map(len)
	   
	custom_db = custom_db[['Length']]
	custom_db['Length'] = custom_db['Length'].astype(np.int16)
		
	for e_length in list(custom_db['Length'].unique()):
		entries = list(custom_db[custom_db['Length'] == e_length].index)
		df_alignment_scores.loc[df_alignment_scores.index.get_level_values(0).isin(entries), 'source_length'] = e_length
		df_alignment_scores.loc[df_alignment_scores.index.get_level_values(1).isin(entries), 'target_length'] = e_length
		
	df_alignment_scores['source_length'] = df_alignment_scores['source_length'].astype(np.int16)
	df_alignment_scores['target_length'] = df_alignment_scores['target_length'].astype(np.int16)
	df_alignment_scores[options.weights] = df_alignment_scores[options.weights].astype(np.int16)
	
	df_alignment_scores[normalized_weight_column_name] = np.vectorize(np.int16)(
	-(np.vectorize(Decimal.log10)
	(
		(Decimal(2) ** (Decimal(-1) * df_alignment_scores[options.weights])) *
		df_alignment_scores['source_length'] * 
		df_alignment_scores['target_length']
	)))
	
	df_alignment_scores[options.weights] = df_alignment_scores[options.weights].astype(float)
	df_alignment_scores[normalized_weight_column_name] = df_alignment_scores[normalized_weight_column_name].astype(float)
	df_alignment_scores = df_alignment_scores.reset_index()
	
	df_alignment_scores.loc[:, ~df_alignment_scores.columns.isin(['source_length', 'target_length'])].to_json('{}.normalized'.format(options.jsonfile), orient='records')

	

if __name__ == '__main__':
	_main()


