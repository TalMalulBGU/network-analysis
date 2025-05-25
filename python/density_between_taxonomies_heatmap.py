import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.patches import Rectangle
import seaborn as sns
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
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--taxonomy_level', help='the taxonomy level', type=str, required=True)
	parser.add_argument('--edges_count_xlsx', help='number of taxonomy to calculate', type=str)
	parser.add_argument('--edges_possible_count_xlsx', help='number of taxonomy to calculate', type=str)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)


	options = parser.parse_args()


	print("protein = " + options.protein)
	print("taxonomy_level = {}".format(options.taxonomy_level))
	print("output = " + options.output)
	print("xgmml_file = {}".format(options.xgmml_file))
	print("edges_count_xlsx = {}".format(options.edges_count_xlsx))
	print("edges_possible_count_xlsx = {}".format(options.edges_possible_count_xlsx))

	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)

	graph = load_network(options.xgmml_file, uniprot_db)
	taxonomy_links = pd.read_excel(options.edges_count_xlsx, index_col=0).astype(float)
	possible_taxonomy_links = pd.read_excel(options.edges_possible_count_xlsx, index_col=0).astype(float)	
	
	df = taxonomy_links / possible_taxonomy_links
	taxomony_in_level_counter = Counter(taxonomy for lst_taxonomy in graph.vs[options.taxonomy_level] for taxonomy in lst_taxonomy)

	xy_label_ticks = []
	df_mean = df.mean(skipna=True).mean(skipna=True)

	for column_taxonomy in df.columns:
		vs_size = taxomony_in_level_counter[column_taxonomy]
		xy_label_ticks.append('{}\n({})'.format(column_taxonomy, vs_size))
		

	fig, ax = plt.subplots(1, 1, figsize=(20,6))
	sns.heatmap(df ,fmt=".5f",annot=True, cmap='RdYlGn',linewidths=0.30,ax=ax, xticklabels=xy_label_ticks, yticklabels=xy_label_ticks, vmin=0, vmax=1)

	for column_taxonomy in df.columns:
		vs_size = taxomony_in_level_counter[column_taxonomy]
		xy_label_ticks.append('{}\n({})'.format(column_taxonomy, vs_size))
		lst_above_mean = df.index[(df >= df_mean)[column_taxonomy]].to_list()
		if len(lst_above_mean) > 0:
			col = df.columns.get_loc(column_taxonomy)
			for index_taxonomy in lst_above_mean:
				row = df.index.get_loc(index_taxonomy)
				ax.add_patch(Rectangle((col, row), 1, 1, fill=False, edgecolor='blue', lw=5))
				


	ax.text(0.9, 0.9, "Mean Size = {:.5f}".format(df_mean), horizontalalignment='center',
			verticalalignment='center', transform=ax.transAxes)
	ax.set_title('Density between taxonomies - {}'.format(options.taxonomy_level))

	plt.savefig(options.output, bbox_inches='tight', transparent=True)

	plt.show()
	
if __name__ == '__main__':
	_main()


