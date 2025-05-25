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
from collections import Counter, OrderedDict
from functools import partial


def read_relevancy(path):
	df = pd.read_excel(path)
	df['relevance'] = df['relevance'].map(lambda r: r == 'yes')
	df['Color'] = df['Color'].map(lambda c: tuple(map(float, str(c).split(','))))
	return df


def create_colored_taxomony_level_degrees_histogram(graph, t_level, ax, title, colors, min_size_threshold = 0):


	taxomonies_in_level = set(taxomony for taxomonies in graph.vs[t_level] for taxomony in taxomonies)
	taxomony_seperated = {}
	ordred_taxomony_seperated = None
	hist_colors = []
	degrees = []
	lables = []

	for taxomony in taxomonies_in_level:
		taxomony_degrees = graph.vs.select(lambda v: taxomony in v[t_level]).degree()
		if len(taxomony_degrees) > min_size_threshold:
			taxomony_seperated[taxomony] = (taxomony_degrees, colors.get(taxomony,('black')))
	
	if len(taxomony_seperated) > 0:
		ordred_taxomony_seperated = OrderedDict(sorted(taxomony_seperated.items()))
		lables = list(ordred_taxomony_seperated.keys())
		degrees, hist_colors = list(map(list, zip(*list(ordred_taxomony_seperated.values()))))
		


	(taxomonies_degree_count ,taxomonies_degree_bins, taxomonies_degree_pathches) = ax.hist(degrees,
																				bins=15, stacked=True,
																				color=hist_colors,
																				label=lables)


	ax.yaxis.set_major_formatter(plticker.PercentFormatter(xmax=(sum(taxomonies_degree_count[-1]))))
	ax.yaxis.set_major_locator(plticker.LinearLocator(10))

	# for b in taxomonies_degree_pathches:
	#	 ax.bar_label(b)

	ax.set_xticks(taxomonies_degree_bins)
	ax.set_title(title)
	ax.set_xlabel('Degrees')
	ax.set_ylabel('Ocurrences')
	ax.legend(ncol=3, bbox_to_anchor=(1.0, 1.0), framealpha=0.5)


	return (taxomonies_degree_count ,taxomonies_degree_bins, taxomonies_degree_pathches)


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
	parser.add_argument('--min_size_threshold', help='minimum number of taxonomy nodes in the graph', type=int, default=0)
	parser.add_argument('--relevance_taxonomy_xlsx', help='the taxonomy relevancy group', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)



	options, unknown_args = parser.parse_known_args()


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
	
	fig, ax = plt.subplots(figsize=(20,8))
	(local_bacteria_phylum_degree_count ,local_bacteria_phylum_degree_bins, local_bacteria_phylum_degree_pathches) = create_colored_taxomony_level_degrees_histogram(graph=graph,
																																				   t_level=options.taxonomy_level,
																																				   ax=ax,
																																				   title='Histogram of Nodes Degrees in {}'.format(options.taxonomy_level),
																																				   colors=dict(zip(df_relevancy.index, df_relevancy['Color'])),
																																				   min_size_threshold=options.min_size_threshold
																																				   )


	plt.savefig(options.output, bbox_inches='tight', transparent=True)

	plt.show()
	
if __name__ == '__main__':
	_main()


