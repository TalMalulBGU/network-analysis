import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import json
import ast
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import ast
import os
import glob
import itertools
from collections import Counter, deque
from decimal import *
from functools import partial


def read_relevancy(path):
	df = pd.read_excel(path)
	df['relevance'] = df['relevance'].map(lambda r: str(r) == 'yes' or str(r).lower() == 'true')
	df['Color'] = df['Color'].map(lambda c: tuple(map(float, [y for y in [x.replace('(','').replace(')','') for x in str(c).split(',')] if len(y) > 0])))
	return df
	
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
	
def community_taxonomy_plot(ax, clusters, graph, title, modularity_score, taxonomy_level, taxonomy_relevancy):

	cluster_organism_dictionary = {}

	for cluster_id, cluster in clusters.items():

		cluster_vertices = graph.vs[cluster]
		if len(cluster_vertices) > 0:
			entry_relevant_toxonomy = cluster_vertices.select(lambda v: any(t in v[taxonomy_level] for t in list(taxonomy_relevancy.index)))
			organism_counter = Counter(t for l in entry_relevant_toxonomy[taxonomy_level] for t in l if t in list(taxonomy_relevancy.index))
			cluster_organism_dictionary[cluster_id] = {'Counter': organism_counter, 'Total_Nodes': len(cluster_vertices), 'Function_Nodes': len(entry_relevant_toxonomy)}



	taxonomies = list(dict([taxonomy for cluster in cluster_organism_dictionary.values() for taxonomy in cluster['Counter'].most_common()[:5]]).keys())
	clusters_taxonomies = []
	buttom = list(np.zeros(len(cluster_organism_dictionary.keys())))
	d_taxonomy_colors = dict(zip(taxonomy_relevancy.index, taxonomy_relevancy.loc[taxonomy_relevancy.index, 'Color']))
	d_taxonomy_colors = dict(sorted(d_taxonomy_colors.items(), key=lambda item: item[1]))

	for idx, taxonomy in enumerate(taxonomies):
		clusters_taxomomy = []
		for cluster_id in cluster_organism_dictionary.keys():
			clusters_taxomomy.append(cluster_organism_dictionary[cluster_id]['Counter'].get(taxonomy, 0))
		clusters_taxonomies.append(clusters_taxomomy)

		ax.bar(list(map(str, cluster_organism_dictionary.keys())), clusters_taxomomy, bottom=buttom, color=d_taxonomy_colors.get(taxonomy,'black'), label=taxonomy)
		buttom = [buttom[i] + v for i, v in enumerate(clusters_taxomomy)]


	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.set_title(title)
	ax.set_xlabel('Cluster ID')
	ax.set_ylabel('Ocurrences')
	ax.set_xticklabels(['{}\n{}/{}'.format(cid, d['Function_Nodes'], d['Total_Nodes']) for cid, d in cluster_organism_dictionary.items()])
	ax.text(1.3, 1.0, "Modularity = {:.2f}".format(modularity_score), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
	ax.legend(ncol=1, bbox_to_anchor=(1.1, 1.0), fancybox=True)


def _main():

	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--alignment_weight', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--algorithm_file', help='name of edge attribute', type=str, required=True)
	parser.add_argument('--taxonomy_level', help='the taxonomy level', type=str, required=True)
	parser.add_argument('--relevance_taxonomy_xlsx', help='the taxonomy relevancy group', required=True, type=str)
	parser.add_argument('--image_title', help='the name of the protien', type=str, required=True)
	parser.add_argument('--min_in_cluster', help='name of edge attribute', type=int, required=False, default=0)
	parser.add_argument('--out_file', help='output file location', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)



	options = parser.parse_args()


	   
	print("protein = " + options.protein)
	print("xgmml_file = " + options.xgmml_file)
	print("alignment_weight = " + str(options.alignment_weight))
	print("algorithm_file = " + str(options.algorithm_file))
	print("min_in_cluster = " + str(options.min_in_cluster))
	print("image_title = {}".format(options.image_title))
	print("taxonomy_level = {}".format(options.taxonomy_level))
	print("relevance_taxonomy_xlsx = {}".format(options.relevance_taxonomy_xlsx))
	print("out_file = {}".format(options.out_file))

	algorithm = None
	uniprot_db = None
	
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)

	graph = load_network(options.xgmml_file, uniprot_db)
		
	
	with open(options.algorithm_file) as json_file:
		data = json.load(json_file)
		graph = graph.subgraph(graph.vs.select(name_in = [name for lst_names in data.values() for name in lst_names]))
		membership = np.zeros(graph.vcount(), dtype=int)
		for cluster_id in data.keys():
			ids = [v.index for v in graph.vs.select(name_in = data[cluster_id])]
			membership[ids] = cluster_id
		algorithm = ig.clustering.VertexClustering(graph, membership, modularity_params={'weights': options.alignment_weight})



	clusters = {cluster_id: cluster for cluster_id, cluster in enumerate(algorithm) if len(cluster) >= options.min_in_cluster}
	df_relevancy = read_relevancy(options.relevance_taxonomy_xlsx)
	df_relevancy = df_relevancy.set_index(options.taxonomy_level)

	fig, ax = plt.subplots(1, 1, figsize=(int(40), 8))
	community_taxonomy_plot(ax, clusters, graph, options.image_title, algorithm.modularity, options.taxonomy_level, df_relevancy)
	plt.show()	
	plt.savefig(options.out_file, transparent=True, dpi=300)

if __name__ == '__main__':
	_main()


