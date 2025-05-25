import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import json
import ast
import os
import glob
import itertools
from decimal import *
from collections import Counter
from functools import partial

def get_graph_taxonomy_modularity(graph, weights, relevancy = None, taxonomy_levels = ['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']):
	

	modularity = {}

	for taxonomy_level in taxonomy_levels:
		lst_vertices_taxonomies = [lst_taxonomy[0] for lst_taxonomy in graph.vs[taxonomy_level]]
		set_vertices_taxonomies = set(lst_vertices_taxonomies)
		cluster_taxonomy = dict(zip(range(len(set_vertices_taxonomies)),set_vertices_taxonomies))
		taxonomy_cluster = dict(zip(set_vertices_taxonomies, range(len(set_vertices_taxonomies))))
		init_vertex_labels = [taxonomy_cluster[t] for t in lst_vertices_taxonomies]
		label_propagation_vertex_clustering = graph.community_label_propagation(weights=weights, initial=init_vertex_labels,fixed = [True] * len(init_vertex_labels))
		modularity[taxonomy_level] = label_propagation_vertex_clustering.modularity
		
	return modularity

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

def read_relevancy(path):
	df = pd.read_excel(path)
	df['relevance'] = df['relevance'].map(lambda r: r == 'yes')
	df['Color'] = df['Color'].map(lambda c: tuple(map(float, str(c).split(','))))
	return df


def get_edges_le(graph, alignment_weight, alignment_score):
	return graph.es.select(lambda e: e[alignment_weight] <= alignment_score)

def get_edges_ge(graph, alignment_weight, alignment_score):
	return graph.es.select(lambda e: e[alignment_weight] >= alignment_score)
	
def edges_function(rule):
	if rule == 'GE':
		return get_edges_ge
	elif rule == 'LE':
		return get_edges_le

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
	parser.add_argument('--algorithm_file', help='output file location', type=str, required=True)
	parser.add_argument('--taxonomy_level', help='output file location', type=str, required=True)
	parser.add_argument('--relevance_taxonomy_xlsx', help='output file location', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--min_in_cluster', help='output file location', type=int, required=False, default=0)
	
	

	options = parser.parse_args()

	print("protein = {}".format(options.protein))
	print("algorithm_file = {}".format(options.algorithm_file))
	print("output = {}".format(options.output))
	print("min_in_cluster = {}".format(options.min_in_cluster))

	with open(options.algorithm_file) as json_file:
		data = json.load(json_file)
		
		xgmml_file = data['xgmml_file']
		algorithm_type = data['type']
		algorithm_clusters = data['clusters']
		
		database = data.get('database', None)
		names = data.get('names', None)
		alignment_score = data.get('alignment_score', None)
		alignment_weight = data.get('alignment_weight', None)
		alignment_score_rule = data.get('alignment_score_rule', None)
		
		uniprot_db = None
		if database:
			print("database = {}".format(database))
			uniprot_db = load_db(database)
			
		graph = load_network(xgmml_file, uniprot_db)
		
		if names:
			print("names = {}".format(names))
			db_cluster_vertices = graph.vs.select(name_in = names)
			if len(db_cluster_vertices) == 0:
				print("names {}, was not found in graph".format(names))
				return
			graph.subgraph(db_cluster_vertices)
			
		if all([item is not None for item in [alignment_score, alignment_weight, alignment_score_rule]]):
			print("alignment_score = {}".format(alignment_score))
			print("alignment_weight = {}".format(alignment_weight))
			print("alignment_score_rule = {}".format(alignment_score_rule))
			es_selection = edges_function(alignment_score_rule)
			graph = graph.subgraph_edges(es_selection(graph, alignment_weight, alignment_score))
		
		membership = np.zeros(graph.vcount(), dtype=int)
		for cluster_id in algorithm_clusters.keys():
			ids = [v.index for v in graph.vs.select(name_in = algorithm_clusters[cluster_id])]
			membership[ids] = cluster_id
		algorithm = ig.clustering.VertexClustering(graph, membership, modularity_params={'weights': alignment_weight})



	clusters = {cluster_id: cluster for cluster_id, cluster in enumerate(algorithm) if len(cluster) >= options.min_in_cluster}
	df_relevancy = read_relevancy(options.relevance_taxonomy_xlsx)
	df_relevancy = df_relevancy.set_index(options.taxonomy_level)
	taxonomy_relevancy = df_relevancy
	attrs = {}
	colors = {}
	taxonomies = {}
	for cluster_id, cluster in clusters.items():
		cluster_vertices = graph.vs[cluster]
		entry_relevant_toxonomy = cluster_vertices.select(lambda v: any(t in v[options.taxonomy_level] for t in list(taxonomy_relevancy.index)))
		organism_counter = Counter(t for l in entry_relevant_toxonomy[options.taxonomy_level] for t in l if t in list(taxonomy_relevancy.index))
		# cluster_organism_dictionary[cluster_id] = {'Counter': organism_counter, 'Total_Nodes': len(cluster_vertices), 'Function_Nodes': len(entry_relevant_toxonomy)}
		# taxonomies = list(dict([taxonomy for cluster in cluster_organism_dictionary.values() for taxonomy in cluster['Counter'].most_common()[:5]]))
		d_taxonomy_colors = dict(zip(taxonomy_relevancy.index, taxonomy_relevancy.loc[taxonomy_relevancy.index, 'Color']))
		d_taxonomy_colors = dict(sorted(d_taxonomy_colors.items(), key=lambda item: item[1]))
		attrs[cluster_id] = []
		colors[cluster_id] = []
		taxonomies[cluster_id] = []
		for taxonomy, count in organism_counter.most_common()[:5]:
			attrs[cluster_id].append(count)
			colors[cluster_id].append(d_taxonomy_colors.get(taxonomy, (0, 0, 0)))
			taxonomies[cluster_id].append(taxonomy)
			
	# make graph

	G = algorithm.cluster_graph().to_networkx()
	# pos = nx.layout.kamada_kawai_layout(G, scale=10)
	# pos = {node: np.array(np.round(cord, 7) for node, cord in pos.items()}
	pos = nx.nx_agraph.graphviz_layout(G)

	legend_latxonomy_set = set()


	fig, ax = plt.subplots(1, 1, figsize=(25, 25))
	# nx.draw_networkx_nodes(G, pos=pos, ax=ax)
	for node in G.nodes:
		# if len(list(G.neighbors(node))) > 0:
		# if sum(attrs[node]) > 20:
		legend_latxonomy_set.update(taxonomies[node])
		patches, texts = ax.pie(
			np.array(attrs[node]), # s.t. all wedges have equal size
			# [1] * len(attrs[node]), # s.t. all wedges have equal size
			center=pos[node], 
			colors = colors[node],
			# labels = taxonomies[node], 
			labeldistance = None,
			radius = 25
		)
		[p.set_zorder(2) for p in patches]
		ax.text(pos[node][0], pos[node][1], str(sum(np.array(attrs[node]))), size=6, ha='left', va='baseline')
	edges = nx.draw_networkx_edges(G, pos=pos, ax=ax)
	# edges.set_zorder(1)

	legend_elements = []
	for taxonomy in sorted(legend_latxonomy_set):
		color = d_taxonomy_colors.get(taxonomy, (0, 0, 0))
		patch = Patch(facecolor=color, label=taxonomy)
		legend_elements.append(patch)

	max_x = max(pos.values(), key=lambda c: abs(c[0]))[0]
	max_y = max(pos.values(), key=lambda c: abs(c[1]))[1]
	maxV = max(max_x, max_y)

	ax.set_ylim(-maxV, maxV)
	ax.set_xlim(-maxV, maxV)

	ax.legend(handles=legend_elements, ncol=3, bbox_to_anchor=(1.1, 1.0), fancybox=True)				
				
	plt.savefig(options.output, transparent=True, dpi=300)

	plt.show()


if __name__ == '__main__':
	_main()


