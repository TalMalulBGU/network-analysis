import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import igraph as ig
import json
import ast
import os


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
	
	
	if db:
		db_missing_entries = set(graph.vs['name']) - set(db.index)
		graph.delete_vertices(graph.vs.select(lambda v: v['name'] in db_missing_entries))
		graph_fragments_entries = set(graph.vs['name']).intersection(set(db[(db['Fragment'])].index))
		graph.delete_vertices(graph.vs.select(lambda v: v['name'] in graph_fragments_entries))

	return graph


def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('protein', help='the name of the protien', type=str)
	parser.add_argument('xgmml', help='output file location', type=str)
	parser.add_argument('output', help='output file location', type=str)
	parser.add_argument('--database', help='database file location for filtering', type=str)
	parser.add_argument('--weights', help='name of edge attribute', type=str, default=None)


	options = parser.parse_args()
	
	my_root_protein_location = '/home/talmalu/thesis/data/' + options.protein + '/'
	my_root_output_location = my_root_protein_location + 'output/'

	print("database = " + options.database)
	print("graph = " + options.xgmml)
	print("protein = " + options.protein)
	print("output = " + options.output)
	print("weights = " + str(options.weights))


	

	uniprot_db = None
	if options.database:
		uniprot_db = load_db(options.database)
		
	graph = load_network(options.xgmml, uniprot_db)

	global_alignment_score_path = '{}{}'.format(my_root_output_location,'gloabl_alignment_score_with_normal.csv')

	df_global_alignment_scores = pd.read_csv(global_alignment_score_path)
	df_global_alignment_scores = df_global_alignment_scores.set_index(['source', 'target'])

	   
	edges_global_alignment_scores = []
	edges_normalized_global_alignment_scores = []

	for edge in graph.es:
		source = graph.vs[edge.source]['name']
		target = graph.vs[edge.target]['name']

		try:
			score = df_global_alignment_scores.loc[(source, target),'global_alignment_score']
			normalized_score = df_global_alignment_scores.loc[(source, target), 'normalized_global_alignment_score']
		except:
			score = df_global_alignment_scores.loc[(target, source),'global_alignment_score']
			normalized_score = df_global_alignment_scores.loc[(target, source), 'normalized_global_alignment_score']

		if score < 0:
			score = 0

		if normalized_score < 0:
			normalized_score = 0
		
		edges_global_alignment_scores.append(score)
		edges_normalized_global_alignment_scores.append(normalized_score)


	graph.es['global_alignment_score'] = edges_global_alignment_scores
	graph.es['normalized_global_alignment_score'] = edges_normalized_global_alignment_scores
	
	
	alignment_scores = []
	lengths = []

	global_alignment_scores = sorted(set(graph.es[options.weights]))


	for alignment_score in global_alignment_scores:
		alignment_score_entries = set([entry for source_target in df_global_alignment_scores[df_global_alignment_scores[options.weights] == alignment_score].index for entry in source_target])
		alignment_score_entries_length = uniprot_db.loc[alignment_score_entries, 'Length']
		alignment_scores.append(alignment_score)
		lengths.append(list(alignment_score_entries_length))

	
	result_dictionary = {}
	for index, alignment_score in enumerate(alignment_scores):
		result_dictionary[int(alignment_score)] = lengths[index]

	json_object = json.dumps(result_dictionary, indent = 4)

	with open(options.output, "w") as outfile:
		outfile.write(json_object)

if __name__ == '__main__':
	_main()


