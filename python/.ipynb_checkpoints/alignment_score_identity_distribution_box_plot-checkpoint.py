import argparse
import matplotlib.pyplot as plt
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
	graph = ig.Graph.from_networkx(nx_g)
	graph.vs['name'] = graph.vs['_nx_name']
	nx_g.clear()

	return graph


def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--out_file', help='output file location', type=str, required=True)
	parser.add_argument('--data_file', help='name of edge attribute', type=str, required=True)


	options = parser.parse_args()
	
	uniprot_db = None
	my_root_protein_location = '/home/talmalu/thesis/data/' + options.protein + '/'
	my_root_output_location = my_root_protein_location + 'output/'
	

	print("protein = " + options.protein)
	print("out_file = " + options.out_file)
	print("data_file = {}".format(options.data_file))



	alignment_scores = []
	identities = []

	with open(options.data_file) as json_file:
		data = json.load(json_file)
		alignment_scores = list(map(int,data.keys()))
		identities = data.values()
		
	fig, ax = plt.subplots(figsize=(120,30))
	fig.suptitle('Distribution of Identity for each Normalized Global Alignment Score')
	ax.boxplot(identities, labels=alignment_scores, patch_artist=True, showfliers=False,
			   whiskerprops=dict(linestyle='--', color='#777777'),
			  medianprops=dict(linewidth=2.5, color='red'),
			  boxprops=dict(linewidth=0, facecolor='purple'))

	ax.axhline(35, color='r', alpha=0.7, linestyle='--')



	ax.set_xlabel('Alignment Score')
	ax.set_ylabel('Identity')
	ax.set_yticks(np.arange(0, 101, 5))



	for tick in ax.get_xticklabels():
		tick.set_rotation(90)



	plt.savefig(options.out_file, bbox_inches='tight')
	plt.show()



if __name__ == '__main__':
	_main()


