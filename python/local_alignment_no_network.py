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
import importlib.util
import sys
from Bio.pairwise2 import align
from Bio.SubsMat import MatrixInfo as matlist
import editdistance
import vaex
import pyarrow
from decimal import *


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


directory_path = "/sise/vaksler-group/IsanaRNA/Tal/python/tutils"  # Replace with the path to your directory
create_modules_from_directory(directory_path)

from tutils.utils import *
from tutils.protein_clustering import *
from tutils.graphics.plots import *
from tutils.graphics.animations import *
from tutils.analyses.igraph_analysis import *
from tutils.analyses.igraph_protein_analysis import *
from tutils.analyses.analysis import *
from tutils.analyses.clustering_analysis import *
from tutils.settings import *
from tutils.records.record import *
from tutils.filters.filters import *
from tutils.graphs import *
from tutils.clustering import *
from tutils.databases import *
from tutils.relevancies import *



def _main():


	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--output_folder', help='output file location', type=str, required=True)
	parser.add_argument('--start', help='starting node index', type=int, required=True)
	parser.add_argument('--step', help='step from strating node', type=int, required=True)
	parser.add_argument('--alignment_name', help='the name for the alignment weight "_score" will be added at the end', type=str, default='alignment')
	parser.add_argument('--gap_open_penalty', help='the cost of opening a gap', type=int, default=0)
	parser.add_argument('--gap_extend_penalty', help='the cost of extending an open gap', type=int, default=-1)


	options = parser.parse_args()


	database = None
	if options.database:
		 database = UniprotDB(options.database)



	matrix = matlist.blosum62
	stoping_index = options.start + options.step if options.start + options.step < len(database.dataframe.index) else len(database.dataframe.index)
	hdf5s = []
	for idx_source in range(options.start, stoping_index, 1):
		alignment_scores = []
		source_name = database.dataframe.index[idx_source]
		output_path = '{}/{}_{}.hdf5'.format(options.output_folder, source_name, idx_source)
		
		if idx_source + 1 < len(database.dataframe.index):
			hdf5s.append(output_path)
		if os.path.exists(output_path):
			continue

		for idx_target in range(idx_source + 1, len(database.dataframe.index)):
			
			target_name = database.dataframe.index[idx_target]
			source_sequence = database.dataframe.loc[source_name,'Sequence']
			target_sequence = database.dataframe.loc[target_name,'Sequence']
			alignment = align.localds(source_sequence, target_sequence, matrix, options.gap_open_penalty, options.gap_extend_penalty, one_alignment_only=True)[0]

			alignment_seqA = alignment.seqA
			alignment_seqB = alignment.seqB
			alignment_start = alignment.start
			alignment_end = alignment.end
			alignment_length = len(alignment.seqA)
			alignment_identity = round(1 - editdistance.eval(alignment.seqA, alignment.seqB) / alignment_length, 4) * 100
			normalized_alignment_score = int(float(-Decimal.log10(Decimal('2') ** (Decimal('-1') * Decimal(str(alignment.score))) * Decimal(str(len(source_sequence))) * Decimal(str(len(target_sequence))))))
			alignment_scores.append({'source': source_name, 'target': target_name, '{}_score'.format(options.alignment_name): alignment.score, '%{}_id'.format(options.alignment_name): alignment_identity , 'normalized_{}_score'.format(options.alignment_name): normalized_alignment_score})


		df = pd.DataFrame(alignment_scores)
		vdf = vaex.from_pandas(df, name='{}_{}_{}'.format(options.protein, source_name, idx_source), copy_index=False)
		vdf.export_hdf5(output_path)
		
	vdf = vaex.open_many(hdf5s)
	vdf.export_hdf5('{}/{}_{}.hdf5'.format(options.output_folder, options.start, stoping_index))
	#json_object = json.dumps(alignment_scores, indent = 4)

	#with open(options.output + '_' + str(options.start) + '_' + str(stoping_index), "w") as outfile:
		#outfile.write(json_object)

if __name__ == '__main__':
	_main()


