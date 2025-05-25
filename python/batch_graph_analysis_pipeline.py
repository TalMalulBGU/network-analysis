import argparse
import json
import os
import glob
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

from tutils.pipelines import GraphAnalysisPipelineTemplate
from tutils.analyses.igraph_analysis import *
from tutils.analyses.analysis import *
from tutils.analyses.clustering_analysis import *
from tutils.settings import *
from tutils.records.record import *
from tutils.filters.filters import *
from tutils.graphs import *
from tutils.clustering import *
from tutils.databases import *

class SimplePipeline(GraphAnalysisPipelineTemplate):
	def __init__(self, xgmml_file, database_file, settings_file, output_folder, graph_filter_iterration, weights):
		super().__init__(xgmml_file, database_file, settings_file, output_folder, graph_filter_iterration, weights)
		
def _main():


	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--output_folder', help='output folder location', type=str, required=True)
	parser.add_argument('--settings_file', help='settings_file', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--itteration', help='database file location for filtering', type=int, required=False, default=0)
	parser.add_argument('--alignment_weights', help='database file location for filtering', type=str, required=False, default=None)
	

	options, unknown_args = parser.parse_known_args()
	
	
		

	print("output_folder = {}".format(options.output_folder))
	print("settings_file = {}".format(options.settings_file))

	kwargs = {}
	while len(unknown_args) > 0:
		k = args.pop(0)
		v = args.pop(0)
		kwargs[k] = v
	
	weights = None
	
	if options.alignment_weights is not None:
		with open(options.alignment_weights) as json_file:
			weights = json.load(json_file)
	
	simple_pipeline = SimplePipeline(options.xgmml_file, options.database, options.settings_file, output_folder = options.output_folder, graph_filter_iterration = options.itteration, weights = weights)	
	simple_pipeline.run()
	
if __name__ == '__main__':
	_main()


