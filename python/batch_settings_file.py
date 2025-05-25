import argparse
import json
import os
import glob
import importlib.util
import sys
from subprocess import run, PIPE


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

from tutils.settings import SettingParser

		
def _main():


	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--output_folder', help='output folder location', type=str, required=True)
	parser.add_argument('--settings_file', help='settings_file', type=str, required=True)
	parser.add_argument('--network_general', help='networkx general file', type=str, required=True)
	parser.add_argument('--network_config', help='network config file', type=str, required=True)
	

	options, unknown_args = parser.parse_known_args()

	print("output_folder = {}".format(options.output_folder))
	print("settings_file = {}".format(options.settings_file))
	print("network_general = {}".format(options.network_general))
	print("network_config = {}".format(options.network_config))
    
	os.environ['START'] = '46'
	os.environ['END'] = '217'
    
	kwargs = {}
	while len(unknown_args) > 0:
		k = args.pop(0)
		v = args.pop(0)
		kwargs[k] = v

	
	runner_batch_settings = '/home/talmalu/thesis/projects/python/NetworkAnalysis/Scripts/bash/runner_batch_settings.sh'
	with open(options.settings_file) as json_file:
		json_config = json.load(json_file)
		parser = SettingParser(json_config)

		filters = parser.parse_filters()

		command = ['bash', runner_batch_settings, '--network_general', options.network_general, '--network_config', options.network_config, '--output_folder', options.output_folder]
		result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
		print(result.stdout, result.stderr)

	

	
if __name__ == '__main__':
	_main()


