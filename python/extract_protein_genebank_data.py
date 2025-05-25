import os
import io
import tempfile
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from xml.etree import ElementTree as ET
import argparse
import json


def query_protein_by_id(protein_id):
	Entrez.email = 'your_email@example.com'  # Enter your email here
	handle = Entrez.efetch(db='protein', id=protein_id, rettype='fasta', retmode='text')
	protein_record = SeqIO.read(handle, 'fasta')
	handle.close()

	features_handle = Entrez.efetch(db='protein', id=protein_id, rettype='gb', retmode='text')
	features_record = SeqIO.read(features_handle, 'gb')
	features_handle.close()
	
	handle = Entrez.efetch(db="popset", id=protein_id, rettype="gb", retmode="xml")
	record = handle.read()
	handle.close()

	
	taxonomy_levels = ['Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
	taxonomies = { taxonomy_level: [] for taxonomy_level in taxonomy_levels }
	
	taxonomy_name = None
	taxid = None
	
	root = ET.fromstring(record)
	for elem in root.iter("GBSeq_organism"):
		taxonomy_name = elem.text
		break
	
	
	for feature in features_record.features:
		if feature.type == 'source' and 'organism' in feature.qualifiers and 'db_xref' in feature.qualifiers:
			if feature.qualifiers['organism'][0].lower() == taxonomy_name.lower():
				for db_xref in feature.qualifiers['db_xref']:
					if db_xref.startswith('taxon:'):
						taxid = db_xref.lstrip('taxon:')
						break
	
	taxonomies['Species'] = [taxonomy_name.capitalize()]
	taxonomy_handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
	taxonomy_record = Entrez.read(taxonomy_handle)
	taxonomy_handle.close()
	taxonomy_dict = dict(taxonomy_record[0])
	for dict_taxonomy_level in taxonomy_dict['LineageEx']:
		taxonomy_level_rank = dict_taxonomy_level['Rank'].capitalize() if not dict_taxonomy_level['Rank'].capitalize() == 'Clade' else 'Kingdom'
		if taxonomy_level_rank not in taxonomy_levels:
			continue
		taxonomy_level_id = dict_taxonomy_level['TaxId']
		taxonomy_level_name = dict_taxonomy_level['ScientificName']
		taxonomies[taxonomy_level_rank].append(taxonomy_level_name)
		
	taxonomies = { taxonomy_level: lst if len(lst) > 0 else ['None'] for taxonomy_level, lst in taxonomies.items() }

	return protein_record, features_record, taxonomies, taxid

def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--proteins_file', help='output file location', type=str, required=True)
	parser.add_argument('--region_names_rule', help='output file location', nargs='*', required=False, default=[])
	

	options, unknown_args = parser.parse_known_args()
	
	proteins = []
	failed_proteins = []
	data_dict = {}

	print(options.proteins_file)
	with open(options.proteins_file, 'r') as file:
		lines = file.readlines()
		proteins = [line.strip() for line in lines]

	
	
	for protein_id in proteins:
		try:
			protein_record, features_record, taxonomies, taxid = query_protein_by_id(protein_id)
		except:
			failed_proteins.append(protein_id)
			continue
		
		gene_name = None
		for feature in features_record.features:
			if feature.type == 'CDS' and 'gene' in feature.qualifiers:
				gene_name = feature.qualifiers['gene'][0]
				break
			elif feature.type == 'gene':
				if 'gene' in feature.qualifiers:
					gene_name = feature.qualifiers['gene'][0]
					break
				elif 'name' in feature.qualifiers:
					gene_name = feature.qualifiers['name'][0]
					break
	   
		gene_name = gene_name.replace(" ", "_")
		data_dict[(protein_record.id, gene_name)] = {}
		data_dict[(protein_record.id, gene_name)]['Organism ID'] = taxid
		data_dict[(protein_record.id, gene_name)]['Gene names'] = [gene_name]
		data_dict[(protein_record.id, gene_name)]['taxonomies'] = taxonomies
		
		data_dict[(protein_record.id, gene_name)]['features'] = []
		for feature in features_record.features:
			### this if is for spesific region Rule
			if 'region_name' in feature.qualifiers and (len(options.region_names_rule) == 0 or any(region_name.lower() in options.region_names_rule for region_name in feature.qualifiers['region_name'])):
				sequence = str(feature.extract(protein_record.seq))
				data_dict[(protein_record.id, gene_name)]['features'].append({
					'region_name': feature.qualifiers['region_name'][0],
					'Length': len(sequence), 
					'Sequence': sequence, 
				})
		
	result_dictionary = { 
		str((protein_id[0], str(protein_id[1]), at_index + 1)): { **data_dict[protein_id]['features'][at_index], 
															'Organism': data_dict[protein_id]['taxonomies']['Species'][0],
															'Organism ID': data_dict[protein_id]['Organism ID'],
															'Gene names': data_dict[protein_id]['Gene names'],
															**data_dict[protein_id]['taxonomies']
														   }
		for protein_id in data_dict.keys()
		for at_index in range(len(data_dict[protein_id]['features']))
	}
	
	json_object = json.dumps(result_dictionary, indent = 4)
	print("for some reason for the proteins bellow retrvie data from server have failed, please try again")
	for protein in failed_proteins:
		print(protein)

	with open(options.output, "w") as outfile:
		outfile.write(json_object)
	
	
if __name__ == '__main__':
	_main()


