import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import json
import ast
from Bio.pairwise2 import align
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import deque
import io
import tempfile

def array_converter(value):
	result = []

	if pd.notnull(value) and value != '':
		if ';' in value:
			result = str(value).strip().split(';')[:-1]
		else:
			result = ast.literal_eval(value)

	return list(result)

def load_db(location):

	uniprot_db = pd.read_csv(location)

	if len(uniprot_db) > 0:
		uniprot_db['Entry'] = uniprot_db['Entry'].astype(str)
		uniprot_db['Sequence'] = uniprot_db['Sequence'].astype(str)

	uniprot_db = uniprot_db.set_index('Entry')
	return uniprot_db



def _main():

	

	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--start', help='starting node index', type=int, required=True)
	parser.add_argument('--size', help='size from strating node', type=int, required=True)
	parser.add_argument('--alignment_name', help='the name for the alignment weight "_score" will be added at the end', type=str, default='alignment')
	parser.add_argument('--gap_open_penalty', help='the cost of opening a gap', type=int, default=0)
	parser.add_argument('--gap_extend_penalty', help='the cost of extending an open gap', type=int, default=-1)


	options = parser.parse_args()

	print("database = {}".format(options.database))
	print("gap_open_penalty = {}".format(options.gap_open_penalty))
	print("gap_extend_penalty = {}".format(options.gap_extend_penalty))
	print("output = {}".format(options.output))
	
	
	netwrok_sequences_db = None
	if options.database:
		 netwrok_sequences_db = load_db(options.database)



	matrix = matlist.blosum62
	alignment_scores = []
	stoping_index = options.start + options.size if options.start + options.size < len(netwrok_sequences_db.index) else len(netwrok_sequences_db.index)
	alignment_scores_index = 0
	pointer_dictionary = {}
	sequence_list = [SeqRecord(Seq(netwrok_sequences_db.loc[v_target,'Sequence']), id=v_target, description='') for v_target in netwrok_sequences_db.index]	

	
	for idx_soruce in range(options.start, stoping_index, 1):
		if idx_soruce + 1 == len(sequence_list):
			break
			
		source_file = tempfile.NamedTemporaryFile(mode='w', delete=False).name
		targets_file = tempfile.NamedTemporaryFile(mode='w', delete=False).name
		
		source_seq = sequence_list[idx_soruce]
		SeqIO.write(source_seq, source_file, "fasta")
		SeqIO.write(sequence_list[idx_soruce + 1:], targets_file, "fasta")
		
		output_source_targets = NcbiblastpCommandline(query=source_file, subject=targets_file, outfmt='5', gapopen=11, gapextend=1, num_descriptions=100000000, num_alignments=100000000, max_hsps=1)()[0]
		output_targets_source = NcbiblastpCommandline(query=targets_file, subject=source_file, outfmt='5', gapopen=11, gapextend=1, num_descriptions=100000000, num_alignments=100000000, max_hsps=1)()[0]
		blast_result_record_source_targets = NCBIXML.parse(io.StringIO(output_source_targets))
		blast_result_record_targets_source = NCBIXML.parse(io.StringIO(output_targets_source))

		for query in list(blast_result_record_source_targets):
			for alignment in query.alignments:
				for hsp in alignment.hsps:
					pointer_dictionary['{}_{}'.format(query.query, alignment.hit_id)] = alignment_scores_index
					alignment_scores_index+=1
					source_name = query.query
					target_name = alignment.hit_id
					source_align_seq = hsp.query
					target_align_seq = hsp.sbjct
					source_start = hsp.query_start
					source_end = hsp.query_end
					target_start = hsp.sbjct_start
					target_end = hsp.sbjct_end
					alignment_length = hsp.align_length
					alignment_identity = round((hsp.identities / hsp.align_length) * 100, 2)
					e_value = hsp.expect
					alignment_score = hsp.bits
					alignment_scores.append({'source': source_name, 'target': target_name, '{}_score'.format(options.alignment_name): alignment_score,
					'source_align_seq': source_align_seq, 'target_align_seq': target_align_seq,
					'source_start': source_start, 'source_end': source_end,
					'target_start': target_start, 'target_end': target_end,
					'alignment_length': alignment_length, '%{}_id'.format(options.alignment_name): alignment_identity,
					'e_value': e_value
					})
		for query in list(blast_result_record_targets_source):
			for alignment in query.alignments:
				for hsp in alignment.hsps:
					alignment_scores_index = pointer_dictionary.get('{}_{}'.format(alignment.hit_id, query.query), -1) 
					if alignment_scores_index != -1 and hsp.bits >= alignment_scores[alignment_scores_index]['{}_score'.format(options.alignment_name)]:
						continue
					source_name = query.query
					target_name = alignment.hit_id
					source_align_seq = hsp.query
					target_align_seq = hsp.sbjct
					source_start = hsp.query_start
					source_end = hsp.query_end
					target_start = hsp.sbjct_start
					target_end = hsp.sbjct_end
					alignment_length = hsp.align_length
					alignment_identity = round((hsp.identities / hsp.align_length) * 100, 2)
					e_value = hsp.expect
					alignment_score = hsp.bits
					alignment_scores[alignment_scores_index] = {'source': source_name, 'target': target_name, '{}_score'.format(options.alignment_name): alignment_score,
					'source_align_seq': source_align_seq, 'target_align_seq': target_align_seq,
					'source_start': source_start, 'source_end': source_end,
					'target_start': target_start, 'target_end': target_end,
					'alignment_length': alignment_length, '%{}_id'.format(options.alignment_name): alignment_identity,
					'e_value': e_value
					}


	json_object = json.dumps(alignment_scores, indent = 4)

	with open(options.output + '_' + str(options.start) + '_' + str(stoping_index), "w") as outfile:
		outfile.write(json_object)

if __name__ == '__main__':
	_main()


