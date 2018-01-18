#!/usr/bin/python

import optparse
from optparse import OptionParser 

import os.path, subprocess

import sets

import string

from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment

## OptionParser
# import arguments and options
usage = "usage: %prog -i identifiers -t phylogenetic_tree -m multiple_sequence_alignment -s new_sequence --cpu [number of processors]"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--IDs", action="store", type="string", dest="IDs", default='', help="Node identifiers file (text)")
parser.add_option("-t", "--tree", action="store", type="string", dest="tree", default='', help="Phylogenetic tree (Newick format)")
parser.add_option("-m", "--msa", action="store", type="string", dest="msa", default='', help="Multiple sequence alignment (FASTA)")
parser.add_option("-s", "--sequence", action="store", type="string", dest="sequence", default='', help="New sequence (FASTA)")
parser.add_option("-c", "--cpu", action="store", type="string", dest="cpu", default='1', help="Parallel threads to use")
(options, args) = parser.parse_args()

def parse_text_delimited(file_name, exclude, index):
	information = []
	Fopen = open(file_name, 'r')
	for line in Fopen.readlines():
		if len(line) > 0:
			if line[0] != exclude:
				sline = string.split(line)
				information.append(sline[index])
	Fopen.close()
	return information

def build_msa(node, sequence_msa_map):
	key = str(node.__hash__())
	file_name = key + '.msa'
	file_handle = open(file_name, 'w')
	terminals = node.get_terminals()
	alignments = [sequence_msa_map[terminal.name] for terminal in
			terminals]
	alignments = MultipleSeqAlignment(alignments)
	AlignIO.write(alignments, file_handle, 'stockholm')
	return file_name

def build_hmm(msa_file_name, threads):
	file_prefix = os.path.splitext(msa_file_name)[0]
	hmm_file_name = ''.join([file_prefix, '.hmm'])
	process_name = "hmmbuild --cpu %s %s %s > %s" % (threads, hmm_file_name, msa_file_name, 'temp.txt') 
	process = subprocess.Popen(process_name, shell = True)
	process.wait()
	return hmm_file_name

def parse_hmm(hmmsearch_file_name):
	evalues = parse_text_delimited(hmmsearch_file_name, '#', 4)
	if len(evalues) > 0:
		return float(evalues[0])
	else:
		return 1.0

def run_hmm(hmm_file_name, sequence_file_name, threads):
	search_file_name = string.split(hmm_file_name, '.')[0] + '_hmmsearch.txt'
	process_name = "hmmsearch --tblout %s --cpu %s %s %s > %s" % (search_file_name, threads, hmm_file_name, sequence_file_name, 'temp.txt') 
	process = subprocess.Popen(process_name, shell = True)
	process.wait()
	return parse_hmm(search_file_name), string.split(hmm_file_name, '.')[0]

def treewalker(root, sequence_msa_map, sequence_file_name, threads, identifiers, maxe, maxnode):
	msa_file_name = build_msa(root, sequence_msa_map)
	hmm_file_name = build_hmm(msa_file_name, threads)
	evalue, node = run_hmm(hmm_file_name, sequence_file_name, threads)

	if root.clades:
		for clade in root.clades:
			if evalue < maxe:
				maxe, maxnode, identifiers = treewalker(clade, sequence_msa_map, sequence_file_name, threads, set_merge(identifiers, [node]), evalue, node)
			else:
				maxe, maxnode, identifiers = treewalker(clade, sequence_msa_map, sequence_file_name, threads, set_merge(identifiers, [node]), maxe, maxnode)

	evalues = [evalue, maxe]
	nodes = [node, maxnode]

	return min(evalues), nodes[evalues.index(min(evalues))], set_merge(identifiers, [node])

def parse_msa(msa_file):
	return parse_text_delimited(msa_file, '#', 0)

def set_merge(set1, set2):
	return list(sets.Set(set1) | sets.Set(set2))

def export_vector(vector_file, vector):
	vector_outfile = open(vector_file, 'w')
	for element in vector:
		vector_outfile.write(element + '\n')
	vector_outfile.close()
	return

def main():
	tree = Phylo.read(options.tree, "newick")
	msa = AlignIO.read(options.msa, "fasta")
	
	sequence_msa_map = {}
	for entry in msa:
		sequence_msa_map[entry.id] = entry
	
	evalue, node, identifiers = treewalker(tree.root, sequence_msa_map, options.sequence, options.cpu, [], 1.0, 'agenehasnoname')

	export_vector(options.IDs, identifiers)

	genes = parse_msa(node + '.msa')
	genes.remove('//')

	print genes
	
	
if __name__ == '__main__':
	main()
