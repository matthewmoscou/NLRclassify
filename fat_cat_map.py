#!/usr/bin/python

## Modules
import optparse
from optparse import OptionParser 

import os.path, subprocess

import string

from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment

## OptionParser
# import arguments and options
usage = "usage: %prog -i identifiers -s new_sequence -c [number of processers]"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--IDs", action="store", type="string", dest="IDs", default='', help="Identifiers of prebuilt FAT-CAT run")
parser.add_option("-m", "--membership", action="store", type="string", dest="membership", default='', help="Clade membership")
parser.add_option("-s", "--sequence", action="store", type="string", dest="sequence", default='', help="New sequence (FASTA)")
parser.add_option("-c", "--cpu", action="store", type="string", dest="cpu", default='1', help="Parallel threads to use")
parser.add_option("-p", "--previous", action="store", type="string", dest="previous", default='', help="Skip evaluations from previous run")
(options, args) = parser.parse_args()

## Functions
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

def parse_msa(msa_file):
	return parse_text_delimited(msa_file, '#', 0)

def parse_dictionary(file_name, index, reference):
	dictionary = {}
	Fopen = open(file_name, 'r')
	for line in Fopen.readlines():
		if len(line) > 0:
			sline = string.split(line)
			dictionary[sline[index]] = sline[reference]
	Fopen.close()
	return dictionary

def parse_FASTA(FASTA_file_name):
	FASTA_file = open(FASTA_file_name, 'r')
	gene_sequence = {}
	for line in FASTA_file.readlines():
		sline = string.split(line)
		if len(line) > 0:
			if line[0] == '>':
				ID = sline[0][1:]
				gene_sequence[ID] = ''
			else:
				gene_sequence[ID] += sline[0]
	FASTA_file.close()
	return gene_sequence

def export_sequence(FASTA_file_name, identifier, sequence):
	Fopen = open(FASTA_file_name, 'w')
	Fopen.write('>' + identifier + '\n')
	Fopen.write(sequence + '\n')
	Fopen.close()

## Main
def main():
	identifiers = parse_text_delimited(options.IDs, '#', 0)
	gene_sequence = parse_FASTA(options.sequence)

	if len(options.membership) > 0:
		clade_membership = parse_dictionary(options.membership, 0, 1)

	if len(options.previous) > 0:
		previously_evaluated = parse_text_delimited(options.previous, '#', 0)
	else:
		previously_evaluated = []
		
	for geneID in gene_sequence.keys():
		if geneID not in previously_evaluated:
			export_sequence('temp.fa', geneID, gene_sequence[geneID])

			maxe = 1.0
			maxnode = 'agenestillhasnoname'

			for ID in identifiers:
				evalue, node = run_hmm(ID + '.hmm', 'temp.fa', options.cpu)

				if evalue < maxe:
					maxe = evalue
					maxnode = node

			if maxnode != 'agenestillhasnoname':
				target_genes = parse_msa(maxnode + '.msa')
				target_genes.remove('//')
			else:
				target_genes = ['No hits']
			
			if len(options.membership) > 0:
				for gene in target_genes:
					if gene in clade_membership.keys():
						print geneID + '\t' + gene + '\t' + clade_membership[gene]
					else:
						print geneID + '\t' + gene
			else:
				for gene in target_genes:
					print geneID + '\t' + gene
	
	
if __name__ == '__main__':
	main()
