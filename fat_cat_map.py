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
parser.add_option("-s", "--sequence", action="store", type="string", dest="sequence", default='', help="New sequence (FASTA)")
parser.add_option("-c", "--cpu", action="store", type="string", dest="cpu", default='1', help="Parallel threads to use")
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

## Main
def main():
	identifiers = parse_text_delimited(options.IDs, '#', 0)

	maxe = 1.0
	maxgene = 'agenestillhasnoname'

	for ID in identifiers:
		evalue, node = run_hmm(ID + '.hmm', options.sequence, options.cpu)

		if evalue < maxe:
			maxe = evalue
			maxnode = node

	genes = parse_msa(maxnode + '.msa')
	genes.remove('//')

	print genes
	
	
if __name__ == '__main__':
	main()
