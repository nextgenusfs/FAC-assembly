#!/usr/bin/env python

import sys
import os
import argparse
import platform
import datetime

__version__ = "0.2.0"

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='filter_seqs.py',
    description = '''Script filters FASTA file.''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', required=True, help='TSV FAC end sequences')
parser.add_argument('-m', '--minlen', type=int, help='Minimum length of contig')
parser.add_argument('-s', '--species', help='Species to keep')
parser.add_argument('-l', '--longest', action='store_true', help='Keep only the longest')
parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
args=parser.parse_args()

def SimpleFastaParser(handle):
    # Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break
    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in Fasta files should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = handle.readline()

        yield title, "".join(lines).replace(" ", "").replace("\r", "")

        if not line:
            return  # StopIteration
    assert False, "Should not reach this line"

def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

longest = {}
SeqRecords = {}
with open(args.input, 'rU') as infile:
	for header, seq in SimpleFastaParser(infile):
		seqlen = len(seq)
		if not args.longest:
			if args.minlen:
				if seqlen < args.minlen:
					continue
			if args.species:
				if not args.species in header:
					continue
			sys.stdout.write('>{:}\n{:}\n'.format(header, softwrap(seq)))
		else: #bit more complicated as have to run through entire fasta file before writing
			SeqRecords[header] = seq
			parts = header.split(';')
			for x in parts:
				if x.startswith('origin='):
					contig = x.replace('origin=', '')
			if not contig in longest:
				longest[contig] = {'length': seqlen, 'id': header}
			else:
				if seqlen > longest[contig]['length']:
					longest[contig] = {'length': seqlen, 'id': header}

if args.longest:
	for k,v in longest.items():
		header = longest[k]['id']
		seq = SeqRecords[header]
		if args.minlen:
			if len(seq) < args.minlen:
				continue
		if args.species:
			if not args.species in header:
				continue
		sys.stdout.write('>{:}\n{:}\n'.format(header, softwrap(seq)))
				
			
	
	