#!/usr/bin/env python

import sys
import os
import argparse
import mappy
import platform
import datetime
from natsort import natsorted

__version__ = "0.0.3"

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='get_full_length_facs.py',
    description = '''Script that parses tsv file, maps end sequences to FAC assembly.''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-a', '--assembly', required=True, help='Assembly in FASTA format')
parser.add_argument('-i', '--input', required=True, help='TSV FAC end sequences')
parser.add_argument('-p', '--pool', nargs='+', help='Name of seq pool')
parser.add_argument('--min_len_complete', type=int, default=75000, help='Name of seq pool')
parser.add_argument('-v', '--debug', action='store_true', help='Debug print some debugging')
parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
args=parser.parse_args()

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def parseInput(input):
    result = {}
    with open(input, 'rU') as infile:
        for line in infile:
            if line.startswith('#') or line.startswith('\n'):
                continue
            line = line.strip()
            cols = line.split('\t')
            if len(cols) < 8:
                continue
            if not cols[2] in result:
                result[cols[2]] = {'organism': cols[0], 'order': cols[1], 'fwd': cols[3], 'rev': cols[4], 'column': cols[5], 'row': cols[6], 'plate': cols[7]}
    return result


print('------------------------------------------')
print('[{:}] Running Python v{:} '.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), platform.python_version()))

#parse the text file
inputDict = parseInput(args.input)
print('[{:}] Loading data for {:,} FACs from {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), len(inputDict), args.input))

trimDict = {}
if args.pool:
    basename = '-'.join(args.pool)
    for k,v in inputDict.items():
        if v['column'] in args.pool or v['row'] in args.pool or v['plate'] in args.pool:
            trimDict[k] = v
    print('[{:}] Limiting search for {:,} FACs from {:} pools'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), len(trimDict), ' '.join(args.pool)))
else:
    basename = os.path.basename(args.assembly).split('.f')[0]
    trimDict = inputDict
print('[{:}] Building minimap2 index from {:,} contigs from {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), countfasta(args.assembly), args.assembly))

a = mappy.Aligner(args.assembly, best_n=1)
PEhits = {}
SEhits = {}
print('[{:}] Mapping Sanger reads to assembly using minimap2'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
mapData = basename+'.map.txt'
with open(mapData, 'w') as rawdata:
	header = ['FAC_name', 'read_num', 'query_start', 'query_end', 'strand', 'contig', 'contig_len', 'ref_start', 'ref_end', 'match_len', 'aln_len', 'map_qual', 'cigar']
	rawdata.write('{:}\n'.format('\t'.join(header)))
	for k,v in trimDict.items():
		#paired
		if len(v['fwd']) > 50 and len(v['rev']) > 50:
			for align in a.map(v['fwd'], v['rev']):
				if int(align.mapq) < 60:
					continue
				rawdata.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(k, align.read_num, align.q_st, align.q_en, align.strand, align.ctg, align.ctg_len, align.r_st, align.r_en, align.mlen, align.blen, align.mapq, align.cigar_str))
				if int(align.ctg_len) < args.min_len_complete:
					continue
				if not k in PEhits:
					PEhits[k] = {align.read_num: (align.ctg, align.strand, align.ctg_len, len(v['fwd']))}
				else:
					if not align.read_num in PEhits[k]:
						PEhits[k][align.read_num] = (align.ctg, align.strand, align.ctg_len, len(v['rev']))
		#since some Sanger Sequences are missing, allow a single hit to be a "full length" if 
		elif len(v['fwd']) > 50:
			for align in a.map(v['fwd']):
				if int(align.mapq) < 60:
					continue
				rawdata.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(k, align.read_num, align.q_st, align.q_en, align.strand, align.ctg, align.ctg_len, align.r_st, align.r_en, align.mlen, align.blen, align.mapq, align.cigar_str))
				if int(align.ctg_len) < args.min_len_complete:
					continue
				if not k in SEhits:
					SEhits[k] = (align.ctg, align.strand, align.ctg_len, len(v['fwd']))
		elif len(v['rev']) > 50:
			for align in a.map(v['rev']):
				if int(align.mapq) < 60:
					continue
				rawdata.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(k, align.read_num, align.q_st, align.q_en, align.strand, align.ctg, align.ctg_len, align.r_st, align.r_en, align.mlen, align.blen, align.mapq, align.cigar_str))
				if int(align.ctg_len) < args.min_len_complete:
					continue
				if not k in SEhits:
					SEhits[k] = (align.ctg, align.strand, align.ctg_len, len(v['rev']))
#loop through PE hits first
completeFACs = {}
orientation = {}
flags = {}
for k,v in PEhits.items():
    if len(v) > 1:
        if v[1][0] == v[2][0] and v[1][1] != v[2][1]:
            if not k in orientation:
                orientation[k] = v[1][1]
            if not v[1][0] in completeFACs:
                completeFACs[v[1][0]] = [k]
            else:
                completeFACs[v[1][0]].append(k)
        else:
            if v[1][0] == v[2][0]: #same contig but orientation is wrong
                if not k in flags:
                    flags[k] = ['primer_orientation_mismatch']
                else:
                    flags[k].append('primer_orientation_mismatch')
                if not v[1][0] in completeFACs:
                    completeFACs[v[1][0]] = [k]
                else:
                    completeFACs[v[1][0]].append(k)             
#now loop through SE hits
for k,v in SEhits.items():
    if len(v) > 1:
        if not k in orientation:
            orientation[k] = v[1]
        if not v[0] in completeFACs:
            completeFACs[v[0]] = [k]
        else:
            completeFACs[v[0]].append(k)
#print(completeFACs)
#print(orientation)
if args.debug:
	for k,v in natsorted(flags.items()):
		print(k,v)
Complete = basename+'_complete_facs.fasta'
Unassigned = basename+'_unassigned_facs.fasta'
complete_count = 0
the_rest = 0
with open(Complete, 'w') as comp:
    with open(Unassigned, 'w') as bad:
        for seq in mappy.fastx_read(os.path.abspath(args.assembly), read_comment=True):
            if seq[0] in completeFACs:
                if seq[0] in orientation:
                    if orientation[seq[0]] == -1:
                        FinalSeq = mappy.revcomp(seq[1])
                    else:
                        FinalSeq = seq[1]
                else:
                    FinalSeq = seq[1]
                comments = seq[3].split(' ')
                facname = completeFACs[seq[0]][0]
                complete_count += len(completeFACs[seq[0]])
                comp.write('>{:};organism={:};{:};\n{:}\n'.format('|'.join(completeFACs[seq[0]]), inputDict[facname]['organism'], ';'.join(comments), FinalSeq))
            else:
            	the_rest += 1
                bad.write('>{:} {:}\n{:}\n'.format(seq[0], seq[3], seq[1]))
print('[{:}] Found {:,} full-length FACs corresponing to {:} unique sequences'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), complete_count, countfasta(Complete)))
print('[{:}] Full-length sequences: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), Complete))
print('[{:}] Unassigned sequences: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), Unassigned))
print('[{:}] Raw alignment data: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), mapData))
print('------------------------------------------')