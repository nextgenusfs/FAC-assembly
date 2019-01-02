#!/usr/bin/env python

import sys
import os
import argparse
import mappy
import platform
import datetime
from natsort import natsorted
from Bio import SeqIO

__version__ = "0.2.0"

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
parser.add_argument('--min_len', type=int, default=20000, help='Minimum contig length to find end primers')
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
                result[cols[2]] = {'organism': cols[0], 'order': cols[1], 'fwd': cols[3],
                                   'rev': cols[4], 'column': cols[5], 'row': cols[6], 'plate': cols[7]}
    return result

def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)
    
print('------------------------------------------')
print('[{:}] Running Python v{:} '.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), platform.python_version()))

#parse the text file
inputDict = parseInput(args.input)
print('[{:}] Loading data for {:,} FACs from {:}'.format(
      datetime.datetime.now().strftime('%b %d %I:%M %p'), len(inputDict), args.input))

trimDict = {}
if args.pool:
    basename = '-'.join(args.pool)
    for k,v in inputDict.items():
        if v['column'] in args.pool or v['row'] in args.pool or v['plate'] in args.pool:
            trimDict[k] = v
    print('[{:}] Limiting search for {:,} FACs from {:} pools'.format(
          datetime.datetime.now().strftime('%b %d %I:%M %p'), len(trimDict), ' '.join(args.pool)))
else:
    basename = os.path.basename(args.assembly).split('.f')[0]
    trimDict = inputDict

print('[{:}] Building minimap2 index from {:,} contigs from {:}'.format(
      datetime.datetime.now().strftime('%b %d %I:%M %p'), countfasta(args.assembly), args.assembly))

a = mappy.Aligner(args.assembly)
PEhits = {}
SEhits = {}
PrimerHits = {}
print('[{:}] Mapping Sanger reads to assembly using minimap2'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
mapData = basename+'.map.txt'
with open(mapData, 'w') as rawdata:
    header = ['FAC_name', 'read_num', 'query_len', 'query_start', 'query_end', 'query_cov', 'strand', 'contig',
              'contig_len', 'ref_start', 'ref_end', 'match_len', 'aln_len', 'pident', 'map_qual', 'cigar']
    rawdata.write('{:}\n'.format('\t'.join(header)))
    for k,v in trimDict.items():
        mapped = None
        if len(v['fwd']) > 50 and len(v['fwd']) > 50:
            mapped = a.map(v['fwd'], v['rev'])
        else:
            if len(v['fwd']) > 50:
                mapped = a.map(v['fwd'], '')
            else:
                mapped = a.map('', v['rev'])
        if not mapped:
            print('No matches for primers in {:}'.format(k))
            continue
            
        for align in mapped:
            pident = align.mlen / float(align.blen) * 100
            #print(dir(align))
            if align.read_num == 1:
                readLen = len(v['fwd'])
            else:
                readLen = len(v['rev'])
            queryCov = (align.q_en - align.q_st) / float(readLen) * 100
            if not align.is_primary or align.mapq < 40 or align.ctg_len < args.min_len or queryCov < 50.00:
                continue
            rawdata.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:.1f}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                          k, align.read_num, readLen, align.q_st, align.q_en, queryCov, align.strand, align.ctg, 
                          align.ctg_len, align.r_st, align.r_en, align.mlen, align.blen, pident, align.mapq, align.cigar_str))
            if not k in PrimerHits:
                PrimerHits[k] = {align.read_num: [(align.ctg, align.strand, align.ctg_len, readLen, align.r_st, align.r_en, pident)]}
            else:
                if not align.read_num in PrimerHits[k]:
                    PrimerHits[k][align.read_num] = [(align.ctg, align.strand, align.ctg_len, readLen, align.r_st, align.r_en, pident)]
                else:
                    PrimerHits[k][align.read_num].append((align.ctg, align.strand, align.ctg_len, readLen, align.r_st, align.r_en, pident))

Classified = {}
print('[{:}] Parsing results, and outputting trimmed contigs'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
SeqRecords = SeqIO.to_dict(SeqIO.parse(args.assembly, 'fasta'))
sortedKeyList = sorted(PrimerHits.keys(), key=lambda s: len(PrimerHits.get(s)), reverse=True)
Counts = {'full':0, 'stitch':0,'single':0}
Species = {'full': {}, 'partial': {}}
organismList = []
with open(basename+'.all-contigs.fasta', 'w') as outfile:
	with open(basename+'.fulllength-contigs.fasta', 'w') as outfile2:
		for k in sortedKeyList:
			v = PrimerHits[k]
			organism = trimDict[k]['organism']
			if not organism in organismList:
				organismList.append(organism)
			if len(v) > 1:
				if v[1][0][0] == v[2][0][0]: #these are full length, so re-orient and write
					if v[1][0][1] == v[2][0][1]: #sanger reads map in wrong orientation
						print('{:} could be missassembled - orientation of {:} Sanger reads is in conflict'.format(v[1][0][0], k))
						print(k,v)
						continue
					elif v[1][0][1] == 1 and v[2][0][1] == -1:
						SeqRec = SeqRecords[v[1][0][0]][v[1][0][4]:v[2][0][5]]
						Seq = str(SeqRec.seq)
					elif v[1][0][1] == -1 and v[2][0][1] == 1: #need to reverse complement after slicing
						SeqRec = SeqRecords[v[1][0][0]][v[2][0][4]:v[1][0][5]]
						SeqRecRev = SeqRec.reverse_complement()
						Seq = str(SeqRecRev.seq)
					outfile.write('>{:};organism={:};origin={:};length={:};full_length=yes;\n{:}\n'.format(k, organism, v[1][0][0], len(Seq), softwrap(Seq)))
					outfile2.write('>{:};organism={:};origin={:};length={:};full_length=yes;\n{:}\n'.format(k, organism, v[1][0][0], len(Seq), softwrap(Seq)))
					Counts['full'] += 1
					if organism in Species['full']:
						Species['full'][organism] += 1
					else:
						Species['full'][organism] = 1
				else: #contigs are different, so stitch them together with 100 N's in between
					Seq = ''
					#add the "forward" sequence
					if v[1][0][1] == 1:
						SeqRec = SeqRecords[v[1][0][0]][v[1][0][4]:]
						Seq += str(SeqRec.seq)
					else:
						SeqRec = SeqRecords[v[1][0][0]][:v[1][0][5]]
						SeqRecRev = SeqRec.reverse_complement()
						Seq += str(SeqRecRev.seq)
					#add the 100 bp linker
					Seq += 'N'*100
					#add "reverse" sequence
					if v[2][0][1] == -1:
						SeqRec = SeqRecords[v[2][0][0]][:v[2][0][5]]
						Seq += str(SeqRec.seq)
					else:
						SeqRec = SeqRecords[v[2][0][0]][v[2][0][4]:]
						SeqRecRev = SeqRec.reverse_complement()
						Seq += str(SeqRecRev.seq)
					outfile.write('>{:};organism={:};origin={:}-{:};length={:};full_length=stitched;\n{:}\n'.format(k, organism, v[1][0][0], v[2][0][0], len(Seq), softwrap(Seq)))
					Counts['stitch'] += 1
					if organism in Species['partial']:
						Species['partial'][organism] += 1
					else:
						Species['partial'][organism] = 1
			else:
				if 1 in v: #forward primer found
					if v[1][0][1] == 1:
						SeqRec = SeqRecords[v[1][0][0]][v[1][0][4]:]
						Seq = str(SeqRec.seq)
					elif v[1][0][1] == -1:
						SeqRec = SeqRecords[v[1][0][0]][:v[1][0][5]]
						SeqRecRev = SeqRec.reverse_complement()
						Seq = str(SeqRecRev.seq)
					Counts['single'] += 1
					outfile.write('>{:};organism={:};origin={:};length={:};full_length=no;\n{:}\n'.format(k, organism, v[1][0][0], len(Seq), softwrap(Seq)))
				elif 2 in v: #reverse found
					if v[2][0][1] == -1:
						SeqRec = SeqRecords[v[2][0][0]][:v[2][0][5]]
						Seq = str(SeqRec.seq)
					elif v[2][0][1] == 1:
						SeqRec = SeqRecords[v[2][0][0]][v[2][0][4]:]
						SeqRecRev = SeqRec.reverse_complement()
						Seq = str(SeqRecRev.seq)
					Counts['single'] += 1
					if organism in Species['partial']:
						Species['partial'][organism] += 1
					else:
						Species['partial'][organism] = 1		 
					outfile.write('>{:};organism={:};origin={:};length={:};full_length=no;\n{:}\n'.format(k, organism, v[2][0][0], len(Seq), softwrap(Seq)))
print('[{:}] Identified {:,} full-length contigs, {:,} stitched contigs, and {:,} with single primer'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), Counts['full'], Counts['stitch'], Counts['single']))
print('[{:}] Contig counts by organism:'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
print('\033[4m  Species			Full	Partial\033[0m')
for sp in organismList:
	if not sp in Species['full']:
		Species['full'][sp] = 0
	if not sp in Species['partial']:
		Species['partial'][sp] = 0
	print('  {:}		{:}	{:}'.format(sp, Species['full'][sp], Species['partial'][sp]))
print('[{:}] Output files:'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
print('  Raw map data: {:}'.format(basename+'.map.txt'))
print('  Full-length:  {:}'.format(basename+'.fulllength-contigs.fasta'))
print('  All contigs:  {:}'.format(basename+'.all-contigs.fasta'))
sys.exit(1)
            