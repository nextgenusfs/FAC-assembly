#!/usr/bin/env python

import sys
import os
from interlap import InterLap
from collections import defaultdict
from natsort import natsorted
from Bio import SeqIO

#antismash GBK file to TSV format
#functions
def getID(input, type):
    #function to get ID from genbank record.features
    locusTag = None
    ID = None
    Parent = None
    if type == 'gene':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
        else:
            try:
                ID = input.qualifiers['gene'][0]
            except KeyError:
                pass
        return locusTag, ID, locusTag
        
    elif type == 'mRNA' or type == 'tRNA' or type == 'ncRNA' or type == 'rRNA' or type == 'exon':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
            if locusTag:
                Parent = locusTag
                try:
                    ID = input.qualifiers['transcript_id'][0]
                except KeyError:
                    pass
            else:
                try:
                    locusTag = input.qualifiers['transcript_id'][0]
                    Parent = locusTag
                except KeyError:
                    pass            
        else:
            try:
                ID = input.qualifiers['transcript_id'][0]
            except KeyError:
                pass
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        else:
            try:
                ID = input.qualifiers['standard_name'][0]
            except KeyError:
                pass        
        return locusTag, ID, Parent
                   
    elif type == 'CDS':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
            if locusTag:
                Parent = locusTag
                try:
                    ID = input.qualifiers['protein_id'][0]
                except KeyError:
                    pass
            else:
                try:
                    locusTag = input.qualifiers['protein_id'][0]
                    Parent = locusTag
                except KeyError:
                    pass       
        else:
            try:
                ID = input.qualifiers['protein_id'][0]
            except KeyError:
                ID = locusTag
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        else:
            try:
                ID = input.qualifiers['standard_name'][0]
            except KeyError:
                pass     
        return locusTag, ID, Parent

def gb_feature_add2dict(f, record, genes):
    '''
    general function to take a genbank feature from flat file and add to funannotate standardized dictionary
    locustag: {
    'contig': contigName
    'type': mRNA/rRNA/tRNA/ncRNA
    'location': (start, end) #integer tuple
    'strand': +/-
    'ids': [transcript/protein IDs] #list
    'mRNA':[[(ex1,ex1),(ex2,ex2)]] #list of lists of tuples (start, end)
    'CDS':[[(cds1,cds1),(cds2,cds2)]] #list of lists of tuples (start, end)
    'transcript': [seq1, seq2] #list of mRNA trnascripts
    'cds_transcript': [seq1, seq2] list of mRNA (no UTRs)
    'protein': [protseq1,protseq2] #list of CDS translations
    'protein_id': [id,id] #from NCBI 
    'codon_start': [1,1] #codon start for translations
    'note': [[first note, second note], [first, second, etc]] #list of lists
    'name': genename
    'product': [hypothetical protein, velvet complex] #list of product definitions
    'go_terms': [[GO:0000001,GO:0000002]] #list of lists
    'db_xref': [[InterPro:IPR0001,PFAM:004384]] #list of lists
    'partialStart': True/False
    'partialStop': True/False
    'source': annotation source
    'pseudo': True/False
    }
    '''
    #get info from features, if there is no locusTag then exit
    if f.type == 'gene' or f.type == 'mRNA' or f.type == 'CDS' or f.type == 'tRNA' or f.type == 'rRNA' or f.type == 'ncRNA' or f.type == 'exon':
        locusTag, ID, Parent = getID(f, f.type)
        if not locusTag:
            return genes
    else:
        return genes
    #check for mismatching funannotate ID locus tag basename
    if ID and '-T' in ID: #then this is from funannotate, okay to modify - this is to capture apparent tbl2asn local error
        if ID.split('-T')[0] != locusTag:  #there is a problem, update locusTag with basename of ID
            locusTag = ID.split('-T')[0]
    #standard information from every feature
    strand = f.location.strand
    if strand == 1:
        strand = '+'
    elif strand == -1:
        strand = '-'
    start = f.location.nofuzzy_start + 1
    end = f.location.nofuzzy_end
    chr = record.id
    num_parts = len(f.location.parts)
    name,Product = (None,)*2
    Fivepartial,Threepartial = (False,)*2
    DBxref = []
    Note = []
    GO = []
    inference = []
    pseudo = False
    if 'pseudo' in f.qualifiers:
        pseudo = True
    #parse each type somewhat differently
    if f.type == 'gene': 
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': None, 'transcript': [], 'cds_transcript': [], 'protein': [], 'source': 'GenBank',
            'codon_start': [], 'ids': [], 'CDS': [], 'mRNA': [], 'strand': strand, 
            'location': (int(start), int(end)), 'contig': chr, 'product': [], 'inference': [],
            'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [], 'partialStop': [], 'protein_id': [], 'pseudo': pseudo}
        else:
            genes[locusTag]['location'] = (int(start), int(end))
            genes[locusTag]['strand'] = strand
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name          
    elif f.type == 'tRNA' or f.type == 'rRNA' or f.type == 'ncRNA':
        feature_seq = str(f.extract(record.seq))
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        try:
            Product = f.qualifiers['product'][0]
            if Product == 'tRNA-OTHER':
                Product = 'tRNA-Xxx'
        except KeyError:
            Product = None
        exonTuples = []
        if num_parts < 2: #only single exon
            exonTuples.append((int(start),int(end)))
        else: #more than 1 exon, so loop through
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                exonTuples.append((int(ex_start),int(ex_end)))
        #now we want to sort the positions I think...
        if strand == '+':
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
            if unicode(f.location.start).startswith('<'):
                Fivepartial = True
            if unicode(f.location.end).startswith('>'):
                Threepartial = True
        else:
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0], reverse=True)
            if unicode(f.location.start).startswith('<'):
                Threepartial = True
            if unicode(f.location.end).startswith('>'):
                Fivepartial = True  
        #update positions
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': f.type, 'transcript': [feature_seq], 'cds_transcript': [], 'protein': [], 'source': 'GenBank',
            'codon_start': [], 'ids': [locusTag+'-T1'], 'CDS': [], 'mRNA': [sortedExons], 'strand': strand, 'inference': [],
            'location': (int(start), int(end)), 'contig': chr, 'product': [Product], 'protein_id': [], 'pseudo': pseudo,
            'db_xref': [DBxref], 'go_terms': [GO], 'note': [Note], 'partialStart': [Fivepartial], 'partialStop': [Threepartial]}
        else:
            genes[locusTag]['mRNA'].append(sortedExons)
            genes[locusTag]['type'] = f.type
            genes[locusTag]['transcript'].append(feature_seq)
            genes[locusTag]['ids'].append(locusTag+'-T'+str(len(genes[locusTag]['ids'])+1))
            genes[locusTag]['db_xref'].append(DBxref)
            genes[locusTag]['note'].append(Note)
            genes[locusTag]['go_terms'].append(GO)
            genes[locusTag]['product'].append(Product)
            genes[locusTag]['partialStart'].append(Fivepartial)
            genes[locusTag]['partialStop'].append(Threepartial)
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name          
    elif f.type == 'mRNA':
        feature_seq = str(f.extract(record.seq))
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        exonTuples = []
        if num_parts < 2: #only single exon
            exonTuples.append((int(start),int(end)))
        else: #more than 1 exon, so loop through
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                exonTuples.append((int(ex_start),int(ex_end)))
        #now we want to sort the positions I think...
        if strand == '+':
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
            if unicode(f.location.start).startswith('<'):
                Fivepartial = True
            if unicode(f.location.end).startswith('>'):
                Threepartial = True
        else:
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0], reverse=True)
            if unicode(f.location.start).startswith('<'):
                Threepartial = True
            if unicode(f.location.end).startswith('>'):
                Fivepartial = True  
        #update positions
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': f.type, 'transcript': [feature_seq], 'cds_transcript': [], 'protein': [], 'source': 'GenBank',
            'codon_start': [], 'ids': [], 'CDS': [], 'mRNA': [sortedExons], 'strand': strand, 'inference': [],
            'location': (int(start), int(end)), 'contig': chr, 'product': [], 'protein_id': [], 'pseudo': pseudo,
            'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [Fivepartial], 'partialStop': [Threepartial]}
        else:
            genes[locusTag]['mRNA'].append(sortedExons)
            genes[locusTag]['type'] = f.type
            genes[locusTag]['transcript'].append(feature_seq)
            genes[locusTag]['partialStart'].append(Fivepartial)
            genes[locusTag]['partialStop'].append(Threepartial)
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    elif f.type == 'exon': #assuming need to overwrite mRNA feature then?
        genes[locusTag]['mRNA'] = []
        genes[locusTag]['transcript'] = []
        feature_seq = str(f.extract(record.seq))
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        exonTuples = []
        if num_parts < 2: #only single exon
            exonTuples.append((int(start),int(end)))
        else: #more than 1 exon, so loop through
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                exonTuples.append((int(ex_start),int(ex_end)))
        #now we want to sort the positions I think...
        if strand == '+':
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
            if unicode(f.location.start).startswith('<'):
                Fivepartial = True
            if unicode(f.location.end).startswith('>'):
                Threepartial = True
        else:
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0], reverse=True)
            if unicode(f.location.start).startswith('<'):
                Threepartial = True
            if unicode(f.location.end).startswith('>'):
                Fivepartial = True  
        #update positions
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': f.type, 'transcript': [feature_seq], 'cds_transcript': [], 'protein': [], 'source': 'GenBank',
            'codon_start': [], 'ids': [], 'CDS': [], 'mRNA': [sortedExons], 'strand': strand, 'inference': [],
            'location': (int(start), int(end)), 'contig': chr, 'product': [], 'protein_id': [], 
            'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [Fivepartial], 'partialStop': [Threepartial], 'pseudo': pseudo}
        else:
            genes[locusTag]['mRNA'].append(sortedExons)
            genes[locusTag]['transcript'].append(feature_seq)
            genes[locusTag]['partialStart'].append(Fivepartial)
            genes[locusTag]['partialStop'].append(Threepartial)             
    elif f.type == 'CDS' and 'codon_start' in f.qualifiers:
        feature_seq = str(f.extract(record.seq))
        if not ID:
            try:
                log.info("putative transcript from %s has no ID\n(%s %s %s)" % (locusTag, locusTag, ID, Parent))
            except NameError:
                print("putative transcript from %s has no ID\n(%s %s %s)" % (locusTag, locusTag, ID, Parent))
            return genes
        try:
            protSeq = f.qualifiers['translation'][0]
        except KeyError:
            try:
                log.debug("%s has no translation" % ID)
            except NameError:
                print("%s has no translation" % ID)
            protSeq = ''
        cdsTuples = []
        phase = int(f.qualifiers['codon_start'][0])
        if num_parts < 2: #only single CDS
            cdsTuples.append((int(start),int(end)))
        else:
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1 
                ex_end = f.location.parts[i].nofuzzy_end
                cdsTuples.append((int(ex_start),int(ex_end)))
        if strand == '+':
            sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0])
        else:
            sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0], reverse=True)     
        #check for annotations
        try:
            Product = f.qualifiers['product'][0]
        except KeyError:
            Product = 'hypothetical protein'
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        #note and dbxref are in a dictionary
        for key,value in f.qualifiers.items():
            if key == 'note':
                notes = value[0].split('; ')
                for n in notes:
                    if n.startswith('GO'):
                        GO.append(n)
                    else:
                        Note.append(n)
            elif key == 'db_xref':
                for ref in value:
                    DBxref.append(ref)
            elif key == 'inference':
            	for x in value:
            		inference.append(x)                    
        #update dictionary
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': 'mRNA', 'transcript': [], 'cds_transcript': [feature_seq], 'protein': [], 'source': 'GenBank',
            'codon_start': [phase], 'ids': [locusTag+'-T1'], 'CDS': [sortedCDS], 'mRNA': [], 'strand': strand, 'inference': inference,
            'location': (int(start), int(end)), 'contig': chr, 'product': [Product], 'protein_id': [ID],
            'db_xref': [DBxref], 'go_terms': [GO], 'note': [Note], 'partialStart': [], 'partialStop': [], 'pseudo': pseudo}
        else:
            genes[locusTag]['protein_id'].append(ID)
            genes[locusTag]['ids'].append(locusTag+'-T'+str(len(genes[locusTag]['ids'])+1))
            genes[locusTag]['CDS'].append(sortedCDS)
            genes[locusTag]['product'].append(Product)
            genes[locusTag]['protein'].append(protSeq)
            genes[locusTag]['cds_transcript'].append(feature_seq)
            genes[locusTag]['codon_start'].append(phase)
            genes[locusTag]['db_xref'].append(DBxref)
            genes[locusTag]['note'].append(Note)
            genes[locusTag]['go_terms'].append(GO)
            genes[locusTag]['inference'] += inference
            if not genes[locusTag]['type']:
            	genes[locusTag]['type'] = 'mRNA'
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    return genes

#load genbank file into dictionary, while looping through also capture antismash cluster locations
Genes = {}
antiSMASH = {}
with open(sys.argv[1], 'rU') as infile:
	for record in SeqIO.parse(infile, 'genbank'):
		for f in record.features:
			if f.type == 'cluster':
				number = len(antiSMASH)+1
				chr = record.id
				start = f.location.start
				end = f.location.end
				antiSMASH[number] = (chr, start, end)
			else: #add to database
				Genes = gb_feature_add2dict(f, record, Genes)

#create interlap object with clusters
Clusters = defaultdict(InterLap)
#now populate interlap object
for k,v in antiSMASH.items():
	contig, start, end = v
	Clusters[contig].add((int(start), int(end), 'cluster{:}'.format(k)))

SeqRecords = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'genbank'))

sys.stdout.write('#GeneID\ttranscriptID\tname\tproduct\ttype\tcontig\tstart\tend\tstrand\tantiSMASH\tinference\tDBxref\tGO\tNote\tgDNA\tmRNA\tProtein\n')
for k,v in natsorted(Genes.items()):
	clusterAnnot = ''
	if v['location'] in Clusters[v['contig']]:
		clusterAnnot = list(Clusters[v['contig']].find(v['location']))[0][2]
	#output to stdout the tsv file
	for i in range(0, len(v['ids'])):
		if not v['name']:
			name = ''
		else:
			name = v['name']
		if v['type'] == 'mRNA':
			try:
				prot = v['protein'][i]
			except IndexError:
				prot = ''
			sys.stdout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
				k, v['ids'][i], v['type'], name, v['product'][i], v['contig'], v['location'][0], v['location'][1], v['strand'], clusterAnnot,
				';'.join(v['inference']), ';'.join(v['db_xref'][i]), ';'.join(v['go_terms'][i]),
				';'.join(v['note'][i]), str(SeqRecords[v['contig']][v['location'][0]:v['location'][1]].seq),
				v['cds_transcript'][i], prot))
		else:
			sys.stdout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
				k, v['ids'][i], v['type'], name, v['product'][i], v['contig'], v['location'][0], v['location'][1], v['strand'], clusterAnnot,
				';'.join(v['inference']), ';'.join(v['db_xref'][i]), ';'.join(v['go_terms'][i]),
				';'.join(v['note'][i]), str(SeqRecords[v['contig']][v['location'][0]:v['location'][1]].seq),
				v['transcript'][i], ''))