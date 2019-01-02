# FAC-assembly
Scripts for assembly of PE Illumina from pooled BAC/FACs 

### Usage:
To run the assembly use the `assemble_facs.py` script:
```
usage: assemble_facs.py [-h] -v VECTOR -1 FORWARD -2 REVERSE [-c CPUS]

Script to quality trim, subtract vector reads, and assemble using unicyler

optional arguments:
  -h, --help                     show this help message and exit
  -v VECTOR, --vector VECTOR     FASTA vector sequence (default: None)
  -1 FORWARD, --forward FORWARD  FASTQ R1 (default: None)
  -2 REVERSE, --reverse REVERSE  FASTQ R2 (default: None)
  -c CPUS, --cpus CPUS           Number of CPUS (default: 8)
  --skip_quality                 Skip quality/adapter trimming (default: False)
```
This would be run like so:
```
$ assemble_facs.py -v pFACvector.fasta -1 C1_S9_L001_R1_001.fastq.gz -2 C1_S9_L001_R2_001.fastq.gz
```
The output will be an assembly file called `C1_unicycler.fasta`.


### Mapping Sanger End Sequences to assembly:

After generating assembly for each pooled run, using a tab-delimited mapping file, can pull out full-length FACs.

The tab-delimited mapping file looks like this (The sanger reads are truncated in this example):
```
#Organism   TeloFAC_Order   FAC_Name    Sanger_fwd_seq  Sanger_rev_seq  Column_Pool Row_Pool    Plate_Pool  Illumina_seqpools
Aspergillus aculeatus   1A1 IGAaBAC1D1  TAGGGCCGGTCGCGGGGGTCCCTG    TTCAATAGCCCTATGGGTCGTGCG    C1  RA  P1  "C1, RA, P1"
```
This data can then be used to map the Sanger reads to the assembly using the `get_full_length_facs.py` script.
```
usage: get_full_length_facs.py [-h] -a ASSEMBLY -i INPUT [-p POOL [POOL ...]]
                               [--min_len_complete MIN_LEN_COMPLETE]

Script that parses tsv file, maps end sequences to FAC assembly.

optional arguments:
  -h, --help                                  show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY            Assembly in FASTA format (default: None)
  -i INPUT, --input INPUT                     TSV FAC end sequences (default: None)
  -p POOL [POOL ...], --pool POOL [POOL ...]  Name of seq pool (default: None)
  --min_len MIN_LEN                           Minimum contig length to find end primers (default: 20000)
```
It can be run using the above example like this:
```
$ get_full_length_facs.py -a C1_unicycler.fasta -i IGteloFACpooldataNF.txt -p C1
```
This will generate 3 files:
1) Raw Mapping data from minimap2
2) FASTA file containing full length sequences
3) FASTA file containing all sequences

The raw data looks like:
```
FAC_name	read_num	query_len	query_start	query_end	query_cov	strand	contig	contig_len	ref_start	ref_end	match_len	aln_len	pident	map_qual	cigar
IGAwBAC1B10	1	1080	183	1078	82.9	-1	ctg00017	90870	89987	90870	875	895	97.76536312849163	60	9M1I20M1I5M1I5M1I12M1I22M1I8M1I19M1I65M1I7M1I71M1I23M1I617M
IGAwBAC1B9	1	985	60	983	93.7	1	ctg00007	130612	0	908	901	928	97.09051724137932	60	4M1I102M5D605M1I50M1I10M1I9M1I13M1I19M1I20M1I8M2I3M1I6M1I8M1I3M1I5M1I6M2I5M1I6M1I11M1I10M
IGAwBAC1B9	2	200	0	200	100.0	-1	ctg00007	130612	99649	99848	199	200	99.5	60	167M1I32M
IGAwBAC1C16	1	160	20	160	87.5	1	ctg00020	49887	23726	23866	140	140	100.0	60	140M
IGAwBAC1D11	1	200	15	200	92.5	1	ctg00020	49887	7514	7699	185	185	100.0	60	185M
IGAwBAC1D15	1	200	19	200	90.5	-1	ctg00008	129517	105945	106126	178	181	98.34254143646409	60	181M
IGAwBAC1D15	2	852	0	732	85.9	1	ctg00008	129517	36	755	705	732	96.31147540983606	60	585M1I21M1I34M1I9M1I5M1I17M1I9M1I11M1I5M1I7M2I10M2I6M
IGAwBAC1
```
And the "classified" FASTA file looks like this (if multiple FACS are overlapping separated by '|'):
```
>>IGAwBAC1B9;organism=Aspergillus wentii;origin=ctg00007;length=99848;full_length=yes;
AGCCCTAAGCCAAGCCCTAAGCCCTAAGCCCTAAGCCCTAAGCCCTAAGCCCTAAGCCCTAA...               
>IGFaBAC4L14;organism=Fusarium solani;origin=ctg00038-ctg00011;length=146644;full_length=stitched;
GGATTCCGGAATGTGGAACGGGGCACAGGCACCGACCCTATCAAGGGGTAAGGGGATTCGGTCA...
>IGPmBAC1A7;organism=Penicillium marneffei;origin=ctg00013;length=109916;full_length=no;
AAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAACCCTAAAC...
```

#### Running the data through a 2-pass method

Since the pools are split over several samples, I suggest to run each pool separately first, then pull out full-length sequences for each one of those pools. You can then run a second or third round of assembly by subtracting the reads that map to the full-length sequences.

A simple example:
```
#assemble C1, RA, and P1 separately
$ assemble_facs.py -v pFACvector.fasta -1 C1_S9_L001_R1_001.fastq.gz -2 C1_S9_L001_R2_001.fastq.gz
$ assemble_facs.py -v pFACvector.fasta -1 RA_S1_L001_R1_001.fastq.gz -2 RA_S1_L001_R2_001.fastq.gz 
$ assemble_facs.py -v pFACvector.fasta -1 P1_S21_L001_R1_001.fastq.gz -2 P1_S21_L001_R2_001.fastq.gz

#now pull out full length for each
$ get_full_length_facs.py -a C1_unicycler.fasta -i IGteloFACpooldataNF.txt -p C1
$ get_full_length_facs.py -a RA_unicycler.fasta -i IGteloFACpooldataNF.txt -p RA
$ get_full_length_facs.py -a P1_unicycler.fasta -i IGteloFACpooldataNF.txt -p P1

#combine all full-length sequences
$ cat *_complete_facs.fasta > full-length_sequences.fasta

#combine the quality trimmed and vector subtracted data for each forward/reverse read
$ cat C1_clean_R1.fastq.gz RA_clean_R1.fastq.gz P1_clean_R1.fastq.gz > cleaned_R1.fastq.gz
$ cat C1_clean_R2.fastq.gz RA_clean_R2.fastq.gz P1_clean_R2.fastq.gz > cleaned_R2.fastq.gz

#we can now re-run assembly but subtract reads that are already part of full-length sequences
$ assemble_facs.py -v full-length_sequences.fasta -1 cleaned_R1.fastq.gz -2 cleaned_R2.fastq.gz

#and then finally pull out full-length from this combined assembly
$ get_full_length_facs.py -a q-trimmed_unicycler.fasta -i IGteloFACpooldataNF.txt --pool C1 RA P1
```

There will still be some ~100 kb inserts that don't seem to match Sanger end sequences - they can be searched against the reference genome to figure out where they came from.  Then there is also an accessory script that you can use to filter the FASTA files, which can be run like this:
```
#to ouput all sequences > 80kb
filter_seqs.py -i P4.all-contigs.fasta -l 80000

#to ouput all sequences > 80kb that are Aspergillus wentii
filter_seqs.py -i P4.all-contigs.fasta -m 80000 -s "Aspergillus wentii"

#to ouput the longest unique contigs
filter_seqs.py -i P4.all-contigs.fasta -l

#to ouput the longest unique contigs from Fusarium solani
filter_seqs.py -i P4.all-contigs.fasta -l -s "Fusarium solani"
```

### Dependencies
* Python 2
* biopython
* trim_galore
* minimap2
* unicycler
* samtools
* mappy (minimap2 python library


