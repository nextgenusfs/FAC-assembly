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
#Organism	TeloFAC_Order 	FAC_Name	Sanger_fwd_seq	Sanger_rev_seq	Column_Pool	Row_Pool	Plate_Pool	Illumina_seqpools
Aspergillus aculeatus	1A1	IGAaBAC1D1	TAGGGCCGGTCGCGGGGGTCCCTG	TTCAATAGCCCTATGGGTCGTGCG	C1	RA	P1	"C1, RA, P1"
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
  --min_len_complete MIN_LEN_COMPLETE         Name of seq pool (default: 75000)
```
It can be run using the above example like this:
```
$ get_full_length_facs.py -a C1_unicycler.fasta -i IGteloFACpooldataNF.txt -p C1
```

### Dependencies
* Python 2
* trim_galore
* minimap2
* unicycler
* samtools
* mappy (minimap2 python library


