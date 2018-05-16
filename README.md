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
  --min_len_complete MIN_LEN_COMPLETE         Name of seq pool (default: 75000)
```
It can be run using the above example like this:
```
$ get_full_length_facs.py -a C1_unicycler.fasta -i IGteloFACpooldataNF.txt -p C1
```
This will generate 3 files:
1) Raw Mapping data from minimap2
2) FASTA file containing full length sequences
3) FASTA file containing unassigned sequences

The raw data looks like:
```
#FAC_name   read_num    query_start query_end   strand  contig  contig_len  ref_start   ref_end match_len   aln_len map_qual    cigar
IGAaBAC10A5 1   100 1005    -1  ctg13   103104  102206  103104  896 905 60  15M1I25M1I24M1I55M1I4M1I12M1I59M1I704M
IGAaBAC10A5 2   0   690 1   ctg13   103104  0   690 689 690 60  690M
AtBAC07B20  1   61  470 -1  ctg8    115790  115382  115790  404 406 60  49M1I359M
AtBAC07B20  2   3   567 1   ctg8    115790  9   574 547 553 60  5M1D559M
AtBAC06G23  1   70  381 1   ctg17   84458   1341    1652    309 311 60  311M
AtBAC06G23  2   0   416 -1  ctg19   46267   45848   46264   411 412 60  416M
```
And the full-length FASTA file looks like this (if multiple FACS are overlapping separated by '|'):
```
>IGPeBAC4N19;organism=Penicillium expansum;length=125086;depth=0.47x;
AGCCCTAAGCCAAGCCCTAAGCCCTAAGCCCTAAGCCCTAAGCCCTAAGCCCTAAGCCCTAA...               
>AtBAC10C18;organism=Aspergillus terreus;length=122271;depth=1.00x;
GGATTCCGGAATGTGGAACGGGGCACAGGCACCGACCCTATCAAGGGGTAAGGGGATTCGGTCA...
>AtBAC01J18|AtBAC07B20;organism=Aspergillus terreus;length=115790;depth=2.33x;
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

#combine the quality trimmed data for each forward/reverse read
$ cat C1_S9_L001_R1_001_val_1.fq.gz RA_S1_L001_R1_001_val_1.fq.gz P1_S21_L001_R1_001_val_1.fq.gz > q-trimmed_R1.fastq.gz
$ cat C1_S9_L001_R2_001_val_2.fq.gz RA_S1_L001_R2_001_val_2.fq.gz P1_S21_L001_R2_001_val_2.fq.gz > q-trimmed_R2.fastq.gz

#we can now re-run assembly but subtract reads that are already part of full-length sequences
$ assemble_facs.py -v full-length_sequences.fasta -1 q-trimmed_R1.fastq.gz -2 q-trimmed_R2.fastq.gz --skip_quality

#and then finally pull out full-length from this combined assembly
$ get_full_length_facs.py -a q-trimmed_unicycler.fasta -i IGteloFACpooldataNF.txt --pool C1 RA P1
```

There will still be some ~100 kb inserts that don't seem to match Sanger end sequences - they can be searched against the reference genome to figure out where they came from.

### Dependencies
* Python 2
* trim_galore
* minimap2
* unicycler
* samtools
* mappy (minimap2 python library


