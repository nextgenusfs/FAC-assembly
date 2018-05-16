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
```
This would be run like so:
```
$ assemble_facs.py -v pFACvector.fasta -1 C1_S9_L001_R1_001.fastq.gz -2 C1_S9_L001_R2_001.fastq.gz
```
The output will be an assembly file called `C1_unicycler.fasta`.


### Dependencies
* Python 2
* trim_galore
* minimap2
* unicycler
* samtools
* mappy (minimap2 python library


