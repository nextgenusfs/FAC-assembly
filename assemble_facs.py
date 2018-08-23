#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
import shutil
import datetime
import platform

__version__ = "0.1.0"
#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='assemble_facs.py',
    description = '''Script to quality trim, subtract vector reads, and assemble using unicyler or shovill''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-v', '--vector', required=True, help='FASTA vector sequence')
parser.add_argument('-1', '--forward', required=True, help='FASTQ R1')
parser.add_argument('-2', '--reverse', required=True, help='FASTQ R2')
parser.add_argument('-c', '--cpus', default=8, type=int, help='Number of CPUS')
parser.add_argument('-a', '--assembler', default='unicycler', choices=['unicycler', 'shovill'],help='Assembler to use')
parser.add_argument('--shovill_assembler', default='skesa', choices=['spades', 'skesa', 'megahit', 'velvet'], help='Shovill assembler to use')
parser.add_argument('--shovill_ram', default=16, type=int, help='RAM in GB')
parser.add_argument('--skip_quality', action='store_true', help='Skip quality/adapter trimming')
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

def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def fasta2bed(input, output):
    with open(output, 'w') as outfile:
        with open(input, 'rU') as infile:
            for Header, Seq in SimpleFastaParser(infile):
                outfile.write('{:}\t{:}\t{:}\n'.format(Header.split(' ')[0], 0, len(Seq)))

def unicycler_version():
    p = subprocess.Popen(['unicycler', '--version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
    vers = p.split(' v')[-1].replace('b', '')
    if tuple(vers.split('.')) > ('0', '4', '0'):
        return True
    else:
        return False
            
#get basename
if '_' in os.path.basename(args.forward):
    base = os.path.basename(args.forward).split('_')[0]
else:
    base = os.path.basename(args.forward).split('.',-1)[0]

#get compression threads
bamcpus = args.cpus // 2

#generate logfile
logfile = base + '_assemble_facs.log'
if os.path.isfile(logfile):
    os.remove(logfile)

print('------------------------------------------')
print('[{:}] Running Python v{:} '.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), platform.python_version()))
dependencies = ['trim_galore', 'minimap2', 'samtools', 'unicycler']
for x in dependencies:
    if not which_path(x):
        print('{:} is not properly installed, install and re-run script'.format(x))
        sys.exit(1)
        
#check input files
files = [args.vector, args.forward, args.reverse]
for x in files:
    if not os.path.isfile(os.path.abspath(x)):
        print('{:} file is not found, full path: {:}'.format(x, os.path.abspath(x)))
        sys.exit(1) 
        
with open(logfile, 'w') as log:
    #quality/adapter filter
    if not args.skip_quality:
        qualR1 = os.path.basename(args.forward).split('.f')[0]
        qualR1 = qualR1+'_val_1.fq.gz'
        qualR2 = os.path.basename(args.reverse).split('.f')[0]
        qualR2 = qualR2+'_val_2.fq.gz'
        if not os.path.isfile(qualR1) and not os.path.isfile(qualR2):
            print('[{:}] Trimming adapters from reads using Trim Galore'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
            cmd = ['trim_galore', '--paired', '--length', '100', '--quality', '10', os.path.abspath(args.forward), os.path.abspath(args.reverse)]
            subprocess.call(cmd, stderr=log, stdout=log)
        else:
            print('[{:}] Using existing Trim-galore output'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))

    else:
        qualR1 = args.forward
        qualR2 = args.reverse

    #map the reads to the vector using minimap2, keep unmapped reads
    print('[{:}] Mapping reads to vector using minimap2'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
    cmd = ['minimap2', '-ax', 'sr', '-t', str(args.cpus), os.path.abspath(args.vector), qualR1, qualR2]
    cleanR1 = base+'_clean_R1.fastq.gz'
    cleanR2 = base+'_clean_R2.fastq.gz'
    if not os.path.isfile(cleanR1) and not os.path.isfile(cleanR2):
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log)
        p2 = subprocess.Popen(['samtools', 'sort', '-@', str(bamcpus), '-'], stdin=p1.stdout, stderr=log, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(['samtools', 'fastq', '-1', cleanR1, '-2', cleanR2, '-@', str(bamcpus), '-f', '12', '-'], stdin=p2.stdout, stderr=log, stdout=subprocess.PIPE)
        p1.stdout.close()
        p2.stdout.close()
        p3.communicate()
    else:
        print('[{:}] Using existing vector-cleaned FASTQ files'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
    
    assembly = base+'.assembly.fasta'
    if args.assembler == 'unicycler':
        #now run unicycler
        print('[{:}] Assembling with Unicycler'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
        if unicycler_version():
            cmd = ['unicycler', '-t', str(args.cpus), '--min_fasta_length', '1000', '--linear_seqs', '60', '-1', cleanR1, '-2', cleanR2, '-o', base+'_unicycler', '--no_rotate']
        else:
            cmd = ['unicycler', '-t', str(args.cpus), '--min_fasta_length', '1000', '-1', cleanR1, '-2', cleanR2, '-o', base+'_unicycler', '--no_rotate']
        subprocess.call(cmd)
        with open(assembly, 'w') as outfile:
            with open(os.path.join(base+'_unicycler', 'assembly.fasta'), 'rU') as infile:
                for line in infile:
                    if line.startswith('>'):
                        line = line.replace('>', '>ctg')
                    outfile.write(line)
        print('[{:}] Unicycler finished, final output assembly: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), assembly))
    elif args.assembler == 'shovill':
        print('[{:}] Assembling with Shovill - using {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), args.shovill_assembler))
        cmd = ['shovill', '--ram', str(args.shovill_ram), '--assembler', args.shovill_assembler, '--outdir', 'shovill_out', '--R1', cleanR1, '--R2', cleanR2, '--cpus', str(args.cpus)]
        if os.path.isdir('shovill_out'):
            cmd.append('--force')
        subprocess.call(cmd)
        with open(assembly, 'w') as outfile:
            with open(os.path.join('shovill_out', 'contigs.fa'), 'rU') as infile:
                for line in infile:
                    if line.startswith('>'):
                        line = line.replace('>contig', '>ctg')
                    outfile.write(line)
        print('[{:}] Shovill finished, final output assembly: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), assembly))
    
    #now map reads back to assembly and get coverage stats
    print('[{:}] Mapping reads to assembly, calculating coverage'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
    lenBED = base+'.assembly.bed'
    fasta2bed(os.path.abspath(assembly), lenBED)
    cmd = ['minimap2', '-ax', 'sr', '-t', str(args.cpus), os.path.abspath(assembly), cleanR1, cleanR2]
    mappedBAM = base+'.mapped2assembly.bam'
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log)
    p2 = subprocess.Popen(['samtools', 'sort', '-@', str(bamcpus), '-o', mappedBAM, '-'], stdout=subprocess.PIPE, stderr=log, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()
    subprocess.call(['samtools', 'index', mappedBAM])
    coverageBed = base+'.mapped_coverage.bed'
    coverageFinal = base+'.coverage-stats.txt'
    cov_cmd = ['samtools', 'bedcov', lenBED, mappedBAM]
    if os.path.isfile(mappedBAM) and os.path.isfile(lenBED): #check files exist before running
        with open(coverageBed, 'w') as output:
            subprocess.call(cov_cmd, stdout=output, stderr=log)
        if os.path.isfile(coverageBed):
            with open(coverageFinal, 'w') as outfile:
                with open(coverageBed, 'rU') as bedfile:
                    for line in bedfile:
                        log.write(line)
                        if not line or line.startswith('\n') or line.count('\t') < 3:
                            continue
                        line = line.strip()
                        cols = line.split('\t')
                        cov = int(cols[3]) / float(cols[2])
                        print('{:} len={:} coverage={:.2f}X'.format(cols[0], cols[2], cov))
                        outfile.write('{:} len={:} coverage={:.2f}X\n'.format(cols[0], cols[2], cov))
    print('[{:}] Logfile: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), logfile))
print('------------------------------------------')
