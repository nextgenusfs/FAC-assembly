#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
import shutil
import datetime
import platform
import socket
import gzip
import io
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

__version__ = "0.2.0"
#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='assemble_facs.py',
    description = '''Script to quality trim, subtract vector reads, and assemble using unicyler or shovill''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-v', '--vector', help='FASTA vector sequence')
parser.add_argument('-s', '--subtract', default='CP014272.1', help='FASTA sequence to map and subtract from reads')
parser.add_argument('-1', '--forward', required=True, help='FASTQ R1')
parser.add_argument('-2', '--reverse', required=True, help='FASTQ R2')
parser.add_argument('-c', '--cpus', default=8, type=int, help='Number of CPUS')
parser.add_argument('-a', '--assembler', default='shovill', choices=['unicycler', 'shovill'],help='Assembler to use')
parser.add_argument('--shovill_assembler', default='spades', choices=['spades', 'skesa', 'megahit', 'velvet'], help='Shovill assembler to use')
parser.add_argument('--shovill_ram', default=16, type=int, help='RAM in GB')
parser.add_argument('--skip_quality', action='store_true', help='Skip quality/adapter trimming')
parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
parser.add_argument('--quiet', action='store_true', help='Keep messages in terminal minimal')
args=parser.parse_args()

def download(url, name):
    file_name = name
    try:
        u = urlopen(url)
        f = open(file_name, 'wb')
        block_sz = 8192
        while True:
            buffer = u.read(block_sz)
            if not buffer:
                break
            f.write(buffer)
        f.close()
    except socket.error as e:
        if e.errno != errno.ECONNRESET:
            raise
        pass

def countfastq(input):
    if input.endswith('.gz'):
        linenum = 0
        gz = gzip.open(input, 'rb')
        f = io.BufferedReader(gz)
        for line in f.readlines():
            linenum += 1
        gz.close()
        count = linenum // 4
    else:
        lines = sum(1 for line in open(input))
        count = int(lines) // 4
    return count

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
    p = subprocess.Popen(['unicycler', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    p = p.decode(encoding='utf-8')
    vers = p.split(' v')[-1].replace('\n', '')
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
if bamcpus > 4:
	bamcpus = 4

#generate logfile
logfile = base + '_assemble_facs.log'
if os.path.isfile(logfile):
    os.remove(logfile)

print('------------------------------------------')
print('[{:}] Running Python v{:} '.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), platform.python_version()))
dependencies = ['minimap2', 'samtools']
dependencies.append(args.assembler)
for x in dependencies:
    if not which_path(x):
        print('{:} is not properly installed, install and re-run script'.format(x))
        sys.exit(1)
        
#check input files
files = [args.forward, args.reverse]
for x in files:
    if not os.path.isfile(os.path.abspath(x)):
        print('{:} file is not found, full path: {:}'.format(x, os.path.abspath(x)))
        sys.exit(1) 

origCount = countfastq(args.forward)
print('[{:}] Loaded {:,} PE reads'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), origCount))
with open(logfile, 'w') as log:
    subtractDB = base+"."+str(os.getpid())+'.subtract.fasta'
    acc_file = args.subtract+".fna"
    if not os.path.isfile(acc_file):
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta' % (args.subtract)
        download(url, acc_file)
    if args.vector and os.path.isfile(acc_file):
        with open(subtractDB, 'wb') as wfd:
            for fname in [args.vector, acc_file]:
                with open(fname,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
    else:
        if not args.vector:
            shutil.copyfile(acc_file, subtractDB)

    #map the reads to the vector using minimap2, keep unmapped reads
    print('[{:}] Mapping reads to subtraction database using minimap2'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
    cmd = ['minimap2', '-ax', 'sr', '-t', str(args.cpus), os.path.abspath(subtractDB), args.forward, args.reverse]
    cleanR1 = base+'_clean_R1.fastq.gz'
    cleanR2 = base+'_clean_R2.fastq.gz'
    if not os.path.isfile(cleanR1) and not os.path.isfile(cleanR2):
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log)
        p2 = subprocess.Popen(['samtools', 'sort', '-@', str(bamcpus), '-'], stdin=p1.stdout, stderr=log, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(['samtools', 'fastq', '-1', cleanR1, '-2', cleanR2, '-@', str(bamcpus), '-f', '12', '-'], stdin=p2.stdout, stderr=log, stdout=subprocess.PIPE)
        p1.stdout.close()
        p2.stdout.close()
        p3.communicate()
        cleanCount = countfastq(cleanR1)
    else:
        print('[{:}] Using existing vector-cleaned FASTQ files'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
        cleanCount = countfastq(cleanR1)
    print('[{:}] {:,} PE reads pass subtraction filtering'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), cleanCount))
    
    assembly = base+'.assembly.fasta'
    if args.assembler == 'unicycler':
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
                qualCount = countfastq(qualR1)
            else:
                print('[{:}] Using existing Trim-galore output'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
                qualCount = countfastq(qualR1)
        else:
            qualR1 = cleanR1
            qualR2 = cleanR2
            qualCount = countfastq(qualR1)
        print('[{:}] {:,} PE reads passed quality trimming'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), qualCount))

        #now run unicycler
        print('[{:}] Assembling with Unicycler'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
        if unicycler_version():
            cmd = ['unicycler', '--mode', 'conservative', '-t', str(args.cpus), '--min_fasta_length', '1000', '--linear_seqs', '60', '-1', qualR1, '-2', qualR2, '-o', base+'_unicycler', '--no_rotate']
        else:
            cmd = ['unicycler', '--mode', 'conservative', '-t', str(args.cpus), '--min_fasta_length', '1000', '-1', qualR1, '-2', qualR2, '-o', base+'_unicycler', '--no_rotate']
        if args.quiet:
            subprocess.call(cmd, stderr=log, stdout=log)
        else:
            subprocess.call(cmd)
        assemblyCount = 0
        with open(assembly, 'w') as outfile:
            with open(os.path.join(base+'_unicycler', 'assembly.fasta'), 'rU') as infile:
                for line in infile:
                    if line.startswith('>'):
                        assemblyCount += 1
                        line = line.replace('>', '>ctg')
                    outfile.write(line)
        print('[{:}] Unicycler finished, {:,} contigs in assembly: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), assemblyCount, assembly))
    elif args.assembler == 'shovill':
        print('[{:}] Assembling with Shovill - using {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), args.shovill_assembler))
        cmd = ['shovill', '--ram', str(args.shovill_ram), '--minlen', '1000', '--assembler', args.shovill_assembler, '--outdir', 'shovill_out', '--R1', cleanR1, '--R2', cleanR2, '--cpus', str(args.cpus)]
        if os.path.isdir('shovill_out'):
            cmd.append('--force')
        if not args.skip_quality:
            cmd.append('--trim')
        if args.quiet:
            subprocess.call(cmd, stderr=log, stdout=log)
        else:
            subprocess.call(cmd)
        assemblyCount = 0
        with open(assembly, 'w') as outfile:
            with open(os.path.join('shovill_out', 'contigs.fa'), 'rU') as infile:
                for line in infile:
                    if line.startswith('>'):
                        assemblyCount += 1
                        line = line.replace('>contig', '>ctg')
                    outfile.write(line)
        print('[{:}] Shovill finished, {:,} contigs in assembly: {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), assemblyCount, assembly))
    
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
