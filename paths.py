"""Module with static paths to avoid circular imports"""
import os
import re
import subprocess

STATIC_PREFIX = './static/'
# batch
FEATURE_FASTA_PATH = STATIC_PREFIX + 'features.fa'
UNS_PATH = STATIC_PREFIX + 'uns.txt'
PLASMID_PREFIX = STATIC_PREFIX + 'plasmid/'
OUTPUT_PREFIX = STATIC_PREFIX + 'output/'
IDT_CSV = OUTPUT_PREFIX + 'idt_import.csv'
PRIMERS_CSV = OUTPUT_PREFIX + 'primers.csv'
PLASMIDS_CSV = OUTPUT_PREFIX + 'plasmids.csv'
PLASMIDS_ZIP = OUTPUT_PREFIX + 'plasmids.zip'
LOCI_ZIP = OUTPUT_PREFIX + 'loci.zip'
# genome
GENOME_PREFIX = STATIC_PREFIX + 'genome/'
GENOME_PATH = GENOME_PREFIX + 'ASM2700v1.fa'
GFF_PATH = GENOME_PREFIX + 'ASM2700v1_genes.gff3'
LOCUS_DICT_PKL_PATH = GENOME_PREFIX + 'locus_dict.pkl'
# tool
TOOL_PREFIX = STATIC_PREFIX + 'tool/'
MAKEBLASTDB_PATH = TOOL_PREFIX + 'ncbi-blast-2.6.0+/bin/makeblastdb'
BLASTN_PATH = TOOL_PREFIX + 'ncbi-blast-2.6.0+/bin/blastn'
PRIMER3_PATH = TOOL_PREFIX + 'primer3-2.3.7/src/primer3_core'
PRIMER3_CONFIG_PATH = TOOL_PREFIX + 'primer3-2.3.7/src/primer3_config/'
SAMTOOLS_PATH = 'samtools'
# tests
TESTDATA_PATH = STATIC_PREFIX + 'testdata/'

# Run some checks on import to fail fast
PREFIX_RE = re.compile('^[A-Z_]+PREFIX$')
if not os.path.exists(OUTPUT_PREFIX):
    subprocess.check_call(['mkdir',OUTPUT_PREFIX])
for path in dir():
    if re.match(PREFIX_RE, path) and path.endswith('PREFIX'):
        if not os.path.exists(eval(path)):
            raise EnvironmentError('{} not found. App needs this directory to run.'.format(eval(path)))
