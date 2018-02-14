"""Wrapper methods to run compiled tools (BLAST, Primer3, samtools)
and use objects from large modules (Bio)
Classes:
    ApeFile
Methods:
    blastn
    make_genome_blastdb
    primer3_pcr
    format_primer3_params
    best_pair
    passes_blastn_filter
    load_primer3_defaults
    faidx"""

from Bio import Alphabet, Blast, Seq, SeqFeature, SeqIO, SeqRecord
import colorsys
import dna
from itertools import chain, cycle, izip_longest
import os
import paths
import re
import subprocess

TOOL_PREFIX = './static/tool/'
PRIMER3_RE = re.compile(
 '^PRIMER_(?P<side>LEFT|RIGHT)_(?P<index>\d+)_(?P<attr>[^=]+)=(?P<value>.*)$'
 )

class ApeFile(SeqRecord.SeqRecord):
    def __init__(self, seq_str, id_str, max_colors=20):
        """Args:
            seq_str (str): sequence to annotate and write to file
            id_str (str): name of the record in the apefile
            max_colors (int): number of spokes in the color wheel
        Attributes: (only those used in this class)
            seq (Seq): sequence object created from input
            features (list of SeqFeature): list of annotations
            colors (list of str): hex codes of equal-lightness colors
                                  for legible highlighting"""
        super(ApeFile, self).__init__(
                Seq.Seq(seq_str, Alphabet.generic_dna), id_str)
        # Each object gets a cyclic iterator of colors,
        # starting with half-saturated green (cool colors)
        self.colors = cycle(('#{:02x}{:02x}{:02x}'.format(
                *map(lambda x: int(x*255),
                colorsys.hls_to_rgb(float(h % max_colors)/max_colors,0.5,0.5)))
                    for h in xrange(max_colors/3,max_colors*4/3)))
    def add_feature(self, search_str, label):
        """Label a feature by literal string match, failing silently.
        Does not label features that wrap around circular sequence
        Args:
            search_str (str): string representing feature
            label (str): feature name to display"""
        f_matches = re.finditer(search_str, str(self.seq), re.IGNORECASE)
        r_matches = ()
        if search_str.lower() != dna.revc(search_str.lower()):
            # Don't label palindromes twice
            r_matches = re.finditer(dna.revc(search_str),
                    str(self.seq), re.IGNORECASE)
        for m, strand in chain(izip_longest(f_matches, [1]),
                izip_longest(r_matches, [-1])):
            if m:
                feature = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(
                    m.start(), m.end(), strand), 'misc_feature')
                color = self.colors.next()
                feature.qualifiers = {
                        'label': [label],
                        'ApEinfo_fwdcolor': [color],
                        'ApEinfo_revcolor': [color],
                        'ApEinfo_graphicformat': [
                            'arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0'
                            ]}
                self.features.append(feature)
        return
    def add_from_fasta(self, fname=paths.FEATURE_FASTA_PATH):
        with open(fname, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                self.add_feature(str(record.seq), record.id)
        return
    def write(self, fname):
        SeqIO.write(self, fname, 'genbank')
        return

def blastn(query):
    """Run BLASTN on genome defined in locus module.
    Args:
        query (str): nucleotide sequence to BLAST
    Returns:
        (tuple of str): (chromosome, start, end, strand)
        sseq (str): aligned sequence in subject (reference)"""
    if len(query) < 100:
        task = 'blastn-short'
    else:
        task = 'megablast'
    cmd_list = [paths.BLASTN_PATH, '-db', paths.GENOME_PATH,
                '-task', task, '-max_hsps', '1',
                '-outfmt','6 sseqid sstart send length pident sseq']
    p = subprocess.Popen(cmd_list,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    # Output is tab-separated fields, one line per alignment
    # sorted by best scoring first
    output = p.communicate(query)[0]
    if not output:
        raise ValueError('No BLASTN hits for query.')
    chromosome, start, end, length, pident, sseq = output.partition(
            '\n')[0].split('\t')
    if int(length) < len(query)/2 or float(pident) < 50:
        raise ValueError('No BLASTN hits spanning 50% of query '
        +'or at least 50% identity in alignment.')
    if start < end:
        strand = '+'
    else:
        # Keep GFF convention of start < end on minus strand
        strand = '-'
        start, end = end, start
    return (chromosome, start, end, strand), sseq.lower()

def make_genome_blastdb():
    """Checks for existing blastdb, else runs makeblastdb on genome."""
    if not os.path.exists(paths.GENOME_PATH + '.nin'):
        cmd_list = [paths.MAKEBLASTDB_PATH, '-dbtype', 'nucl',
                '-in', paths.GENOME_PATH]
        subprocess.check_call(cmd_list)
    return

def primer3_pcr(template, size_range, blast=False, options=None):
    """Runs a primer3 task
    Args:
        template (Dna): template sequence for priming.
                        SEQUENCE_TEMPLATE required for primer search.
        size_range (list of int): (min, max) product size.
                                  PRIMER_PRODUCT_SIZE_RANGE required.
        options (dict): key:value pairs of optional primer3 parameters
            SEQUENCE_TARGET=start,length: primer pair must flank this sequence
            SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=
            left_start,left_length,right_start,right_length
        blast (bool): whether to filter for genome mispriming by BLASTN.
                      default=False"""
    params = PRIMER3_DEFAULTS.copy()
    params['SEQUENCE_TEMPLATE'] = template.seq
    params['PRIMER_PRODUCT_SIZE_RANGE'] = '-'.join(map(str, size_range))
    if options:
        params.update(options)
    p = subprocess.Popen([paths.PRIMER3_PATH],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output = p.communicate(format_primer3_params(params))[0]
    return dna.Pcr(best_pair(output, blast), template)

def primer3_clone(template):
    """PCR with fixed endpoints, PRIMER_TASK=pick_cloning_primers"""
    return primer3_pcr(template,
            size_range=[str(len(template.seq))] * 2,
            options={'PRIMER_TASK':'pick_cloning_primers'})
         
def format_primer3_params(params):
    """Format parameter dictionary for primer3 input"""
    # Must terminate input with empty assignment
    return '\n'.join(['='.join(map(str, (k,v)))
        for k, v in params.iteritems()] + ['='])

def best_pair(output, blast):
    """Return highest scoring primer pair from primer3 output"""
    pair = (dna.Primer(''), dna.Primer(''))
    first_pair = None
    # dicts decode primer3 parameters into dna.Primer attributes
    side = {'LEFT':0,'RIGHT':1}
    attrs = {'SEQUENCE':'base','GC_PERCENT':'gc','TM':'tm'}
    index = 0
    for line in output.split('\n'):
        m = re.search(PRIMER3_RE, line)
        if m:
            d = m.groupdict()
            if int(d['index']) > index and all([p.seq for p in pair]):
                if not blast or passes_blast_filter(pair):
                    # Pair is populated and passes filters
                    return pair
                else:
                    if index == 0:
                        # In case none pass filter,
                        # fall back on the first result
                        first_pair = pair
                    index += 1
                    pair = (dna.Primer(''), dna.Primer(''))
            if d['side'] in side and d['attr'] in attrs:
                pair[side[d['side']]].__setattr__(attrs[d['attr']],d['value'])
    if not first_pair:
        raise ValueError('Failed to find primer attributes in Primer3 output.')
    return first_pair

def passes_blastn_filter(pair):
    """Returns bool: primer pair predicted to mis-amplify in genome."""
    #TODO
    pass
    
def load_primer3_defaults():
    d = dict(
    PRIMER_PICK_LEFT_PRIMER=1,
    PRIMER_PICK_INTERNAL_OLIGO=0,
    PRIMER_PICK_RIGHT_PRIMER=1,
    PRIMER_PICK_ANYWAY=1,
    PRIMER_MIN_TM=50,
    PRIMER_MAX_TM=65
    )
    d['PRIMER_THERMODYNAMIC_PARAMETERS_PATH']=paths.PRIMER3_CONFIG_PATH
    return d

def faidx(region):
    """Wrapper for samtools faidx for FASTA random access.
    Will make .fai index only if it doesn't exist."""
    fasta = subprocess.check_output([paths.SAMTOOLS_PATH,
        'faidx', paths.GENOME_PATH, region])
    return ''.join(fasta.strip().split('\n')[1:]).lower()

# Run on import
make_genome_blastdb()
PRIMER3_DEFAULTS = load_primer3_defaults()
