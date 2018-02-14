"""Paths and objects for interacting with genome and annotation
Classes:
    Locus
Methods:
    load_locus_dict
    build_locus_dict
Attributes:
    COORD_RE (re.SRE): regular expression for parsing GFF to coordinates
Run on import to supply LOCUS_DICT for lookups"""

import copy
import dna
import os
import paths
import pickle
import re
import tool

COORD_RE = re.compile('\t'.join([
                        '(?P<chromosome>\d+)','\w+','gene',
                        '(?P<start>\d+)',
                        '(?P<end>\d+)','\S',
                        '(?P<strand>[+-])','\S',
                        'ID=gene:(?P<identifier>[^;]+);']))

class Locus(dna.Dna):
    """Dna with chromosomal coordinates
    Properties:
        seq (str): sequence of target
    Attributes:
        precise (bool): if identifier not supplied, disable search radius
        chromosome (str): chromosome identifier in annotation
        start (int): coordinate start, 1-indexed inclusive for .gff, .fa
        end (ind): coordinate end, 1-indexed inclusive for .gff, .fa
        strand (str): + or - strand of feature
        up (str): upstream sequence (5' context)
        down (str): downstream sequence (3' context)
        """
    def __init__(self, target, context=0):
        """Args:
        target (dna.Dna): Dna defining either locus identifier or seq
        context (int): length in bp of flanking regions to extract"""
        self.context = context
        self.identifier = target.name
        if self.identifier in LOCUS_DICT:
            self.precise = False
            self.search_by_identifier()
        else:
            self.base = target.seq
            self.precise = True
            self.search_by_seq()
        self.name = '{}:{}-{}({})'.format(
                self.chromosome, self.start, self.end, self.strand)
    def search_by_identifier(self):
        """Get coordinates from gene identifier"""
        self.chromosome, self.start, self.end, self.strand = LOCUS_DICT[
                self.identifier]
        self.search_common()
        return
    def search_by_seq(self):
        """Run blastn on target sequence.
        Return Locus corresponding to top hit."""
        # Alignment may not be exact, so overwrite
        # the target sequence with the alignment in the reference
        locus_tuple, self.base = tool.blastn(self.base)
        self.chromosome, self.start, self.end, self.strand = locus_tuple  
        self.search_common()
        return
    def search_common(self):
        """Common cleanup operations for all searches"""
        self.start, self.end = map(int, (self.start, self.end))
        self.seq_from_coords()
        return
    def validate_strand(self):
        if self.strand not in ['+','-']:
            raise AttributeError(
                    'Locus coordinates strand must be one of + or -')
        return
    def seq_from_coords(self):
        """Update base, upstream, downstream sequences given locus coordinates.
        Use samtools faidx FASTA random access to save time and space."""
        # samtools faidx uses 1-indexed inclusive coordinates
        region = ''.join(map(str, [self.chromosome, ':',
                                   self.start - self.context, 
                                   '-',
                                   self.end + self.context]))
        region_seq = tool.faidx(region)
        base_len = self.end - self.start + 1
        self.up = region_seq[:self.context]
        self.base = region_seq[self.context:self.context + base_len]
        self.down = region_seq[self.context + base_len:]
        if self.strand == '-':
            self.up, self.base, self.down = map(dna.revc, (
                self.down, self.base, self.up))
        return

def load_locus_dict():
    """Returns serialized dict object. If not found, build from gff file."""
    if os.path.exists(paths.LOCUS_DICT_PKL_PATH):
        with open(paths.LOCUS_DICT_PKL_PATH, 'rb') as f:
            return pickle.load(f)
    else:
        return build_locus_dict()

def build_locus_dict():
    """Returns dict of loci, also writes object to .pkl
    format = {identifier: (chromosome, start, end, strand)}"""
    locus_dict = {}
    with open(paths.GFF_PATH, 'r') as f:
        for line in f:
            match = re.search(COORD_RE, line)
            if match:
                groups = list(match.groups())
                locus_dict[groups[-1]] = groups[:-1]
    # Need to open in binary mode for most efficient pickling protocols
    with open(paths.LOCUS_DICT_PKL_PATH, 'wb') as f:
        # -1 uses the latest available protcol
        pickle.dump(locus_dict, f, protocol=-1)
    return locus_dict

# Run on import
LOCUS_DICT = load_locus_dict()
