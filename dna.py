"""DNA objects
Attributes:
    alphabet: allowed characters for DNA
    dna_re: regex for allowed characters
    comp_table: string translation table for complement
Methods:
    revc: reverse complement for strings, used by objects
    validate_seq
    validate_fragments
    merge
Classes:
    Dna
        CutSite
        Primer
        Assembly
            Gibson
            Ligation
        Fragment    
            Digest
            Pcr
                SoePcr"""

import re
import string

alphabet='AaCcGgTt'
dna_re = re.compile('^['+alphabet+']*$')
comp_table = string.maketrans(alphabet,'TtGgCcAa')
revc = lambda s: s[::-1].translate(comp_table)

def validate_seq(seq):
    """Check Dna object sequence property data. Enforce lowercase."""
    if not isinstance(seq, str):
        raise TypeError('seq must be str.')
    if re.match(dna_re,seq) is None:
        raise ValueError('seq must only consist only of '+alphabet+'.')
    return seq.lower()

def validate_fragments(fragments):
    if any([not isinstance(f, Dna) for f in fragments]):
        raise ValueError('fragments must be iterable of Dna.')
    return fragments

def merge(a, b, min_overlap=0, max_overlap=80):
    """Merge A,B deleting from B the end of A that repeats in B.
    Note: O(n^2) time but n restricted by max_overlap.
    Raises ValueError if A and B do overlap by at least min_overlap
    or at most max_overlap."""
    len_match = 0
    i_max = min([max_overlap]+map(len,(a,b)))
    i = 0
    while i < i_max:
        i += 1
        # slice out of bounds OK. Case-aware for general strings.
        if a[-i:] == b[:i]:
            len_match = i
    if len_match < min_overlap:
        raise ValueError(
                'Require min_overlap={}, observed overlap={}'.format(
                    min_overlap,len_match))
    return a + b[len_match:]

class Dna(object):
    """Named DNA sequence with reverse complement"""
    def __init__(self, base, name='', circular=False):
        """Args: 
            base (str): required, must be minimal DNA alphabet sequence
            name (str): annotation, default None
            circular (bool): sequence string topology
        Properties:
            seq (str): gets sequence, cannot set;
                       in base class, seq is just base
            rep_seq (str): tandem repeat of sequence for
                           circular string search
            search_seq (str): allows literal string search
                              regardless of circular"""
        self._base = validate_seq(base)
        self._name = str(name)
        self.circular = circular
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        self._name = str(value)
    @property
    def base(self):
        return self._base
    @base.setter
    def base(self,value):
        self._base = validate_seq(value)
    @property
    def seq(self):
        return self._base
    @property
    def rep_seq(self):
        return self.seq*2
    @property
    def search_seq(self):
        if self.circular:
            return self.rep_seq
        return self.seq

class CutSite(Dna):
    """Restriction enzyme cut site"""
    def __init__(self, site, name=''):
        """Args:
            site (str): restriction site nucleotide sequence
            name (str): enzyme name
        Properties:
            palindromic (bool): True if palindrome by reverse complement"""
        super(CutSite, self).__init__(base=site, name=name, circular=False)
    @property
    def palindromic(self):
        return self.seq == revc(self.seq)

class Primer(Dna):
    """Allows separation of base and extension strings. Convention 5' to 3'"""
    def __init__(self, base, ext='', name=''):
        """Args:
            base (str): required, the template-binding portion of the primer
            ext (str): extension to the base of the primer, default ''
            name (str): annotation, default None
        Properties:
            ext (str): see Args
            seq (str): derived as ext + base; cannot set
        Attributes:
            gc (str): %GC from Primer3
            tm (str): Melting temperature from Primer3"""
        super(Primer, self).__init__(base=base, name=name)
        self._ext = validate_seq(ext)
        self.gc = self.tm = None
    @property
    def ext(self):
        return self._ext
    @ext.setter
    def ext(self,value):
        self._ext = validate_seq(value)
    @property
    def seq(self):
        return self.ext + self.base

class Fragment(Dna):
    """Base class for template extracted by pair of flanking sequences.
    Properties:
        f (Dna): forward (starting) member of pair
        r (Dna): reverse (ending) member of pair
        min_re (_sre.SRE_PATTERN): compiled regular expression
                                   for sequence extraction
        **seq must be implemented in subclasses
        name (str): fragment name, derived from other attributes."""
    def __init__(self, pair, template, name=''):
        """Args:
            pair (list of Dna): pair of flanking sequences
            template (Dna): sequence from which to extract fragment"""
        self.pair = pair
        if not isinstance(template, Dna):
            raise ValueError('Template must be Dna.')
        self.template = template
        super(Fragment, self).__init__(base='', circular=False, name=name)
    @property
    def instruction(self):
        return 'Fragment of {} with {} and {}'.format(
                self.template.name, self.f.name, self.r.name)
    @property
    def f(self):
        return self.pair[0]
    @property
    def r(self):
        return self.pair[1]
    @property
    def min_re(self):
        return re.compile(
                self.f.base+'(['+alphabet+']*?)'+revc(self.r.base),
                re.IGNORECASE)
    @property
    def seq(self):
        raise NotImplementedError(
                'seq must be implemented by subclasses of Fragment.')
    def extract(self, n_allowed, no_ext=False):
        """Returns template sequence extracted by flanking pair
        if it uniquely exists. Uses base to search, but returns
        full flankings seqs surrounding internal seq.
        Args:
            n_allowed (int): number of bindings of each flanking element
                             allowed. Defined in subclasses."""
        for p in self.pair:
            n_bindings = sum([self.template.search_seq.count(s)
                for s in (p.base, revc(p.base))])
            if n_bindings > n_allowed:
                raise ValueError(
                        'Pair flanks multiple sequences on template.')
        # Nongreedy regex to return shortest matches
        # (greedy will span both repeats). Search both strands.
        internals = self.min_re.findall(self.template.search_seq)
        internals += self.min_re.findall(revc(self.template.search_seq))
        if not internals:
            raise ValueError('Pair not found on template.')
        # If multiple products, return that on the top strand (first found)
        # to keep pair orientation
        if no_ext:
            return self.f.base + internals[0] + revc(self.r.base)
        else:
            return self.f.seq + internals[0] + revc(self.r.seq)

class Digest(Fragment):
    """Restriction digest product DNA sequence with enzyme pair and template.
    """
    def __init__(self, cutsite_pair, template, name=''):
        """Args:
            cutsite_pair (list of CutSite): beginning, ending cut sites
            template (Dna): sequence to be digested
            name (str): identifier of the Digest
        Properties:
            seq (str): the sequence cut out of the template
                       between start and end cut sites."""
        if len(cutsite_pair) != 2 or any([
            not isinstance(c, CutSite) for c in cutsite_pair]):
            raise ValueError('CutSite pair must be iterable of 2 CutSites.')
        super(Digest, self).__init__(
                pair=cutsite_pair, template=template, name=name)
    @property
    def instruction(self):
        return 'Digest {} with {} and {}'.format(
                self.template.name, self.f.name, self.r.name)
    @property
    def seq(self):
        # Cut site found >N times will produce unwanted fragments
        # For circular template, N=2. For linear, N=1.
        # For palindromic cutsite, N=2. Nonpalindromic N=1.
        # Multiplying the factors from circular and palindromic
        #     gives the total N allowed.
        # For a circular template searching palindromic cutsites,
        #     N=(1+1)*(1+1)=4.
        palindromic = any([c.palindromic for c in self.pair])
        return self.extract((self.template.circular+1)*(palindromic+1))

class Pcr(Fragment):
    """PCR product DNA sequence with primer pair and template"""
    def __init__(self, primer_pair, template, name=''):
        """Args:
            primer_pair (list of Primer): forward and reverse primers
            template (Dna)
            name (str): identifier of the Pcr
        Properties:
            f (Primer)
            r (Primer)
            seq (str): the sequence amplified from the template
                       with the primer_pair. """
        if len(primer_pair) != 2 or any([
            not isinstance(p, Primer) for p in primer_pair]):
            raise ValueError('Primer pair must be iterable of 2 primers.')
        super(Pcr, self).__init__(
                pair=primer_pair, template=template, name=name)
    @property
    def instruction(self):
        return 'PCR from {} with {} and {}'.format(
                self.template.name, self.f.name, self.r.name)
    @property
    def seq(self):
        # A primer binding the search sequence >N times
        # has >=50% mispriming, possibly multiple products.
        # For circular template, N=2. For linear, N=1.
        # Primer cannot bind in both directions.
        return self.extract(self.template.circular + 1)
    @property
    def seq_no_ext(self):
        # No extensions on primers, only base (use for sequence annotation)
        return self.extract(self.template.circular + 1, no_ext=True)

class SoePcr(Pcr):
    """SOE-PCR product DNA sequence with primer pair and fragments"""
    def __init__(self, primer_pair, fragments, name=''):
        """Args:
        primer_pair (list of Primer): used to amplify() sequence
        fragments (list of Dna): sequences to splice
        name (str): identifier for the resulting fragment"""
        self.fragments = validate_fragments(fragments)
        template = Assembly(fragments)
        super(SoePcr, self).__init__(primer_pair, template, name=name)
    @property
    def instruction(self):
        return 'SOE-PCR from {} with {} and {}'.format(
                self.template.instruction, self.f.name, self.r.name)

class Assembly(Dna):
    """Assembly of Dna fragments. Maintains references to names
    for instructions.
    Attributes:
        circularize (bool): whether to circularize the assembly
        min_overlap (int): bound passed to merge subroutine
        max_overlap (int): same as above"""
    circularize = False
    min_overlap = 0
    max_overlap = 100
    def __init__(self, fragments):
        """Args:
            fragments (list of Dna): fragments to assemble"""
        self.fragments = validate_fragments(fragments)
        super(Assembly, self).__init__(
                self.assemble(), circular=self.circularize)
        # Create a copy for modification 
    def assemble(self):
        """Assemble fragment ends with overlapping sequence. Returns seq."""
        fragments = list(self.fragments)
        if any([d.circular for d in fragments]):
            raise ValueError('Cannot assemble with circular Dna.')
        while len(fragments)>1:
            # pop the top 2 and merge in order
            fragments.append(
                    Dna(merge(fragments.pop(-2).seq, fragments.pop().seq,
                              self.min_overlap, self.max_overlap)))
        if self.circularize:
            final = fragments.pop().seq
            half = len(final)/2
            # cut in half and merge halves together at ends, 180 rotation
            merged = merge(final[half:], final[:half],
                    self.min_overlap, self.max_overlap)
            # Rotate by another 180 degrees.
            # Resulting appearance is deletion at start.
            half = (len(final)+1)/2
            fragments.append(Dna(merged[half:]+merged[:half]))
        return fragments.pop().seq
    @property
    def instruction(self):
        return ','.join([f.name for f in self.fragments])

class Ligation(Assembly):
    """__init__ from Assembly, just change class attributes and instruction"""
    circularize = True
    min_overlap = 6
    max_overlap = 8
    @property
    def instruction(self):
        return 'Ligate ' + super(Ligation, self).instruction

class Gibson(Assembly):
    """__init__ from Assembly, just change class attributes and instruction"""
    circularize = True
    min_overlap = 30
    @property
    def instruction(self):
        return 'Gibson assemble ' + super(Gibson, self).instruction
