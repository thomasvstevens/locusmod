"""Each form submission creates a Batch.
Each Batch contains a Design list.
Attributes and Methods of a Batch instance
suffice to build output pages and files.
Module globals are static file paths and dictionaries
for decoding selected form options.

Classes
    Batch
    AssemblyMethod
        SoePcrMethod
        GibsonMethod
    Design
        KoDesign
        NtagDesign
        CtagDesign
    (multiple inheritance of AssemblyMethod, Design)
    SoePcrKoDesign
    SoePcrCtagDesign
    GibsonKoDesign
    GibsonNtagDesign
    GibsonCtagDesign"""

from Bio import SeqIO
import copy
from cStringIO import StringIO
import dna
import itertools
import locus
import paths
import re
import subprocess
import tool

# Run on module import
with open(paths.UNS_PATH, 'r') as f:
    UNS = f.read().strip().lower().split('\n')
# Module globals: maps to decode form selectoptions
MAP_ASSEMBLY_METHOD = {'SOE-PCR':'SoePcr',
                       'Gibson':'Gibson'}
MAP_MOD_TYPE = {'Knockout':'Ko',
                'Knockdown':'Ctag',
                'Overexpression':'Ntag',
                'Ntag':'Ntag',
                'Ctag':'Ctag'}
MAP_PLASMID_PATH = dict(zip(
            ('Lox-ZeoR','Lox-NourR','Lox-G418R','Lox-HygroR','pUC-AmpR'),
            [paths.PLASMID_PREFIX + s for s in 
            ('RMp4694.gb','RMp4930.gb','RMp4957.gb','RMp4954.gb','RMp1087.ape')
            ]))
MAP_MOD_IN_FRAME = {'Knockout':False,
                    'Knockdown':False,
                    'Overexpression':False,
                    'Ntag':True,
                    'Ctag':True}
# Module globals: small dna objects and dna-related variables
M_BASE_PRIMER_PAIR = (dna.Primer('ctgattctgtggataaccgtagtc', name='RMo11228'),
                     dna.Primer('gaattggttaattggttgtaacacat', name='RMo11229'))
M_GIBSON_PRIMER_PAIR = copy.deepcopy(M_BASE_PRIMER_PAIR)
M_GIBSON_PRIMER_PAIR[0].ext = UNS[4]
M_GIBSON_PRIMER_PAIR[0].name = 'RMo11230'
M_GIBSON_PRIMER_PAIR[1].ext = dna.revc(UNS[5])
M_GIBSON_PRIMER_PAIR[1].name = 'RMo11231'
VECTOR = 'pUC-AmpR'
V_GIBSON_PRIMER_PAIR = (dna.Primer(ext=UNS[-1],
                                   base='ggcggtaatacggttatcca',
                                   name='RMo11232'),
                        dna.Primer(ext=dna.revc(UNS[0]),
                                   base='gcggaacccctatttgttta',
                                   name='RMo11233'))
V_CUTSITE_PAIR = (dna.CutSite('ggcgcgcc', name='AscI'),
                  dna.CutSite('cctgcagg', name='SbfI'))
INTEGRATION_CUTSITES = [dna.CutSite('gcggccgc', name='NotI'),
                        dna.CutSite('gtttaaac', name='PmeI'),
                        dna.CutSite('gcgatcgc', name='AsiSI'),
                        dna.CutSite('ttaattaa', name='PacI')]
SPACER_OUTER = 'gtcta' # padding to aid restriction digestion
SPACER_INNER = 'taaat' # AT region between GC rich restriction sites
GIBSON_OVERLAP = 40
# Module globals: Output formatting
HEADER_PRIMERS = 'Description,Sequence'
HEADER_IDT_PLATE = 'Well Position,Name,Sequence'
HEADER_IDT_TUBE = 'Name,Sequence,Scale,Purification'
HEADER_PLASMIDS='Description,Notes,Status,Resistance,Origin,APE File,Sequence'

class Batch(object):
    """Interface with application front-end.
    Encapsulates a batch of designs from a single form submission.
    Each is a physical plasmid for locus modification.
    Each is made by conjoining the following 6 sequence slots:
        U: Upstream aka 5'homology arm; used in all designs
        T (C): C-terminal custom sequence; used only in C-tagging
        M: Marker (dominant); used in all designs
        T (N): N-terminal custom sequence; used only in Overexpression
                and N-tagging designs
        D: Downstream aka 3'homology arm; used in all designs
        V: Vector backbone; used in all designs
    The accessible modifications are: 
        U_M_DV = Deletion (Knockout). Insertion = M
        U_MTDV = Overexpression (Promoter swap) and N-tagging. Insertion = MT
        UTM_DV = Knockdown (DAmP) and C-tagging. Insertion = TM.
                 C is empty for Knockdown"""
    def __init__(self, assembly_method, mod_type,
                 marker_name, target_text, tag_text=''):
        """Args:
            assembly_method (str): reaction scheme for assembling plasmid
                                    (Gibson or SOE-PCR)
            mod_type (str): modification type, constrained by assembly_method
                            (Knockout, etc.)
            marker (str): name that references a static marker
                          sequence file on server
            target_text (str): textarea input of gene IDs
                               or FASTA of nucleotide sequences
            tag_text (str): textarea input of tag nucleotide sequence
                            (for tagging mod_types)
        Attributes:
            designs (list of Design)"""
        self.assembly_method = MAP_ASSEMBLY_METHOD[assembly_method]
        self.mod_type = MAP_MOD_TYPE[mod_type]
        self.in_frame = MAP_MOD_IN_FRAME[mod_type]
        self.marker_name = marker_name
        self.targets = self.parse(target_text)
        self.tag = self.parse(tag_text)[0]
        self.common_ops()
        self.designs = []
        for target in self.targets:
            # Create hybrid instance dynamically
            args = (target, [self.M, self.T, self.V], self.in_frame)
            design = eval(self.assembly_method+self.mod_type+'Design')(*args)
            design.flanking_pcrs()
            design.colony_pcr()
            design.extend_fragments()
            design.add_cutsites()
            design.assemble()
            self.designs.append(design)
    def parse(self, text):
        """Parse textarea input to Dna objects.
        Must be either gene IDs, one per line, or FASTA records."""
        if not text:
            return [None]
        if text[0] == '>':
            # try to parse as FASTA records
            return [dna.Dna(str(r.seq), name=r.name)
                    for r in SeqIO.parse(StringIO(text), 'fasta')]
        else:
            return [dna.Dna('',name=t) for t in text.strip().split('\n')]
    def common_ops(self):
        """Performs Dna operations common to all Designs in a Batch"""
        basename = lambda f: f.name.split('/')[-1].split('.')[0]
        with open(MAP_PLASMID_PATH[self.marker_name]) as f:
            name = basename(f)
            self.M = dna.Dna(str(SeqIO.read(f,'genbank').seq),
                        name=name, circular=True)
        with open(MAP_PLASMID_PATH[VECTOR]) as f:
            name = basename(f)
            self.V = dna.Dna(str(SeqIO.read(f,'genbank').seq),
                        name=name, circular=True)
        if self.assembly_method=='SoePcr':
            self.M = dna.Pcr(M_BASE_PRIMER_PAIR, self.M, name='pcr_M')
            # Global is in insert order. Reverse for vector order.
            self.V = dna.Digest(V_CUTSITE_PAIR[::-1], self.V, name='dig_V')
        elif self.assembly_method=='Gibson':
            self.M = dna.Pcr(M_GIBSON_PRIMER_PAIR, self.M, name='pcr_M')
            self.V = dna.Pcr(V_GIBSON_PRIMER_PAIR, self.V, name='pcr_V')
        else:
            raise AttributeError('Unrecognized assembly method {}'.format(
                                    self.assembly_method))
        if self.tag:
            self.T = tool.primer3_clone(self.tag)
            self.T.label = self.tag.name
            self.T.name = 'pcr_T'
        else:
            self.T = None
        return 
    def write_primers_csv(self):
        """Write csv of all unique primers in design for LIMS upload."""
        self.primers = []
        num = 1
        with open(paths.PRIMERS_CSV, 'w') as f:
            f.write(HEADER_PRIMERS + '\n')
            for design in self.designs:
                # Obtain pcrs unique to this design by filtering on type
                # and excluding common fragments.
                # Need to compare type, as isinstance will check superclass. 
                pcrs = (frag for frag in design.operations
                        if type(frag) is dna.Pcr
                        and frag not in [self.M, self.V])
                for pcr in pcrs:
                    pcr.name = 'pcr_{}'.format(num)
                    num += 1
                    for primer, direction in zip(pcr.pair, ('F','R')):
                        primer.description = '{} {} {}'.format(
                                design.target.name, pcr.name, direction)
                        f.write(','.join([primer.description, primer.seq]))
                        f.write('\n')
                        self.primers.append(primer)
        if len(self.primers) > 96:
            raise AttributeError(
                    'Only 96 primers allowed per batch (per plate)')
        return
    def assign_numbers(self, start_id='A1'):
        """User input defines the starting numerical id to apply to primers.
        For now, maximum number of primers per batch is 96.
        Args:
            start_id (str): starting primer identifier received by mass import.
                            Default: 'A1'"""
        if not self.primers:
            raise AttributeError(
                'Need to write_primers_csv() and import to receive start_id.')
        prefix, number = split_start_id(start_id)
        rows = 'ABCDEFGH'
        columns = range(1, 13)
        if prefix in rows and number in columns:
            self.is_plate = True
            id_gen = itertools.chain(
                    itertools.product(prefix, xrange(number, 13)),
                    itertools.product(rows.split(prefix)[-1], columns)
                    )
        else:
            self.is_plate = False
            id_gen = ((prefix, int(number) + i) for i in xrange(96))
        for primer in self.primers:
            primer.name = '{}{}'.format(*id_gen.next())
        return
    def write_idt_csv(self):
        """Write csv of all unique primers in design for IDT order."""
        if self.is_plate:
            header = HEADER_IDT_PLATE
        else:
            header = HEADER_IDT_TUBE
        with open(paths.IDT_CSV, 'w') as f:
            if self.is_plate:
                f.write(HEADER_IDT_PLATE + '\n')
                for p in self.primers:
                    f.write(','.join([p.name, p.description, p.seq]))
                    f.write('\n')
            else:
                f.write(HEADER_IDT_TUBE + '\n')
                for p in self.primers:
                    # Scale up long primers to avoid IDT errors.
                    # Leave purification blank.
                    scale = '25nm' if len(p.seq) <= 60 else '100nm'
                    f.write(','.join([p.name, p.seq, scale, '']))
                    f.write('\n')
        return 
    def write_plasmids_zip(self, start_id='RMp1'):
        """Write csv of plasmids for LIMS upload,
        archived together with .zip of plasmids"""
        archive_contents = [paths.PLASMIDS_CSV]
        prefix, number = split_start_id(start_id)
        id_gen = ((prefix, int(number) + i) for i in xrange(len(self.designs)))
        with open(paths.PLASMIDS_CSV, 'w') as f:
            f.write(HEADER_PLASMIDS + '\n')
            for design in self.designs:
                ape_fname = '{}{}.gb'.format(*id_gen.next())
                design.assembly.name = ape_fname
                type_str = re.findall('\.(\w+).*', str(design))[0]
                desc = '{}. YeastR:{}. Locus:{}.'.format(
                        type_str, self.marker_name, design.target.name)
                f.write(','.join([desc, '', 'Under Construction', 'Amp',
                                 'pUC', ape_fname, '']))
                f.write('\n')
                design.annotate_plasmid(ape_fname).write(
                        paths.OUTPUT_PREFIX + ape_fname)
                archive_contents.append(paths.OUTPUT_PREFIX + ape_fname)
        subprocess.check_call(['zip', paths.PLASMIDS_ZIP] + archive_contents)
        return
    def write_loci_zip(self):
        """Write an archive of all the locus maps for a design,
        annotated with primers."""
        archive_contents = []
        for design in self.designs:
            ape_fname = design.target.name + '_locus.gb'
            # SeqIO will not write "long" identifiers. Call it locus_map.
            design.annotate_locus('locus_map').write(
                    paths.OUTPUT_PREFIX + ape_fname)
            archive_contents.append(paths.OUTPUT_PREFIX + ape_fname)
        subprocess.check_call(['zip', paths.LOCI_ZIP] + archive_contents)
        return
    def print_operations(self):
        """Prints (name, instruction) pairs for common batch operations,
        then each design in the batch."""
        print('BATCH OPERATIONS')
        for op in [self.M, self.V]:
            if op:
                print('{}\t{}'.format(op.name, op.instruction))
        for design in self.designs:
            print('DESIGN {}'.format(design.assembly.name))
            for op in design.operations:
                if op and op not in [self.M, self.V]:
                    print('{}\t{}'.format(op.name, op.instruction))

class AssemblyMethod(object):
    """Base class for specific assembly methods.
    Use as an interface for hybrid classes.
    Properties:
        seq (str): requires completion of assemble(). alias for assembly.seq
        operations (list of Dna): exhaustive list for instructions,
                                  not for assembly"""
    def __init__(self, fragments=None, **kw):
        """Args:
            fragments (list of Dna): supplied in order of intended assembly"""
        # super() is required to call next __init__() method in the
        # multiple resolution order. object.__init__() is no-op
        super(AssemblyMethod, self).__init__()
        # Defaults unset to allow sibling class to supply fragments attribute
        if fragments is not None:
            self.fragments = fragments
    @property
    def seq(self):
        return self.assembly.seq
    @property
    def operations(self):
        raise NotImplementedError('Subclass must override operations.')
    def extend_fragments(self):
        raise NotImplementedError('Subclass must override extend_fragments.')
    def add_cutsites(self):
        """Just before assembly, check for integration 8-cutters and
        append sites flanking cassette.
        NEVER add to vector, as this results in a high frequency of
        Gibson recircularization (8p homology)."""
        # Integration only here. For cloning, SoePcr will extend.
        for site in INTEGRATION_CUTSITES:
            if not any([site.seq in f.search_seq for f in self.fragments]):
                # Append integration cutsite to interior of extension of
                # homology PCR primers
                self.U.f.ext = self.U.f.ext + site.seq
                # Map to top strand, concatenate, then reverse complement again
                self.D.r.ext = dna.revc(''.join(
                    map(dna.revc, [site.seq, self.D.r.ext])
                    ))
                return
        raise Exception(
            'All supplied integration cutsites are '
            +'found in the assembly fragments.')
    def assemble(self):
        """Assembles plasmid from fragments.
        Default is naive linear string assembly."""
        fragments = dna.validate_fragments(self.fragments)
        self.assembly = dna.Assembly(self.fragments)
        return

class SoePcrMethod(AssemblyMethod):
    """Joins exactly 4 fragments: 3 by SOE-PCR, 1 by 8-cutter ligation.
    To maintain common marker primers,
    overlaps are simply the marker base primers.
    Attributes:
        soepcr (SoePcr)
        insert (Digest)
        vector (Digest)"""
    def __init__(self, fragments=None, **kw):
        """Args:
            fragments (list of Dna): supplied in U,M,D,V order"""
        super(SoePcrMethod, self).__init__(fragments)
    def validate_fragments(self):
        if len(self.fragments) != 4:
            raise ValueError('SoePcrMethod assembly is only compatible '
            +'with exactly 4 fragments.')
        return dna.validate_fragments(self.fragments)
    def extend_fragments(self):
        """Extend U.r with M.f, D.f with M.r identically,
        accounting for strand."""
        U,M,D,_ = self.fragments
        U.r.ext += dna.revc(M.f.seq)
        D.f.ext += dna.revc(M.r.seq)
        return
    def add_cutsites(self):
        super(SoePcrMethod, self).add_cutsites()
        #Also check for cloning 8-cutters. Extend with sites including spacer.
        U,_,D,_ = self.fragments
        U.f.ext = ''.join([SPACER_OUTER, V_CUTSITE_PAIR[0].seq,
                           SPACER_INNER, U.f.ext])
        # Map to top strand, concatenate, then reverse complement again
        D.r.ext = dna.revc(''.join(map(dna.revc, [D.r.ext, SPACER_INNER,
                                    V_CUTSITE_PAIR[1].seq, SPACER_OUTER])))
        return
    def assemble(self):
        """Assembles in 2 phases: SOE-PCR then Digestion-Ligation"""
        U,M,D,V = self.fragments
        self.soepcr = dna.SoePcr((U.f,D.r), [U,M,D], name='pcr_S')
        self.insert_digest = dna.Digest(
                V_CUTSITE_PAIR, self.soepcr, name='dig_S')
        # V is initialized as a Digest in Batch.__init__
        self.assembly = dna.Ligation([self.insert_digest, V])
        return
    @property
    def operations(self):
        return self.fragments + [
                self.soepcr, self.insert_digest, self.assembly, self.colony
                ]

class GibsonMethod(AssemblyMethod):
    """Gibson assembly joins several fragments with homologous overlaps,
    leveraging the Unique Nucleotide Sequences (UNS) from Torella et al. 2014.
    Attributes:
        uns_list (list of str): read from UNS constants file"""
    def __init__(self, fragments=None, **kw):
        super(GibsonMethod, self).__init__(fragments)
        # each design contains its own copy of the module UNS list
        self.uns_list = UNS[:]
    def validate_fragments(self):
        if len(self.fragments) > len(self.uns_list):
            raise ValueError('Cannot assemble with more fragments '
            +'than available connecting sequences (UNS).') 
        return dna.validate_fragments(self.fragments)
    def extend_fragments(self):
        """While any fragment is missing extensions in its primers,
        fill from previous and next.
        If previous extended but not next, use the next unused UNS.
        Encountering UNS removes from list.
        Fragment identities not required.
        If existing overlap has UNS, must be on top strand."""
        # Abuse list as circular linked list with modulo-n indexing
        n = len(self.fragments)
        if not any([p.ext for f in self.fragments for p in f.pair]):
            # If none of the fragments have extensions,
            # this is a generic design. Cycle through UNSs, O(n) time.
            for i in range(n):
                uns = self.uns_list.pop(0)
                self.fragments[i-1].r.ext += dna.revc(uns)
                self.fragments[i].f.ext += uns
            return
        # Maximally examines each fragment once, O(n) time.
        # There must be at least one fragment with extensions:
        # start at that fragment.
        i = (i for i,f in enumerate(self.fragments)
               for p in f.pair if p.ext).next()
        prv = self.fragments[(i-1)%n]
        cur = self.fragments[i]
        while not all([p.ext for f in self.fragments for p in f.pair]):
            i = self.fragments.index(cur)
            nxt = self.fragments[(i+1)%n]
            if cur.f.ext:
                # Propagate to previous node
                if prv.r.ext:
                    # Enforce assembly integrity:
                    if cur.f.ext != dna.revc(prv.r.ext):
                        try:
                            dna.merge(prv.seq, cur.seq, dna.Gibson.min_overlap)
                        except ValueError:
                            print('Input fragment overlaps do not assemble.')
                            raise
                else:
                    prv.r.ext = dna.revc(cur.f.ext)
            else:
                if not prv.r.ext:
                    # Assign the next unused UNS
                    prv.r.ext = dna.revc(staged_uns)
                # Propagate previous extension to current.
                cur.f.ext = dna.revc(prv.r.ext)
            # Remove used UNS and point staged_uns to next unused UNS
            # (circular linked list)
            if cur.f.ext and cur.f.ext in self.uns_list:
                u = self.uns_list.index(cur.f.ext)  
                self.uns_list.remove(cur.f.ext)
                staged_uns = self.uns_list[u%len(self.uns_list)]
            # Advance
            prv = cur
            cur = nxt
    def assemble(self):
        self.assembly = dna.Gibson(self.fragments)
        return
    @property
    def operations(self):
        return self.fragments + [self.assembly, self.colony]

class Design(object):
    """Base Class for Designs corresponding to distinct modification types.
    Both upstream (5') and downstream (3') ends have in/outside search radii.
    Attributes:
        locus (Locus): locus of modification identified by target
    Properties:
        fragments (list of Dna): only dna fragments for assembly
    Example 5' end:
          search              search  
     ___|____ ___|__________|out_ _in|_______ ...
    | search |      opt_arm      | locus.seq  ...
    |          locus.up          |"""
    opt_arm = 1000 #optimal homology arm length
    up_out = 100 #radius outside feature 5' end to pick primers
    up_in = 100 #radius inside feature 5' end to pick primers
    down_in = 150 #radius inside feature 3' end to pick primers
    down_out = 50 #radius outside feature 3' end to pick primers
    # by default, extra context is one additional primer search radius
    search = max(map(sum, ((up_out, up_in), (down_out, down_in)))) 
    def __init__(self, target, common_frags, in_frame=False):
        """Args:
            target (Dna): 
            common_frags (list of Dna): Dna fragments shared by Batch.
            in_frame (bool): True if tag should create a protein fusion
                            (maintain START or STOP codon)"""
        # super() required to call next __init__() method in the
        # multiple resolution order. (object.__init__() is no-op)
        super(Design, self).__init__()
        self.in_frame = in_frame
        self.target = target
        # Single retrieval of outermost sequences
        self.locus = locus.Locus(target, context=self.opt_arm + 2*self.search)
        self.M, self.T, self.V = common_frags
        # T overhangs are unique to the design, though Pcr common to batch.
        self.T = copy.deepcopy(self.T)
        self.U = self.D = self.colony = None
    @property
    def fragments(self):
        raise NotImplementedError('Subclass must determine fragment order.')
    def flanking_pcrs(self):
        """Precise flanking primer design, triggered by FASTA target input.
        By default, primer3 is 0-indexed.
        Must override for interior primer search radius."""
        size_range = (self.opt_arm, self.opt_arm + self.search)
        options = {}
        options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str, 
                    [len(self.locus.up) - 1 - self.opt_arm - (self.search - 1),
                    self.search,
                    len(self.locus.up) - 1 - (self.search - 1),
                    self.search]))
        # This should override the search radius above
        options['SEQUENCE_FORCE_RIGHT_START'] = len(self.locus.up) - 1
        template = dna.Dna(self.locus.up, name='P.p. gDNA')
        self.U = tool.primer3_pcr(template, size_range, options=options)
        options = {}
        options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str,
                    [0, self.search, self.opt_arm, self.search]))
        # This should override the search radius above
        options['SEQUENCE_FORCE_LEFT_START'] = 0
        template = dna.Dna(self.locus.down, name='P.p. gDNA')
        self.D = tool.primer3_pcr(template, size_range, options=options)
        return
    def colony_pcr(self):
        if not all([self.U, self.D]):
            raise AttributeError(
                  'Must first generate U and D fragments with flanking_pcrs().'
                  )
        full = self.locus.up + self.locus.seq + self.locus.down
        spacer = 30
        start = full.index(self.U.f.base) - spacer
        end = full.index(dna.revc(self.D.r.base)) + len(self.D.r.base) + spacer
        options = {'SEQUENCE_EXCLUDED_REGION':','.join(
                    map(str, (start, end - start))
                    )}
        # Also set this region as a target to prevent N/Ctag from sidestepping
        options['SEQUENCE_TARGET'] = options['SEQUENCE_EXCLUDED_REGION']
        size_range = (end - start, end - start + 2*self.search)
        self.colony = tool.primer3_pcr(dna.Dna(full, name='STRAIN'),
                size_range=size_range, options=options)
    def annotate_plasmid(self, id_str=''):
        """Genbank:ApE feature labeling from feature library as FASTA.
        Returns ApeFile, not saved as attribute of Design.
        raises AttributeError unless seq is inherited and assembly completed.
        """
        ape = tool.ApeFile(self.seq, id_str)
        ape.add_from_fasta()
        # annotate design variable portions: homology arms and tag (U,D,T).
        # exclude ext from Pcr.
        for frag, label in zip((self.U, self.D),('5','3')):
            ape.add_feature(frag.seq_no_ext, "{} {}'HA".format(
                self.target.name, label))
        if self.T:
            ape.add_feature(self.T.seq_no_ext, self.T.label)
        return ape
    def annotate_locus(self, id_str=''):
        """Annotate unmodified locus with primer names.
        Returns ApeFile, not saved as attribute of Design."""
        full = self.locus.up + self.locus.seq + self.locus.down
        ape = tool.ApeFile(full, id_str, max_colors=7)
        # annotate target
        ape.add_feature(self.locus.seq, '{} {}'.format(
            self.target.name, self.locus.name))
        # annotate U, D primers and colony pcr primers
        for frag in (self.U, self.D, self.colony):
            for primer in frag.pair:
                ape.add_feature(primer.base, primer.name)
        return ape

class KoDesign(Design):
    """assembly_method can be Gibson or SOEPCR"""
    def __init__(self, target, common_frags, in_frame=False):
        super(KoDesign, self).__init__(target, common_frags, in_frame)
        # retrieve extra context to allow for exterior colony PCR
    @property
    def fragments(self):
        return [self.U, self.M, self.D, self.V]
    def flanking_pcrs(self):
        if self.locus.precise:
            super(KoDesign, self).flanking_pcrs()
            return
        else:
            size_range = (self.opt_arm, self.opt_arm + self.search)
            options = {}
            options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str, 
                        [len(self.locus.up) - self.opt_arm - self.up_out,
                        self.up_out + self.up_in,
                        len(self.locus.up) - self.up_out,
                        self.up_out + self.up_in]))
            template = dna.Dna(
                    self.locus.up + self.locus.seq[:self.up_in],
                    name='P.p. gDNA')
            self.U = tool.primer3_pcr(template, size_range, options=options)
            options = {}
            options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str,
                        [0,
                        self.down_in + self.down_out,
                        self.opt_arm,
                        self.down_in + self.down_out]))
            template = dna.Dna(
                    self.locus.seq[-self.down_in:] + self.locus.down,
                    name='P.p. gDNA')
            self.D = tool.primer3_pcr(template, size_range, options=options)
            return

class NtagDesign(Design):
    """N-terminal tagging and overexpression.
    assembly_method must be Gibson."""
    # Search leeway only in upstream, outside end
    up_in = 0
    down_in = 0
    down_out = 0
    search = 200
    def __init__(self, target, common_frags, in_frame=False):
        super(NtagDesign, self).__init__(target, common_frags, in_frame)
    @property
    def fragments(self):
        return [self.U, self.M, self.T, self.D, self.V]
    def flanking_pcrs(self):
        if self.locus.precise:
            super(NtagDesign, self).flanking_pcrs()
            return
        size_range = (self.opt_arm, self.opt_arm + self.search)
        options = {}
        options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str, 
                    [len(self.locus.up) - 1 - self.opt_arm - (self.search - 1),
                    self.search,
                    len(self.locus.up) - 1 - (self.up_out - 1),
                    self.up_out]))
        template = dna.Dna(self.locus.up, name='P.p. gDNA')
        self.U = tool.primer3_pcr(template, size_range, options=options)
        options = {}
        options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str,
                    [0, self.search, self.opt_arm, self.search]))
        options['SEQUENCE_FORCE_LEFT_START'] = 0
        if self.in_frame:
            # Fusion to first codon after ATG
            options['SEQUENCE_FORCE_LEFT_START'] = 3
        # Append downstream in case the locus sequence is shorter than opt_arm
        template = dna.Dna(self.locus.seq + self.locus.down, name='P.p. gDNA')
        self.D = tool.primer3_pcr(template, size_range, options=options)
        return
    def extend_fragments(self):
        """First, extends primers to splice in_frame sequences between
        tag and D homology arm. Then calls superclass method to complete."""
        # Use seq_no_ext to avoid duplication of overlap when ext modified
        # revc to top strand, concatenate, then revc back to bottom
        if self.T is not None:
            self.T.r.ext = dna.revc(''.join(
                [dna.revc(self.T.r.ext),
                 self.D.seq_no_ext[:GIBSON_OVERLAP / 2]]))
            self.D.f.ext = ''.join(
                [self.T.seq_no_ext[- GIBSON_OVERLAP / 2:], self.D.f.ext])
        super(NtagDesign, self).extend_fragments()
        return

class CtagDesign(Design):
    """C-terminal tagging and knockdown by DAmP.
    assembly_method must be Gibson"""
    # Search leeway only in downstream, outside end
    up_out = 0
    up_in = 0
    down_in = 0
    # Increase the down_out search radius since terminators are low Tm
    down_out = 100
    search = 200
    def __init__(self, target, common_frags, in_frame=False):
        super(CtagDesign, self).__init__(target, common_frags, in_frame)
    @property
    def fragments(self):
        # T may be None for Knockdown (no tag Pcr required)
        return [frag for frag in (self.U, self.T, self.M, self.D, self.V)
                    if frag is not None]
    def flanking_pcrs(self):
        if self.locus.precise:
            super(CtagDesign, self).flanking_pcrs()
            return
        size_range = (self.opt_arm, self.opt_arm + self.search)
        options = {}
        # Include upstream in template in case locus seq shorter than opt_arm
        template = dna.Dna(self.locus.up + self.locus.seq, name='P.p. gDNA')
        options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str, 
                    [len(template.seq) - 1 - self.opt_arm - (self.search - 1),
                    self.search,
                    len(template.seq) - 1 - (self.search - 1),
                    self.search]))
        options['SEQUENCE_FORCE_RIGHT_START'] = len(template.seq) - 1
        if self.in_frame:
            # Fusion to last codon before STOP
            options['SEQUENCE_FORCE_RIGHT_START'] = len(template.seq) - 1 - 3
        self.U = tool.primer3_pcr(template, size_range, options=options)
        options = {}
        options['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ','.join(map(str,
                    [0,
                    self.down_out,
                    self.opt_arm,
                    self.search]))
        # Append downstream in case the locus sequence is shorter than opt_arm
        template = dna.Dna(self.locus.down, name='P.p. gDNA')
        self.D = tool.primer3_pcr(template, size_range, options=options)
        return
    def extend_fragments(self):
        """First, extends primers to splice in_frame sequences between
        tag and U homology arm. Then calls superclass method to complete."""
        # Use seq_no_ext to avoid duplication of overlap when ext modified
        # revc to top strand, concatenate, then revc back to bottom
        if self.T is not None:
            self.U.r.ext = dna.revc(''.join(
                [dna.revc(self.U.r.ext), self.T.seq_no_ext[:GIBSON_OVERLAP]]))
            self.T.f.ext = self.U.seq_no_ext[-GIBSON_OVERLAP:] + self.T.f.ext
        super(CtagDesign, self).extend_fragments()
        return

# Multiple inheritance.
# Safer to depth-first resolve Design before AssemblyMethod.
SoePcrKoDesign = type('SoePcrKoDesign', (KoDesign, SoePcrMethod), {})
SoePcrCtagDesign = type('SoePcrCtagDesign', (CtagDesign, SoePcrMethod), {})
GibsonKoDesign = type('GibsonKoDesign', (KoDesign, GibsonMethod), {})
GibsonNtagDesign = type('GibsonNtagDesign', (NtagDesign, GibsonMethod), {})
GibsonCtagDesign = type('GibsonCtagDesign', (CtagDesign, GibsonMethod), {})

# General purpose methods
def split_start_id(start_id):
    groups_list = re.findall('(.*?)(\d+)', start_id)
    if groups_list:
        prefix, number = groups_list[0]
        number = int(number)
    else:
        prefix = start_id
        number = 1
    return prefix, number
