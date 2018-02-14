"""Tests Locus Mod objects and use cases"""

import batch
import dna
import locus
import paths
import re
import tool
import unittest

class DnaTestCase(unittest.TestCase):
    """Test dna.py"""
    def setUp(self):
        self.s = 'GAGAtact'
        self.d = dna.Dna(self.s, name='test')
        self.asci, self.sbfi = batch.V_CUTSITE_PAIR
        self.bsai = dna.CutSite('ggtctc', 'BsaI')
    def test_revc(self):
        self.assertEquals(dna.revc(self.s),'agtaTCTC')
    def test_Dna(self):
        self.assertRaises(TypeError, dna.Dna, [self.s])
        self.assertRaises(ValueError, dna.Dna, self.s + 'n')
        with self.assertRaises(ValueError):
            self.d.base = self.d.seq + 'n'
        with self.assertRaises(AttributeError):
            dna.Dna(self.s).seq = 'atg'
        self.assertNotIn('tga',self.d.search_seq.lower())
        self.d.circular = True
        self.assertIn('tga',self.d.search_seq.lower())
    def test_CutSite(self):
        self.assertTrue(self.asci.palindromic)
        self.assertFalse(self.bsai.palindromic)
    def test_Primer(self):
        self.assertEquals(dna.Primer('TAA',ext='atg').seq, 'atgTAA'.lower())
        with self.assertRaises(ValueError):
            dna.Primer('TAA').ext = 'NNN'
    def test_Digest(self):
        a = dna.CutSite('ata')
        c = dna.CutSite('cgc')
        d = dna.Dna('AataAAAgcgA')
        x = dna.Digest((a,c), d)
        self.assertEquals(x.seq, 'ataAAAgcg'.lower())
        # Nonpalindromic, searches opposite strand like primer pair
        x.pair = (c,a)
        self.assertEquals(x.seq, 'cgcTTTtat'.lower())
        # CutSite C tandem palindrome
        d.base = 'AataAAAgcgcA'
        with self.assertRaises(ValueError):
            dna.Digest((a,c), d).seq
        # Palindromic
        a = dna.CutSite('aatt')
        c = dna.CutSite('ccgg')
        d = dna.Dna('AaattAAAccggA')
        x = dna.Digest((a,c), d)
        self.assertEquals(x.seq, 'aattAAAccgg'.lower())
        d.circular = True
        # Now two products made, but returns top strand match only
        self.assertEquals(x.seq, 'aattAAAccgg'.lower())
        # Flipping strand prioritizes the shorter match around the circuit
        d.base = dna.revc(d.seq)
        self.assertEquals(x.seq, 'aattTTccgg'.lower())
    def test_Pcr(self):
        non = dna.Primer('tact')
        pal = dna.Primer('aatt')
        d = dna.Dna('AAAtactGGGaattAAA')
        pcr = dna.Pcr((non,pal), d)
        # palindromic primer has too many bindings
        with self.assertRaises(ValueError):
            pcr.seq
        new = dna.Primer('aattcc')
        pcr.pair = (non, new)
        self.assertEqual(pcr.seq, 'tactgggaatt')
        # Swapping pair returns reverse complement of product
        pcr.pair = pcr.pair[::-1]
        self.assertEqual(pcr.seq, dna.revc('tactgggaatt'))
        # Reversing each primer's orientation will raise exception
        # unless circular
        pcr.pair = [dna.Primer(dna.revc(p.seq)) for p in pcr.pair]
        with self.assertRaises(ValueError):
            pcr.seq
        pcr.template.circular = True
        self.assertEquals(pcr.seq, 'ggaattaaaaaatact')
    def test_SoePcr(self):
        # w-x-y will assemble, but will not SOE-PCR
        # due to interior primer binding
        w = dna.Dna('aaaccc')
        x = dna.Dna('cccttt')
        y = dna.Dna('tttggg')
        self.assertEquals(dna.Assembly((w,x,y)).seq, 'aaaccctttggg')
        f = dna.Primer('aaa')
        r = dna.Primer('ccc')
        soe = dna.SoePcr((f,r), [w,x,y])
        with self.assertRaises(ValueError):
            soe.seq
        # increasing specificity allows SOE-PCR
        f.base = 'aaac'
        r.base = 'ccca'
        self.assertEqual(soe.seq, 'aaaccctttggg')
    def test_Assembly(self):
        w = dna.Dna('aaaccc')
        x = dna.Dna('cccttt')
        y = dna.Dna('tttggg')
        z = dna.Dna('gggaaa')
        a = dna.Assembly((w,x,y,z))
        self.assertFalse(a.circular)
        a.circularize = True
        a.min_overlap = 0
        # string appearance in linear form should be
        # deletion of starting overlap
        self.assertEquals(a.assemble(), 'ccctttgggaaa')
    def test_Ligation(self):
        a8 = 'a'*8
        c8 = 'c'*8
        ac = dna.Dna(a8 + 'AC' + c8)
        ca = dna.Dna(c8 + 'CA' + a8)
        # Assembly will greedily merge homology,
        # while Ligation does between 6 and 8
        self.assertEquals(dna.Assembly([ac, ca]).seq,
                (a8 + 'A' + c8 + 'CA' + a8).lower())
        self.assertEquals(dna.Ligation([ac, ca]).seq,
                ('AC' + c8 + 'CA' + a8).lower())

class LocusTestCase(unittest.TestCase):
    """Test sequence retrieval in locus.py and tool.py:
    tool.blastn, tool.faidx"""
    def setUp(self):
        # short genes on each strand. Seqs from Ensembl in static/testdata/
        self.context = 100
        self.plus = {}
        self.plus['id'] = 'PAS_chr2-1_0751' 
        self.plus['coord'] = '2:1427323-1427517(+)'
        self.minus = {}
        self.minus['id']  = 'PAS_chr1-4_0504'
        self.minus['coord'] = '1:2322810-2322980(-)'
        self.dicts = (self.plus, self.minus)
        for d in self.dicts:
            d['locus'] = locus.Locus(dna.Dna('',name=d['id']),
                    context=self.context)
            with open(
                paths.TESTDATA_PATH + d['id'] + '.seq', 'r') as f:
                d['seq_file'] = f.read().strip()
            with open(
                paths.TESTDATA_PATH + d['id'] + '_context100.seq', 'r') as f:
                d['context_file'] = f.read().strip()
    def test_search_by_identifier(self):
        # Cross-reference locus retrieval to files from Ensembl webpage"""
        for d in self.dicts:
            self.assertEquals(d['locus'].name, d['coord'])
            self.assertEquals(d['locus'].seq.lower(), d['seq_file'].lower())
            full = d['locus'].up + d['locus'].seq + d['locus'].down
            self.assertEquals(full.lower(), d['context_file'].lower())
    def test_search_by_seq(self):
       #  Ensure that searching by sequence using blastn gives same results"""
        for d in self.dicts:
            seqloc = locus.Locus(dna.Dna(d['locus'].seq), context=self.context)
            self.assertEquals(seqloc.name, d['coord'])
            self.assertEquals(seqloc.seq.lower(), d['seq_file'].lower())
            full = d['locus'].up + d['locus'].seq + d['locus'].down
            self.assertEquals(full.lower(), d['context_file'].lower())
    def test_blastn(self):
        p = self.plus['locus']
        # If indel in the sequence, retrieves corrected sequence
        indel = p.seq[:90] + p.seq[100:]
        seqloc = locus.Locus(dna.Dna(indel))
        self.assertEquals(seqloc.name, self.plus['coord'])
        self.assertEquals(seqloc.seq.lower(), self.plus['seq_file'].lower())
        with self.assertRaises(ValueError):
            # Nonsense 
            locus.Locus(dna.Dna('tact'*1000))
        with self.assertRaises(ValueError):
            # Paste together in wrong order
            locus.Locus(dna.Dna(p.down + p.seq + p.up))

class Primer3TestCase(unittest.TestCase):
    """Test the primer3 options accessible from tool.py.
    Warning: primer3 may not be deterministic.
    Test logic rather than sequence equality."""
    pass

class SoePcrKoDesignTestCase(unittest.TestCase):
    """Existing design strategy."""
    def setUp(self):
        self.b = batch.Batch('SOE-PCR','Knockout','Lox-ZeoR','PAS_chr2-1_0666')
        self.d = self.b.designs[0]
    def test_assemble(self):
        calculated = sum(
                     map(len, (self.d.U.seq, self.d.M.seq, self.d.D.seq))
                     ) - len(self.d.U.r.ext) - len(self.d.D.f.ext)
        self.assertEquals(len(self.d.soepcr.seq), calculated)
    def test_write(self):
        #batch_script(self.b)
        pass

class GibsonKoDesignTestCase(unittest.TestCase):
    """New design strategy."""
    def setUp(self):
        self.b = batch.Batch('Gibson','Knockout','Lox-ZeoR','PAS_chr2-1_0666')
    def test_assemble(self):
        # TODO needs thought on automated checking. Visually validated in ApE.
        pass
    def test_write(self):
        #batch_script(self.b)
        pass

class BatchTestCase(unittest.TestCase):
    """Multiple designs in one batch"""
    def setUp(self):
        self.b = batch.Batch(
            'Gibson','Knockout','Lox-ZeoR','PAS_chr1-4_0370\nPAS_chr3_1028')
    def test_write(self):
        batch_script(self.b)

class PreciseKoTestCase(BatchTestCase):
    """FASTA input"""
    def setUp(self):
        with open(paths.TESTDATA_PATH + 'ku70_precise.fa') as f:
            self.b = batch.Batch('Gibson','Knockout','Lox-ZeoR',
                    target_text=f.read())
    def test_precise(self):
        d = self.b.designs[0]
        pat = re.compile(''.join(
                [d.U.seq_no_ext, '(.*)', d.D.seq_no_ext]))
        groups = re.findall(pat, d.locus.up + d.locus.seq + d.locus.down)
        self.assertEquals(groups[0], d.target.seq)

class KnockdownTestCase(BatchTestCase):
    """No target, but interrupt 3'UTR with marker"""
    def setUp(self):
        self.b = batch.Batch(
            'Gibson','Knockdown','Lox-ZeoR','PAS_chr1-4_0370')
    def test_knockdown(self):
        d = self.b.designs[0]
        pat = re.compile(''.join(
                [d.U.seq_no_ext, '(.*)', d.D.seq_no_ext]))
        groups = re.findall(pat, d.locus.up + d.locus.seq + d.locus.down)
        self.assertTrue(groups)
        self.assertTrue(0 <= len(groups[0]) <= d.down_out)

class OverexpressionTestCase(BatchTestCase):
    """Target is a promoter (with Kozak)"""
    def setUp(self):
        with open(paths.TESTDATA_PATH + 'pThi11.fa') as f:
            self.b = batch.Batch('Gibson','Overexpression','Lox-ZeoR',
                    target_text='PAS_chr2-1_0666', tag_text=f.read())
    def test_Overexpression(self):
        d = self.b.designs[0]
        # Check that homology arms delete less than upper bound bp
        ha_pat = re.compile(''.join(
                [d.U.seq_no_ext, '(.*)', d.D.seq_no_ext]))
        groups = re.findall(ha_pat, d.locus.up + d.locus.seq + d.locus.down)
        self.assertTrue(groups)
        self.assertTrue(0 <= len(groups[0]) <= d.up_out)
        # Check that promoter tags start
        orf_pat = re.compile(''.join([d.T.seq_no_ext + 'atg']))
        match = re.search(orf_pat, d.seq)
        self.assertTrue(match)

class NtagDesignTestCase(BatchTestCase):
    """Target is N-terminal tag including promoter-ATG"""
    def setUp(self):
        with open(paths.TESTDATA_PATH + 'Ndegron.fa') as f:
            self.b = batch.Batch('Gibson','Ntag','Lox-ZeoR',
                    target_text='PAS_chr2-1_0666', tag_text=f.read())
    def test_Ntag(self):
        d = self.b.designs[0]
        # Check that homology arms delete less than upper bound bp
        ha_pat = re.compile(''.join(
                [d.U.seq_no_ext, '(.*)', d.D.seq_no_ext]))
        groups = re.findall(ha_pat, d.locus.up + d.locus.seq + d.locus.down)
        self.assertTrue(groups)
        self.assertTrue(0 <= len(groups[0]) <= d.up_out)
        # Check that fusion is to first amino acid past start
        orf_pat = re.compile(''.join([d.T.seq_no_ext + d.target.seq[3:]]))
        match = re.search(orf_pat, d.seq)
        self.assertTrue(match)

class CtagDesignTestCase(BatchTestCase):
    """Target is C-terminal tag with STOP-terminator"""
    def setUp(self):
        with open(paths.TESTDATA_PATH + 'Dasher.fa') as f:
            self.b = batch.Batch('Gibson','Ctag','Lox-ZeoR',
                    target_text='PAS_chr2-1_0666', tag_text=f.read())
    def test_Ctag(self):
        d = self.b.designs[0]
        # Check that homology arms delete less than upper bound bp
        ha_pat = re.compile(''.join(
                [d.U.seq_no_ext, '(.*)', d.D.seq_no_ext]))
        groups = re.findall(ha_pat, d.locus.up + d.locus.seq + d.locus.down)
        self.assertTrue(groups)
        self.assertTrue(0 <= len(groups[0]) <= d.down_out)
        # Check that fusion is to last amino acid before start
        orf_pat = re.compile(''.join([d.target.seq[:-2] + d.T.seq_no_ext]))
        match = re.search(orf_pat, d.seq)
        self.assertTrue(match)

def batch_script(b):
    b.write_primers_csv()
    b.assign_numbers('RMo1')
    b.write_idt_csv()
    b.write_plasmids_zip()
    b.write_loci_zip()
    b.print_operations()

if __name__=='__main__':
    unittest.main(module='tests',verbosity=2, exit=False)
