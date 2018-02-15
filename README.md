locusmod on Flask
======
**locusmod** is a web app that automates design of plasmids that can **mod**ify a strain genotype at any **locus** of interest.

Supported Use Cases:
1. Knockout
2. Knockdown (DAmP)
3. Overexpression (Promoter swap)
4. N-terminal tagging
5. C-terminal tagging

External Data Requirements
------
**locusmod** needs the following external data to translate user input into plasmid designs:
```bash
static/
      genome/
            ASM2700v1.fa # GS115 genome
            ASM2700v1_genes.gff3 # GS115 gene annotation
      plasmid/
            # plasmid files and other lengthy inputs used for design
      testdata/
            # sequence files for unit testing
      tool/
            # see instructions in next section
```

Output Data (Downloads)
------
```bash
static/output/
              primers.csv
              idt_import.csv
              plasmids.csv
              plasmids.zip
              loci.zip
```

Python Dependencies
------
Python-2.7.10
Flask-0.12.2
biopython-1.64

Compiled Dependencies
------

0. Navigate to directory for binaries:
```bash
mkdir -p static/tool
cd static/tool
```

1. [Primer3](http://primer3.sourceforge.net/primer3_manual.htm#installLinux)
```bash
# download via browser
https://sourceforge.net/projects/primer3/files/primer3/2.3.7/primer3-2.3.7.tar.gz
tar -zxvf primer3-2.3.7.tar.gz 
cd primer3-2.3.7/src
make all # creates primer3_core, ntdpal, olgotm, and long_seq_tm_test
make test # plan on 10-30min to run all tests
```
2. [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
```bash
curl "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz" > "ncbi-blast-2.6.0+-x64-linux.tar.gz"
# you could also verify the md5 digital signature
tar -zxvf "ncbi-blast-2.6.0+-x64-linux.tar.gz"
```
3. [Samtools](http://www.htslib.org/download/)
```bash
sudo apt-get install samtools
```
