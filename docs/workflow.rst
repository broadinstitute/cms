Sample workflow
=============================================================

**Data formatting**

CMS requires the user to provide population genetic (i.e., within-species) diversity data, including genotype phase and allele polarity. 

	-If your dataset is unphased, you can preprocess it using a program like `Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`_ or `PLINK <https://pngu.mgh.harvard.edu/~purcell/plink/>`_. 
	-The identity of the ancestral allele is typically determined by comparison to outgroups at orthologous sites. Inferred ancestral sequence is available for a number of species through `Ensembl <http://ensembl.org>`_ via their `ftp <ftp://ftp.ensembl.org/pub/release-84/fasta/ancestral_alleles/>`_.
	-In most cases, the user will want to provide a genetic recombination map. If this is unavailable, CMS will assume uniform recombination rates when calculating haplotype scores. Human recombination maps are available from the `HapMap Project <http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/>`_.
	-CMS works with `TPED <http://varianttools.sourceforge.net/Format/Tped>`_ datafiles, and includes support to convert from `VCF <http://samtools.github.io/hts-specs/VCFv4.3.pdf>`_ using the command line tool scans.py.


**Calculating selection statistics**
	-> scans.py (should I build off of this or make new scripts?)

**Demographic modeling**
	(-> make a wrapper method, model.py? use argparse?)

**Combining scores**

**Identifying regions**

**Localizing signals**