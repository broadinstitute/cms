Sample workflow
=============================================================

CMS provides a computational framework to explore signals of natural selection in diploid organisms. This document describes how to do so at an abstract level. {Link to a sample results paper?} For a more technical description, see {command line tools}.

Preliminary considerations
-------------------

- Sample size
- Defining populations (thin to same sample size )
- Coverage and SNP density (multiple testing)
- Population structure/relatedness
- Autosome vs etc. (gw vs per-chrom)


Data formatting
-------------------

CMS requires the user to provide population genetic (i.e., within-species) diversity data, including genotype phase and allele polarity. 

- If your dataset is **unphased**, you can preprocess it using a program like `Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`_ or `PLINK <https://pngu.mgh.harvard.edu/~purcell/plink/>`_. 
- The identity of the **ancestral allele** at each site is typically determined by comparison to outgroups at orthologous sites. Inferred ancestral sequence is available for a number of species through e.g. `Ensembl <http://ensembl.org>`_ via their `ftp <ftp://ftp.ensembl.org/pub/release-84/fasta/ancestral_alleles/>`_. You can use `VCFtools <https://github.com/vcftools/vcftools.github.io>`_ to populate the "AA=" section of your VCF's INFO field.
- In most cases, the user will want to provide a **genetic recombination map**. If this is unavailable, CMS will assume uniform recombination rates when calculating haplotype scores. Human recombination maps are available from the `HapMap Project <http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/>`_.
- CMS works with `TPED <http://varianttools.sourceforge.net/Format/Tped>`_ datafiles, and includes support to convert from `VCF <http://samtools.github.io/hts-specs/VCFv4.3.pdf>`_ using the command line tool scans.py. {LINK}


Calculating selection statistics
-------------------
	-> scans.py (should I build off of this or make new scripts?)


Demographic modeling
-------------------
	model.py? use argparse?
		define neutral regions
		calculate target stats
		fit model
			define topology


Combining scores
-------------------


Identifying regions
-------------------


Localizing signals
-------------------