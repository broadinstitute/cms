Sample workflow
=============================================================

CMS provides a computational framework to explore signals of natural selection in diploid organisms. This section describes how to do so at an abstract level. 

Preliminary considerations
-------------------

CMS is a computational and statistical framework for exploring the evolution of populations within a species at a genomic level. To that end, the user must first provide a dataset containing genotype calls for individuals in at least one putative selected population and at least one putative 'outgroup.' CMS 2.0 is designed to be flexible with respect to the number and configuration of input populations -- that is, given input of however many populations, the user can easily calculate CMS scores for any configuration of these populations. Nonetheless, CMS still relies on the user to define these populations appropriately. 

- To determine or confirm appropriate population groupings, identify outliers, etc., we recommend that users first characterize their dataset using such methods as likeliness clustering (e.g. `STRUCTURE <http://pritchardlab.stanford.edu/structure.html>`_), principal components analysis, or phylogenetic methods (see e.g. `SNPRelate <https://github.com/zhengxwen/SNPRelate>`_)
- Each population should be randomly thinned to the same number of individuals, none of whom should be related within the past few generations.
- Larger samples are generally preferable (e.g. 50-100+ diploid individuals per population). However, depending on such factors as the landscape of recombination in the species, the quality/density of genotype data, and the extent of neutral genetic divergence between represented populations, it may be possible to leverage smaller datasets. As CMS necessitates the generation of a demographic model for the given dataset, the user is advised to use their model to generate simulated data with which to perform power estimations. In the case of the 1000 Genomes Phase 3 dataset, we make our models available HERE {LINK?}.
- Denser genotype data is generally preferable. When performing genome-wide analyses, users should consider multiple testing.

Data formatting
-------------------

CMS requires the user to provide population genetic (i.e., within-species) diversity data, including genotype phase and allele polarity. 

- If your dataset is **unphased**, you can preprocess it using a program like `Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`_ or `PLINK <https://pngu.mgh.harvard.edu/~purcell/plink/>`_. 
- The identity of the **ancestral allele** at each site is typically determined by comparison to outgroups at orthologous sites. Inferred ancestral sequence is available for a number of species through e.g. `Ensembl <http://ensembl.org>`_ via their `ftp <ftp://ftp.ensembl.org/pub/release-84/fasta/ancestral_alleles/>`_. You can use `VCFtools <https://github.com/vcftools/vcftools.github.io>`_ to populate the "AA=" section of your VCF's INFO field.
- In most cases, the user will want to provide a **genetic recombination map**. If this is unavailable, CMS will assume uniform recombination rates when calculating haplotype scores. Human recombination maps are available from the `HapMap Project <http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/>`_.
- CMS works with `TPED <http://varianttools.sourceforge.net/Format/Tped>`_ datafiles, and includes support to convert from `VCF <http://samtools.github.io/hts-specs/VCFv4.3.pdf>`_ using the command line tool scans.py. {LINK}

Demographic modeling
-------------------

CMS combines several semi-independent component tests for selection in a Bayesian or Machine Learning framework. In the former case, a demographic model for the species in question is critical in order to furnish posterior distributions of scores for said component tests under alternate hypotheses of neutrality or selection. Put otherwise: a demographic model is a (conjectural) descriptive historical account of our dataset, including population sizes and migration rates across time, that can be used to generate simulated data that 'looks like' our original dataset. We then simulate many scenarios of selection in order and calculate the distributions of component scores for adaptive, linked, and neutral variants. These distributions form the basis of our Bayesian classifier. We can also circumvent the need to define posterior score distributions by using simulated data as training data for Machine Learning implementations of CMS {LINK}. 

Our modeling framework is designed to accomodate an arbitrary number of populations in a tree of arbitrary complexity; as such, it is designed to be exploratory, allowing users to iteratively perform optimizations while visualizing the effect on model goodness-of-fit. For rigorous demographic inference (i.e., in the case of a model with known topology and tractably few parameters), users may consider programs such as `dadi <https://bitbucket.org/gutenkunstlab/dadi>`_ or `diCal <https://sourceforge.net/projects/dical2/>`_. 

Following `Schaffner et al 2005 <http://www.ncbi.nlm.nih.gov/pubmed/16251467>`_ , our framework calculates a range of population summary statistics as target values, and defines error as the Root Mean Square discrepancy between target and simulated values. These summary statistics are calculated by bootstrap estimate from user-specified putative neutral regions. For human populations, the `Neutral Regions Explorer <http://nre.cb.bscb.cornell.edu/nre/>`_ is a useful resource.

The user must specify tree topology and ranges for parameter values. These can be added and removed as desired through the script params.py (??!). After target values have been estimated and model topology defined, the user can iteratively search through subsets of parameter-space using models.py (?!?) with a masterfile specifying search input. 

Calculating selection statistics
-------------------

CMS packages a number of previously described population genetic tests for recent positive selection. Haplotype scores are calculated using `selscan <https://github.com/szpiech/selscan/>`_. 

Combining scores
-------------------

CMS 2.0 allows users to define CMS scores flexibly with respect to (i) number and identity of putative selected/neutral populations, (ii) assumed demographic model, (iii) input component scores, (iv) method of score combination. In each case the user should motivate their choices and consider how robust a putative signal of selection is to variation or arbitrariness in these factors.

Identifying regions
-------------------

CMS is motivated by the need to resolve signals of selection -- that is, to identify genetic variants that confer adaptive phenotypes. Because selective events can alter patterns of population genetic diversity across large genomic regions, we take a two-step approach to this goal: we first identify putative selected regions (using CMS, another framework, prior knowledge, etc.), and then examine each region with CMS to identify a tractable list of candidate variants for further scrutiny.

Localizing signals
-------------------

Once regions are defined, we can reapply our composite framework in order to thin our list of candidate variants for further scrutiny and prioritize those sites that have the strongest evidence of selection (or other compelling evidence, e.g. overlap with known or predicted functional elements).