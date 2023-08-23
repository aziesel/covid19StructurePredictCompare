# covid19StructurePredictCompare
Data files produced in association with our 2023 paper, Unveiling Hidden Structural Patterns in the SARS-CoV-2 Genome: Computational Insights and Comparative Analysis.


Index:

RNAz_flanked/: RNAz and LocARNA analyses on the flanked viral genome conserved regions.

RNAz_flanked_sto/: RNAz flanked sequences, converted to Stockholm format and analyzed with CaCoFold.

genomes/: the viral genomes analyzed in this paper.

multiz_tba_files/: output of MULTIZ and TBA identifying conserved regions between viral genomes.

Li sequence parsing.ipynb: a Jupyter notebook to remove gaps from the published sequence in Li et al. (2021).

RNAdistance_utilities.ipynb: a Jupyter notebook of code snippets to calculate a variety of string distances between different predicted dot-bracket formatted structures.

RNAz_to_FA.pl: a Perl script to convert RNAz output into FASTA-formatted files for each result.

RNAz_to_sto.pl: a Perl script to convert RNAz output into Stockholm-formatted files for each result.

corona_RNAz_out: results of the first round of RNAz analysis on the conserved viral regions.

corona_RNAz_out2: results of the second round of RNAz analysis on the flanked, conserved viral regions.

corona_mafft_newick: a MAFFT-based phylogenetic tree for the genomes analyzed, in modified Newick format.

flankedseqs.maf: MAF-formatted flanked conserved regions for analysis with RNAz.

statistical_analysis_of_distances.ipynb: a Jupyter notebook to calculate differences in distances between pairs of data sets to identify statistically significant differences.
