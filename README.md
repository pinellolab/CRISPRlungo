# CRISPRlungo
Genome Editing Analysis of Long UNidirectional sequencing for GenOme Rearrangements


Input is provided via a settings file with entries:

```
# general required input
fastq_r1 (input fastq, e.g. test.fastq)
fastq_r2 (optional r2 input)
genome (genome fasta file, e.g. /REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa, must be bowtie2-indexed)

# specify cut sites
cut_sites (space-separated list of known on- and off-target cut sites, e.g. chr2_72933869 chr5_45358965)
# or search for cut sites based on guides using Cas-OFFinder
on-target_guides (space-separated list of guides used)
PAM (pam of guides, e.g. NGG)
casoffinder_num_mismatches (number of mismatches for casoffinder (e.g. 3))

#optional arguments
primer_site (location of primer site, e.g. chr2_72933850)

fragment_size (size of fragments to create for unaligned reads, e.g. 20)
fragment_step_size (step between fragments for unaligned reads, e.g. 10)

crispresso2_min_count (min number to run crispresso2, e.g. 10)

bowtie2_threads	(number of threads to run bowtie2 with)

```
A tab separates the setting name from the setting value(s).
