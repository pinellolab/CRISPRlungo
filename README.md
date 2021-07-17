# CRISPRlungo
Genome Editing Analysis of Long UNidirectional sequencing for GenOme Rearrangements


Input is provided via a settings file with entries:

```
# general required input
fastq_r1 (input fastq, e.g. test.fastq)
fastq_r2 (optional r2 input)
bowtie2_genome (genome fasta file, e.g. /REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa, must be bowtie2-indexed)

# specify cut sites
cut_sites (space-separated list of known on- and off-target cut sites, e.g. chr2_72933869 chr5_45358965)
# or search for cut sites based on guides using Cas-OFFinder
on-target_guides (space-separated list of guides used)
PAM (pam of guides, e.g. NGG)
casoffinder_num_mismatches (number of mismatches for casoffinder (e.g. 3))

#deduplication
dedup_input (True/False, if True, reads will be deduplicated based on UMI before processing. Default: False)
umi_regex (Regex to match valid UMIs -- default: NNNWNNWNN)


#optional arguments
primer_site (location of primer site, e.g. chr2_72933850)

fragment_size (size of fragments to create for unaligned reads, e.g. 20)
fragment_step_size (step between fragments for unaligned reads, e.g. 10)

root (root file path to write output at)

crispresso2_min_count (min number to run crispresso2, e.g. 10)

n_processes (the number of processes to spawn for parallelization)

#program defaults
bowtie2_command: command to run bowtie2 (Default: bowtie2)
casoffinder_command: command to run casoffinder (Default: cas-offinder)
crispresso_command: command to run CRISPResso (Default: CRISPResso)
samtools_command: command to run samtools (Default: samtools)


```
A tab separates the setting name from the setting value(s).
