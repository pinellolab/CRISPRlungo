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
on_target_cut_sites (space-separated list of known on-target cut sites, e.g. chr2_72933869. Translocations at these sites will be reported in a separate file. Default: None)

# or search for cut sites based on guides using Cas-OFFinder
on-target_guides (space-separated list of guides used)
PAM (pam of guides, e.g. NGG)
casoffinder_num_mismatches (number of mismatches for casoffinder (e.g. 3))
cleavage_offset (number of bp from the end of the guide where cleavage is expected. Default: -3)

# specify primer information by sequence:
primer_seq (sequence of primer)
# or specify primer information by site:
primer_site (e.g. chr2_72933869)

add_non_primer_cut_targets (True/False, if true, when generating custom targets, only sites that are translocations with the primer will be included. Default: False)

#alignment stringency
arm_min_seen_bases (number of bp that are required to be aligned to each arm of custom targets. Custom targets are given by two arms from different cut sites. Default: 5)
arm_min_matched_start_bases (number of bp that are required to match exactly at the start of each arm for a valid alignment. Default: 5)

#deduplication
dedup_input (True/False, if True, reads will be deduplicated based on UMI before processing. Default: False)
dedup_input_based_on_aln_pos_and_UMI (True/False, if True, reads will be duplicated based on UMI and alignment position. This will be performed after alignment. Default: True)
umi_regex (Regex to match valid UMIs -- default: NNNWNNWNN)

#fragmentation step
fragment_size (size of fragments to create for unaligned reads, e.g. 20)
fragment_step_size (step between fragments for unaligned reads, e.g. 10)

#CRISPResso analysis
run_crispresso_genome_sites (True/False, if True, CRISPResso will be run on frequently-aligned genomic locations. If False, CRISPResso will be run only on cut sites and custom targets. Default: False)
crispresso_min_aln_score (Minimum alignment score for inclusion in CRISPResso analysis. Default: 20)
crispresso_min_count (min number of reads at a site to run crispresso2. Default: 10)

#other arguments
root (root file path to write output at)

n_processes (the number of processes to spawn for parallelization)

alignment_extension (the number of bp to extend beyond the average read length to create custom target seqeuences. Default: 50)

#program defaults
bowtie2_command: command to run bowtie2 (Default: bowtie2)
casoffinder_command: command to run casoffinder (Default: cas-offinder)
crispresso_command: command to run CRISPResso (Default: CRISPResso)
samtools_command: command to run samtools (Default: samtools)

```
A tab separates the setting name from the setting value(s).
