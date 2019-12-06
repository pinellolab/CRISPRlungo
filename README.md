# CRISPRlungo
Genome Editing Analysis of Long UNidirectional sequencing for GenOme Rearrangements


Input is provided via a settings file with entries:


genome (genome fasta file, e.g. /REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa)

bowtie2_genome (bowtei2 genome root, e.g. /REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome)

fastq (input fastq, e.g. test.fastq)

fragment_size (size of fragments to create, e.g. 20)

fragment_step_size (step between fragments, e.g. 10)

crispresso2_min_count (min number to run crispresso2, e.g. 10)

primer_site (location of primer site, e.g. chr2_72933850)

cut_sites (location of cut sites, space-separated e.g. chr2_72933869 chr5_45358965)


A tab separates the setting name from the setting value(s).
