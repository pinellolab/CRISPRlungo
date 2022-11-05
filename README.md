# CRISPRlungo
Genome Editing Analysis of Long UNidirectional sequencing for GenOme Rearrangements

CRISPRlungo enables the analysis of single-anchor amplicon sequencing data through quantifying complex genome editing outcomes, identifying novel cut points, and using a biologically-aware alignment method to precisely measure small insertions and deletions. 

A list of CRISPRlungo parameters can be accessed by runnning `python CRISPRlungo.py -h`. Parameters can be specified via the command line or a settings file. Parameters provided in the settings file will override settings provided on the command line in case of conflicts. 

If provided, the settings file should be passed as the first argument to CRISPRlungo. The settings file may contain comments (lines starting with a '#' character). Setting names should be tab-separated from setting values.

## Example run:
- Download the file [cuts.fa](tests/testInputs/cuts.fa)
- Prepare a Bowtie2 index of your genome or download the [test genome folder](tests/genome.chr11_5225364_5225863)
- In the same directory, create a settings file called `settings.txt` with the contents:
```
fastq_r1        input.fq
genome  {path to genome}/Bowtie2Index/genome.fa or genome.chr11_5225364_5225863/genome.fa
guide_sequences CTTAGGGAACAAAGGAACCT
n_processes     20
primer_seq      TTGCAATGAAAATAAATGTTT
```
Run CRISPRlungo with the command `python CRISPRlungo.py settings.txt`. View the example output by viewing the settings.txt.CRISPRlungo.html file.


## CRISPRlungo Parameters

### Positional arguments:
```
  settings_file         Tab-separated settings file (default: None)
```

### Optional arguments:
```
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --debug               Print debug output (default: False)
  --root ROOT, --name ROOT
                        Output directory file root (default: None)
  --keep_intermediate   If true, intermediate files are not deleted (default:
                        False)
  --write_discarded_read_info
                        If true, a file with information for discarded reads
                        is produced (default: False)
  --suppress_plots      If true, no plotting will be performed (default:
                        False)
  --guide_sequences [GUIDE_SEQUENCES ...]
                        Spacer sequences of guides (multiple guide sequences
                        are separated by spaces). Spacer sequences must be
                        provided without the PAM sequence, but oriented so the
                        PAM would immediately follow the provided spacer
                        sequence (default: [])
  --cuts [CUTS ...], --cut_sites [CUTS ...]
                        Cut sites in the form chr1:234 (multiple cuts are
                        separated by spaces) (default: [])
  --on_target_cut_sites [ON_TARGET_CUT_SITES ...]
                        On-target cut sites in the form chr1:234 (multiple
                        cuts are separated by spaces) (default: [])
  --cut_classification_annotations [CUT_CLASSIFICATION_ANNOTATIONS ...]
                        User-customizable annotations for cut products in the
                        form: chr1:234:left:Custom_label (multiple annotations
                        are separated by spaces) (default: [])
  --cleavage_offset CLEAVAGE_OFFSET
                        Position where cleavage occurs, for in-silico off-
                        target search (relative to end of spacer seq -- for
                        Cas9 this is -3) (default: -3)
  --genome GENOME       Genome sequence file for alignment. This should point
                        to a file ending in ".fa", and the accompanying index
                        file (".fai") should exist. (default: None)
  --bowtie2_genome BOWTIE2_GENOME
                        Bowtie2-indexed genome file. (default: None)
  --fastq_r1 FASTQ_R1   Input fastq r1 file. Reads in this file are primed
                        from the provided primer sequence (default: None)
  --fastq_r2 FASTQ_R2   Input fastq r2 file (default: None)
  --fastq_umi FASTQ_UMI
                        Input fastq umi file (default: None)
  --novel_cut_merge_distance NOVEL_CUT_MERGE_DISTANCE
                        Novel cut sites discovered within this distance (bp)
                        from each other (and not within
                        homology_cut_merge_distance to a known/provided cut
                        site or a site with homology to guide_sequences) will
                        be merged into a single cut site. Variation in the cut
                        sites or in the fragments produced may produce
                        clusters of cut sites in a certain region. This
                        parameter will merge novel cut sites within this
                        distance into a single cut site. (default: 100)
  --homology_cut_merge_distance HOMOLOGY_CUT_MERGE_DISTANCE
                        Novel cut sites discovered within this distance (bp)
                        with a known/provided/homologous site will be merged
                        to that site. Homologous sites are defined as those
                        that have homology to guide_sequences. Novel cut sites
                        farther than homology_cut_merge_distance will be
                        merged into novel cut sites based on the parameter
                        novel_cut_merge_distance. (default: 10000)
```

#### In silico off-target search parameters:
```
  --PAM PAM             PAM for in-silico off-target search (default: None)
  --casoffinder_num_mismatches CASOFFINDER_NUM_MISMATCHES
                        If greater than zero, the number of Cas-OFFinder
                        mismatches for in-silico off-target search. If this
                        value is zero, Cas-OFFinder is not run (default: 0)
```

#### Primer and filtering parameters and settings:
```
  --primer_seq PRIMER_SEQ
                        Sequence of primer (default: None)
  --primer_in_r2        If true, the primer is in R2. By default, the primer
                        is required to be present in R1.
                        (default: False)
  --min_primer_aln_score MIN_PRIMER_ALN_SCORE
                        Minimum primer/origin alignment score for trimming.
                        (default: 40)
  --min_primer_length MIN_PRIMER_LENGTH
                        Minimum length of sequence required to match between
                        the primer/origin and read sequence (default: 30)
  --min_read_length MIN_READ_LENGTH
                        Minimum length of read after all filtering (default:
                        30)
  --transposase_adapter_seq TRANSPOSASE_ADAPTER_SEQ
                        Transposase adapter sequence to be trimmed from reads
                        (default: CTGTCTCTTATACACATCTGACGCTGCCGACGA)
```

#### Alignment cutoff parameters:
```
  --arm_min_matched_start_bases ARM_MIN_MATCHED_START_BASES
                        Number of bases that are required to be matching (no
                        indels or mismatches) at the beginning of the read on
                        each "side" of the alignment. E.g. if
                        arm_min_matched_start_bases is set to 5, the first and
                        last 5bp of the read alignment would have to match
                        exactly to the aligned location. (default: 10)
  --arm_max_clipped_bases ARM_MAX_CLIPPED_BASES
                        Maximum number of clipped bases at the beginning of
                        the alignment. Bowtie2 alignment marks reads on the
                        beginning or end of the read as "clipped" if they do
                        not align to the genome. This could arise from CRISPR-
                        induced insertions, or bad alignments. We would expect
                        to see clipped bases only on one side. This parameter
                        sets the threshold for clipped bases on both sides of
                        the read. E.g. if arm_max_clipped_bases is 0, read
                        alignments with more than 0bp on the right AND left
                        side of the alignment would be discarded. An alignment
                        with 5bp clipped on the left and 0bp clipped on the
                        right would be accepted. An alignment with 5bp clipped
                        on the left and 3bp clipped on the right would be
                        discarded. (default: 0)
  --ignore_n            If set, "N" bases will be ignored. By default (False)
                        N bases will count as mismatches in the number of
                        bases required to match at each arm/side of the read
                        (default: False)
  --discard_reads_with_poor_alignment
                        If set, reads with poor alignment (fewer than
                        --arm_min_matched_start_bases mismatches at the
                        alignment ends or more than --arm_max_clipped_bases on
                        both sides of the read) are discarded from final
                        analysis and counts (default: False)
```
#### CRISPResso settings:
```
  --crispresso_min_count CRISPRESSO_MIN_COUNT
                        Min number of reads required to be seen at a site for
                        it to be analyzed by CRISPResso (default: 50)
  --crispresso_max_indel_size CRISPRESSO_MAX_INDEL_SIZE
                        Maximum length of indel (as determined by genome
                        alignment) for a read to be analyzed by CRISPResso.
                        Reads with indels longer than this length will not be
                        analyzed by CRISPResso, but the indel length will be
                        reported elsewhere. (default: 50)
  --crispresso_min_aln_score CRISPRESSO_MIN_ALN_SCORE
                        Min alignment score to reference sequence for
                        quantification by CRISPResso (default: 20)
  --crispresso_quant_window_size CRISPRESSO_QUANT_WINDOW_SIZE
                        Number of bp on each side of a cut to consider for
                        edits (default: 1)
  --run_crispresso_on_novel_sites
                        If set, CRISPResso analysis will be performed on novel
                        cut sites. If false, CRISPResso analysis will only be
                        performed on user-provided on- and off-targets
                        (default: False)

```
#### Pipeline parameters:
```
  --cutadapt_command CUTADAPT_COMMAND
                        Command to run cutadapt (default: cutadapt)
  --samtools_command SAMTOOLS_COMMAND
                        Command to run samtools (default: samtools)
  --bowtie2_command BOWTIE2_COMMAND
                        Command to run bowtie2 (default: bowtie2)
  --crispresso_command CRISPRESSO_COMMAND
                        Command to run crispresso (default: CRISPResso)
  --casoffinder_command CASOFFINDER_COMMAND
                        Command to run casoffinder (default: cas-offinder)
  --n_processes N_PROCESSES
                        Number of processes to run on (may be set to "max")
                        (default: 1)
```
#### UMI parameters:
```
  --dedup_input_on_UMI  If set, input reads will be deduplicated based on UMI
                        before alignment (default: False)
  --dedup_input_based_on_aln_pos_and_UMI
                        If set, perform deduplication based on alignment
                        position and UMI. (default: False)
  --umi_regex UMI_REGEX
                        String specifying regex that UMI must match (default:
                        NNWNNWNNN)
```
#### R1/R2 support settings:
``` 
  --r1_r2_support_max_distance R1_R2_SUPPORT_MAX_DISTANCE
                        Max distance between r1 and r2 for the read pair to be
                        classified as "supported" by r2 (default: 10000)
  --discard_reads_without_r2_support
                        If set, reads without r2 support will be discarded
                        from final analysis and counts (default: False)


```
