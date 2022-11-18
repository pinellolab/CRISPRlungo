import argparse
from collections import defaultdict
import ctypes
import gzip
import json
import logging
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle,Patch
import multiprocessing as mp
import math
import numpy as np
import os
import re
import subprocess
import sys
from lib import ssw_lib
from lib.pyCircos import pyCircos_lib as pc
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPResso2Align

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

__version__ = "v0.1.9"

def processCRISPRlungo(settings):
    """Run the CRISPRlungo pipeline

    Args:
        settings (dict): Dictionary of setting_name->setting_value to use for the run.
    """
    logger = logging.getLogger('CRISPRlungo')
    final_file = settings['root'] + '.summary.txt'
    if settings['can_use_previous_analysis'] and os.path.exists(final_file):
        logger.info('Successfully completed!')
        return final_file

    #data structures for plots for report
    summary_plot_objects = []  # list of PlotObjects for plotting

    assert_dependencies(
            cutadapt_command=settings['cutadapt_command'],
            samtools_command=settings['samtools_command'],
            bowtie2_command=settings['bowtie2_command'],
            crispresso_command=settings['crispresso_command'],
            casoffinder_command=settings['casoffinder_command']
            )

    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, primer_is_genomic, av_read_length, num_reads_input = prep_input(
            root = settings['root']+'.primerInfo',
            primer_seq = settings['primer_seq'],
            min_primer_length = settings['min_primer_length'],
            guide_seqs = settings['guide_sequences'],
            cleavage_offset = settings['cleavage_offset'],
            fastq_r1 = settings['fastq_r1'],
            samtools_command = settings['samtools_command'],
            genome = settings['genome'],
            bowtie2_command = settings['bowtie2_command'],
            bowtie2_genome = settings['bowtie2_genome'],
            can_use_previous_analysis = settings['can_use_previous_analysis']
            )
    logger.info('%d reads in input'%num_reads_input)

    if settings['PAM'] is not None and settings['casoffinder_num_mismatches'] > 0:
        casoffinder_cut_sites,casoffinder_cut_annotations = get_cut_sites_casoffinder(
                root = settings['root']+'.casoffinder',
                genome=settings['genome'],
                pam=settings['PAM'],
                guides=settings['guide_sequences'],
                cleavage_offset=settings['cleavage_offset'],
                num_mismatches=settings['casoffinder_num_mismatches'],
                casoffinder_command=settings['casoffinder_command'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )
        for cut_site in casoffinder_cut_sites:
            if cut_site not in cut_sites:
                cut_sites.append(cut_site)
                cut_annotations[cut_site] = casoffinder_cut_annotations[cut_site]


    curr_r1_file = settings['fastq_r1'] #keep track of current input files (through UMI adding, etc.)
    curr_r2_file = settings['fastq_r2']

    if settings['primer_in_r2']:
        curr_r1_file = settings['fastq_r2']
        curr_r2_file = settings['fastq_r1']

    #if umis are provided, add them to the fastqs
    if settings['fastq_umi']:
        umi_r1, umi_r2 = add_umi_from_umi_file(
            root = settings['root']+'.addUMI',
            fastq_r1 = curr_r1_file,
            fastq_r2 = curr_r2_file,
            fastq_umi = settings['fastq_umi'],
            can_use_previous_analysis = settings['can_use_previous_analysis']
        )
        curr_r1_file = umi_r1
        curr_r2_file = umi_r2

    post_dedup_count = num_reads_input
    if settings['dedup_input_on_UMI']:
        (dedup_r1, dedup_r2, dedup_tot_read_count, dedup_count_with_regex, post_dedup_count, post_dedup_read_count
                ) = dedup_input_file(
                        root = settings['root']+'.dedup',
                        fastq_r1 = curr_r1_file,
                        fastq_r2 = curr_r2_file,
                        umi_regex = settings['umi_regex'],
                        can_use_previous_analysis = settings['can_use_previous_analysis']
        )
        curr_r1_file = dedup_r1
        curr_r2_file = dedup_r2

    filtered_on_primer_fastq_r1, filtered_on_primer_fastq_r2, post_filter_on_primer_read_count, filter_on_primer_plot_obj = filter_on_primer(
            root = settings['root']+'.trimPrimers',
            fastq_r1 = curr_r1_file,
            fastq_r2 = curr_r2_file,
            origin_seq = origin_seq,
            min_primer_aln_score = settings['min_primer_aln_score'],
            min_primer_length = settings['min_primer_length'],
            min_read_length = settings['min_read_length'],
            transposase_adapter_seq = settings['transposase_adapter_seq'],
            n_processes = settings['n_processes'],
            cutadapt_command = settings['cutadapt_command'],
            keep_intermediate = settings['keep_intermediate'],
            suppress_plots = settings['suppress_plots'],
            can_use_previous_analysis = settings['can_use_previous_analysis'],
            )

    if filter_on_primer_plot_obj is not None:
        filter_on_primer_plot_obj.order = 3
        summary_plot_objects.append(filter_on_primer_plot_obj)

    #perform alignment
    (genome_mapped_bam
            ) = align_reads(
                root = settings['root']+'.genomeAlignment',
                fastq_r1 = filtered_on_primer_fastq_r1,
                fastq_r2 = filtered_on_primer_fastq_r2,
                bowtie2_reference = settings['bowtie2_genome'],
                bowtie2_command = settings['bowtie2_command'],
                bowtie2_threads = settings['n_processes'],
                samtools_command = settings['samtools_command'],
                keep_intermediate = settings['keep_intermediate'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )
    if not settings['keep_intermediate']:
        if os.path.exists(filtered_on_primer_fastq_r1):
            logger.debug('Deleting intermediate file ' + filtered_on_primer_fastq_r1)
            os.remove(filtered_on_primer_fastq_r1)
        if filtered_on_primer_fastq_r2 is not None and os.path.exists(filtered_on_primer_fastq_r2):
            logger.debug('Deleting intermediate file ' + filtered_on_primer_fastq_r2)
            os.remove(filtered_on_primer_fastq_r2)

    (final_assignment_file,final_read_ids_for_crispresso,final_cut_counts,cut_classification_lookup, final_read_count, discarded_read_counts, classification_read_counts, classification_indel_read_counts,
        chr_aln_plot_obj,tlen_plot_obj,deduplication_plot_obj,tx_order_plot_obj,tx_count_plot_obj,tx_circos_plot_obj,classification_plot_obj,classification_indel_plot_obj,
        origin_indel_hist_plot_obj,origin_inversion_hist_plot_obj,origin_deletion_hist_plot_obj,origin_indel_depth_plot_obj,
        r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,discarded_reads_plot_obj) = make_final_read_assignments(
                root = settings['root']+'.final',
                genome_mapped_bam = genome_mapped_bam,
                origin_seq = origin_seq,
                cut_sites = cut_sites,
                cut_annotations = cut_annotations,
                cut_classification_annotations = settings['cut_classification_annotations'],
                guide_seqs = settings['guide_sequences'],
                cleavage_offset = settings['cleavage_offset'],
                min_primer_length = settings['min_primer_length'],
                genome = settings['genome'],
                r1_r2_support_max_distance = settings['r1_r2_support_max_distance'],
                novel_cut_merge_distance = settings['novel_cut_merge_distance'],
                known_cut_merge_distance = settings['known_cut_merge_distance'],
                origin_cut_merge_distance = settings['origin_cut_merge_distance'],
                arm_min_matched_start_bases = settings['arm_min_matched_start_bases'],
                arm_max_clipped_bases = settings['arm_max_clipped_bases'],
                genome_map_resolution = 1000000,
                crispresso_max_indel_size = settings['crispresso_max_indel_size'],
                suppress_dedup_on_aln_pos_and_UMI_filter = settings['suppress_dedup_on_aln_pos_and_UMI_filter'],
                suppress_r2_support_filter = settings['suppress_r2_support_filter'],
                suppress_poor_alignment_filter = settings['suppress_poor_alignment_filter'],
                write_discarded_read_info = settings['write_discarded_read_info'],
                samtools_command = settings['samtools_command'],
                keep_intermediate = settings['keep_intermediate'],
                suppress_plots = settings['suppress_plots'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )

    if chr_aln_plot_obj is not None:
        chr_aln_plot_obj.order=20
        summary_plot_objects.append(chr_aln_plot_obj)

    if tlen_plot_obj is not None:
        tlen_plot_obj.order=15
        summary_plot_objects.append(tlen_plot_obj)

    if tlen_plot_obj is not None:
        deduplication_plot_obj.order=16
        summary_plot_objects.append(deduplication_plot_obj)

    if classification_plot_obj is not None:
        classification_plot_obj.order=36
        summary_plot_objects.append(classification_plot_obj)

    if classification_indel_plot_obj is not None:
        classification_indel_plot_obj.order=2
        summary_plot_objects.append(classification_indel_plot_obj)

    if tx_order_plot_obj is not None:
        tx_order_plot_obj.order= 38
        summary_plot_objects.append(tx_order_plot_obj)

    if tx_count_plot_obj is not None:
        tx_count_plot_obj.order= 40
        summary_plot_objects.append(tx_count_plot_obj)

    if tx_circos_plot_obj is not None:
        tx_circos_plot_obj.order= 41
        summary_plot_objects.append(tx_circos_plot_obj)

    if origin_indel_hist_plot_obj is not None:
        origin_indel_hist_plot_obj.order=25
        summary_plot_objects.append(origin_indel_hist_plot_obj)

    if origin_inversion_hist_plot_obj is not None:
        origin_inversion_hist_plot_obj.order=26
        summary_plot_objects.append(origin_inversion_hist_plot_obj)

    if origin_deletion_hist_plot_obj is not None:
        origin_deletion_hist_plot_obj.order=27
        summary_plot_objects.append(origin_deletion_hist_plot_obj)

    if origin_indel_depth_plot_obj is not None:
        origin_indel_depth_plot_obj.order=28
        summary_plot_objects.append(origin_indel_depth_plot_obj)

    if r1_r2_support_plot_obj is not None:
        r1_r2_support_plot_obj.order=29
        summary_plot_objects.append(r1_r2_support_plot_obj)
    if r1_r2_support_dist_plot_obj is not None:
        r1_r2_support_dist_plot_obj.order=30
        summary_plot_objects.append(r1_r2_support_dist_plot_obj)

    if discarded_reads_plot_obj is not None:
        discarded_reads_plot_obj.order=3
        summary_plot_objects.append(discarded_reads_plot_obj)

    #crispresso_infos: dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command
    crispresso_infos = prep_crispresso2(
                root = settings['root']+'.CRISPResso_r1',
                input_fastq_file = curr_r1_file,
                read_ids_for_crispresso = final_read_ids_for_crispresso,
                cut_counts_for_crispresso = final_cut_counts,
                cut_classification_lookup = cut_classification_lookup,
                cut_annotations = cut_annotations,
                av_read_length = av_read_length,
                origin_seq = origin_seq,
                cleavage_offset = settings['cleavage_offset'],
                genome = settings['genome'],
                genome_len_file = settings['genome']+'.fai',
                crispresso_min_count = settings['crispresso_min_count'],
                crispresso_min_aln_score = settings['crispresso_min_aln_score'],
                crispresso_quant_window_size = settings['crispresso_quant_window_size'],
                run_crispresso_on_novel_sites = settings['run_crispresso_on_novel_sites'],
                samtools_command=settings['samtools_command'],
                crispresso_command=settings['crispresso_command'],
                n_processes = settings['n_processes']
                )
    crispresso_results, crispresso_classification_plot_obj = run_and_aggregate_crispresso(
                root = settings['root']+".CRISPResso",
                crispresso_infos = crispresso_infos,
                final_assignment_file = final_assignment_file,
                n_processes = settings['n_processes'],
                keep_intermediate = settings['keep_intermediate'],
                suppress_plots = settings['suppress_plots'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )
    if crispresso_classification_plot_obj is not None:
        crispresso_classification_plot_obj.order=40
        summary_plot_objects.append(crispresso_classification_plot_obj)

    final_summary_file, final_summary_plot_obj = make_final_summary(
            root = settings['root']+'.summary',
            num_reads_input = num_reads_input,
            post_dedup_count = post_dedup_count,
            post_filter_on_primer_read_count = post_filter_on_primer_read_count,
            final_read_count = final_read_count,
            discarded_read_counts = discarded_read_counts,
            classification_read_counts = classification_read_counts,
            classification_indel_read_counts = classification_indel_read_counts,
            suppress_plots = settings['suppress_plots']
            )

    if final_summary_plot_obj is not None:
        final_summary_plot_obj.order= 1
        summary_plot_objects.append(final_summary_plot_obj)

    if not settings['suppress_plots']:
        make_report(report_file=settings['root']+".html",
                report_name = 'Report',
                crisprlungo_folder = '',
                crispresso_run_names = crispresso_results['run_names'],
                crispresso_sub_html_files = crispresso_results['run_sub_htmls'],
                summary_plot_objects = summary_plot_objects,
                )

    logger.info('Successfully completed!')
    return final_summary_file

    # FINISHED

def parse_settings(args):
    """
    Parses settings from the command line
        First parses from the settings file, then parses from command line

    Args:
        args: command line arguments

    Returns:
        settings: dict of parsed settings
    """
    parser = argparse.ArgumentParser(description='CRISPRlungo: Analyzing unidirectional sequencing of genome editing', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version="%(prog)s "+__version__)
    parser.add_argument('settings_file', nargs='*', help='Tab-separated settings file')

    parser.add_argument('--debug', action='store_true', help='Print debug output')
    parser.add_argument('--root','--name', type=str, default=None, help='Output directory file root')
    parser.add_argument('--keep_intermediate',action='store_true',help='If true, intermediate files are not deleted')
    parser.add_argument('--write_discarded_read_info',action='store_true',help='If true, a file with information for discarded reads is produced')
    parser.add_argument('--suppress_plots',action='store_true',help='If true, no plotting will be performed')


    parser.add_argument('--guide_sequences', nargs='+', help='Spacer sequences of guides (multiple guide sequences are separated by spaces). Spacer sequences must be provided without the PAM sequence, but oriented so the PAM would immediately follow the provided spacer sequence', action='extend', default=[])
    parser.add_argument('--cut_classification_annotations', nargs='+', help='User-customizable annotations for cut products in the form: chr1:234:left:Custom_label (multiple annotations are separated by spaces)', action='extend', default=[])
    parser.add_argument('--cleavage_offset', type=int, help='Position where cleavage occurs, for in-silico off-target search (relative to end of spacer seq -- for Cas9 this is -3)', default=-3)

    parser.add_argument('--genome', help='Genome sequence file for alignment. This should point to a file ending in ".fa", and the accompanying index file (".fai") should exist.', default=None)
    parser.add_argument('--bowtie2_genome', help='Bowtie2-indexed genome file.',default=None)

    parser.add_argument('--fastq_r1', help='Input fastq r1 file. Reads in this file are primed from the provided primer sequence', default=None)
    parser.add_argument('--fastq_r2', help='Input fastq r2 file', default=None)
    parser.add_argument('--fastq_umi', help='Input fastq umi file', default=None)

    parser.add_argument('--novel_cut_merge_distance', type=int, help='Novel cut sites discovered within this distance (bp) from each other (and not within known_cut_merge_distance to a known/provided cut site or a site with homology to guide_sequences) will be merged into a single cut site. Variation in the cut sites or in the fragments produced may produce clusters of cut sites in a certain region. This parameter will merge novel cut sites within this distance into a single cut site.', default=50)
    parser.add_argument('--known_cut_merge_distance', type=int, help='Novel cut sites discovered within this distance (bp) with a known/provided/homologous site (that is not the origin) will be merged to that site. Homologous sites are defined as those that have homology to guide_sequences. Novel cut sites farther than known_cut_merge_distance will be merged into novel cut sites based on the parameter novel_cut_merge_distance.', default=50)
    parser.add_argument('--origin_cut_merge_distance', type=int, help='Reads aligned within this distance (bp) to the origin site will be merged to that origin.', default=10000)
    parser.add_argument('--short_indel_length_cutoff', type=int, help='For reads aligned to a cut site, indels this size or shorter are classified as "short indels" while indels larger than this size are classified as "long indels" ', default=50)

    #for finding offtargets with casoffinder
    ot_group = parser.add_argument_group('In silico off-target search parameters')
    ot_group.add_argument('--PAM', type=str, help='PAM for in-silico off-target search', default=None)
    ot_group.add_argument('--casoffinder_num_mismatches', type=int, help='If greater than zero, the number of Cas-OFFinder mismatches for in-silico off-target search. If this value is zero, Cas-OFFinder is not run', default=0)

    #specify primer filtering information
    p_group = parser.add_argument_group('Primer and filtering parameters and settings')
    p_group.add_argument('--primer_seq', type=str, help='Sequence of primer',default=None)
    p_group.add_argument('--primer_in_r2', help='If true, the primer is in R2. By default, the primer is required to be in R1.',action='store_true')
    p_group.add_argument('--min_primer_aln_score', type=int, help='Minimum primer/origin alignment score for trimming.',default=40)
    p_group.add_argument('--min_primer_length', type=int, help='Minimum length of sequence required to match between the primer/origin and read sequence',default=30)
    p_group.add_argument('--min_read_length', type=int, help='Minimum length of read after all filtering',default=30)
    p_group.add_argument('--transposase_adapter_seq', type=str, help='Transposase adapter sequence to be trimmed from reads',default='CTGTCTCTTATACACATCTGACGCTGCCGACGA')


    #min alignment cutoffs for alignment to each arm/side of read
    a_group = parser.add_argument_group('Alignment cutoff parameters')
    a_group.add_argument('--arm_min_matched_start_bases', type=int, help='Number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each "side" of the alignment. E.g. if arm_min_matched_start_bases is set to 5, the first and last 5bp of the read alignment would have to match exactly to the aligned location.', default=10)
    a_group.add_argument('--arm_max_clipped_bases', type=int, help='Maximum number of clipped bases at the beginning of the alignment. Bowtie2 alignment marks reads on the beginning or end of the read as "clipped" if they do not align to the genome. This could arise from CRISPR-induced insertions, or bad alignments. ' + \
        'We would expect to see clipped bases only on one side. This parameter sets the threshold for clipped bases on both sides of the read.  E.g. if arm_max_clipped_bases is 0, read alignments with more than 0bp on the right AND left side of the alignment would be discarded. An alignment with 5bp clipped on the left and 0bp clipped on the right would be accepted. An alignment with 5bp clipped on the left and 3bp clipped on the right would be discarded.', default=0)
    a_group.add_argument('--ignore_n', help='If set, "N" bases will be ignored. By default (False) N bases will count as mismatches in the number of bases required to match at each arm/side of the read', action='store_true')
    a_group.add_argument('--suppress_poor_alignment_filter', help='If set, reads with poor alignment (fewer than --arm_min_matched_start_bases matches at the alignment ends or more than --arm_max_clipped_bases on both sides of the read) are included in final analysis and counts. By default they are excluded.', action='store_true')

    #CRISPResso settings
    c_group = parser.add_argument_group('CRISPResso settings')
    c_group.add_argument('--crispresso_min_count', type=int, help='Min number of reads required to be seen at a site for it to be analyzed by CRISPResso', default=50)
    c_group.add_argument('--crispresso_max_indel_size', type=int, help='Maximum length of indel (as determined by genome alignment) for a read to be analyzed by CRISPResso. Reads with indels longer than this length will not be analyzed by CRISPResso, but the indel length will be reported elsewhere.', default=50)
    c_group.add_argument('--crispresso_min_aln_score', type=int, help='Min alignment score to reference sequence for quantification by CRISPResso', default=20)
    c_group.add_argument('--crispresso_quant_window_size', type=int, help='Number of bp on each side of a cut to consider for edits', default=1)
    c_group.add_argument('--run_crispresso_on_novel_sites', help='If set, CRISPResso analysis will be performed on novel cut sites. If false, CRISPResso analysis will only be performed on user-provided on- and off-targets', action='store_true')

    #sub-command parameters
    p_group = parser.add_argument_group('Pipeline parameters')
    p_group.add_argument('--cutadapt_command', help='Command to run cutadapt', default='cutadapt')
    p_group.add_argument('--samtools_command', help='Command to run samtools', default='samtools')
    p_group.add_argument('--bowtie2_command', help='Command to run bowtie2', default='bowtie2')
    p_group.add_argument('--crispresso_command', help='Command to run crispresso', default='CRISPResso')
    p_group.add_argument('--casoffinder_command', help='Command to run casoffinder', default='cas-offinder')
    p_group.add_argument('--n_processes', type=str, help='Number of processes to run on (may be set to "max")', default='1')


    #umi settings
    u_group = parser.add_argument_group('UMI parameters')
    u_group.add_argument('--dedup_input_on_UMI', help='If set, input reads will be deduplicated based on UMI before alignment', action='store_true')
    u_group.add_argument('--suppress_dedup_on_aln_pos_and_UMI_filter', help='If set, reads that are called as deduplicates based on alignment position and UMI will be included in final analysis and counts. By default, these reads are excluded.', action='store_true')
    u_group.add_argument('--umi_regex', type=str, help='String specifying regex that UMI must match', default='NNWNNWNNN')

    #R1/R2 support settings
    r_group = parser.add_argument_group('R1/R2 support settings')
    r_group.add_argument('--r1_r2_support_max_distance', type=int, help='Max distance between r1 and r2 for the read pair to be classified as "supported" by r2', default=10000)
    r_group.add_argument('--suppress_r2_support_filter', help='If set, reads without r2 support will be included in final analysis and counts. By default these reads are excluded.', action='store_true')


    cmd_args = parser.parse_args(args[1:])

    settings = {}
    settings_file_args = {}

    # try to read in args from the settings file
    if (len(cmd_args.settings_file) > 0):
        for s_file in cmd_args.settings_file:
            with open(s_file, 'r') as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    line_els = line.split("#")[0].rstrip('\n').split("\t")
                    if line_els == ['']:
                        continue
                    if (len(line_els) < 2):
                        raise Exception('Cannot parse line "' + line + '"\nA tab must separate the key and value in the settings file')
                    key = line_els[0].strip()
                    val = line_els[1].strip()
                    if key in settings_file_args and settings_file_args[key] != val:
                        raise Exception('The parameter ' + key + ' is provided twice')
                    settings_file_args[key] = val

    settings['root'] = cmd_args.root
    if 'root' in settings_file_args:
        settings['root'] = settings_file_args['root']
        settings_file_args.pop('root')
    if settings['root'] is None:
        if len(cmd_args.settings_file) > 0:
            settings['root'] = cmd_args.settings_file[0] + ".CRISPRlungo"
        else:
            settings['root'] = "CRISPRlungo"
    #create output directory if appropriate
    output_dir = os.path.dirname(settings['root'])
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)


    settings['debug'] = cmd_args.debug
    if 'debug' in settings_file_args:
        settings['debug'] = (settings_file_args['debug'].lower() == 'true')
        settings_file_args.pop('debug')

    logger = logging.getLogger('CRISPRlungo')
    logging_level = logging.INFO
    if settings['debug']:
        logging_level=logging.DEBUG

    logger.setLevel(logging.DEBUG)

    if not logger.hasHandlers():
        log_formatter = logging.Formatter("%(asctime)s:%(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")

        sh = logging.StreamHandler()
        sh.setFormatter(log_formatter)
        sh.setLevel(logging_level)
        logger.addHandler(sh)

        fh = logging.FileHandler(settings['root']+".log")
        fh.setFormatter(log_formatter)
        logger.addHandler(fh)


    logger.info('CRISPRlungo ' + __version__)

    logger.info('Parsing settings file')

    settings['keep_intermediate'] = cmd_args.keep_intermediate
    if 'keep_intermediate' in settings_file_args:
        settings['keep_intermediate'] = (settings_file_args['keep_intermediate'].lower() == 'true')
        settings_file_args.pop('keep_intermediate')

    settings['write_discarded_read_info'] = cmd_args.write_discarded_read_info
    if 'write_discarded_read_info' in settings_file_args:
        settings['write_discarded_read_info'] = (settings_file_args['write_discarded_read_info'].lower() == 'true')
        settings_file_args.pop('write_discarded_read_info')

    settings['suppress_plots'] = cmd_args.suppress_plots
    if 'suppress_plots' in settings_file_args:
        settings['suppress_plots'] = (settings_file_args['suppress_plots'].lower() == 'true')
        settings_file_args.pop('suppress_plots')

    settings['cut_classification_annotations'] = []
    annotations_arr = cmd_args.cut_classification_annotations
    if 'cut_classification_annotations' in settings_file_args:
        annotations_arr = settings_file_args['cut_classification_annotations'].split(" ")
        if "[" in settings_file_args['cut_classification_annotations']:
            annotations_arr = [x.strip() for x in settings_file_args['cut_classification_annotations'].strip("[] ").replace('"','').replace("'","").split(",")]
        settings_file_args.pop('cut_classification_annotations')

    for arg in annotations_arr:
        arg_els = arg.split(":")
        if len(arg_els) != 4: # chr:pos:direction:annotation
            parser.print_usage()
            raise Exception('Error: cut classification annotations must contain 4 elements (e.g. chr1:234:left:my_annotation). Please check the value "' + arg + '"')
        settings['cut_classification_annotations'].append(arg.replace("__"," ")) #replace __ with space

    settings['PAM'] = cmd_args.PAM
    if 'PAM' in settings_file_args:
        settings['PAM'] = settings_file_args['PAM']
        settings_file_args.pop('PAM')

    guide_arr = cmd_args.guide_sequences
    settings['guide_sequences'] = []
    if 'guide_sequences' in settings_file_args:
        guide_arr = settings_file_args['guide_sequences'].split(" ")

        if "[" in settings_file_args['guide_sequences']:
            guide_arr = [x.strip() for x in settings_file_args['guide_sequences'].strip("[] ").replace('"','').replace("'","").split(",")]
        settings_file_args.pop('guide_sequences')

    for guide in guide_arr:
        wrong_nt=CRISPRessoShared.find_wrong_nt(guide)
        if wrong_nt:
            raise Exception('The guide sequence %s contains bad characters: %s'%(guide,'"'+'", "'.join(wrong_nt)+'"'))
        settings['guide_sequences'].append(guide)


    settings['casoffinder_num_mismatches'] = cmd_args.casoffinder_num_mismatches
    if 'casoffinder_num_mismatches' in settings_file_args:
        settings['casoffinder_num_mismatches'] = int(settings_file_args['casoffinder_num_mismatches'])
        settings_file_args.pop('casoffinder_num_mismatches')

    settings['cleavage_offset'] = cmd_args.cleavage_offset
    if 'cleavage_offset' in settings_file_args:
        settings['cleavage_offset'] = int(settings_file_args['cleavage_offset'])
        settings_file_args.pop('cleavage_offset')

    settings['primer_seq'] = cmd_args.primer_seq
    if 'primer_seq' in settings_file_args:
        settings['primer_seq'] = settings_file_args['primer_seq']
        settings_file_args.pop('primer_seq')

    settings['primer_in_r2'] = cmd_args.primer_in_r2
    if 'primer_in_r2' in settings_file_args:
        settings['primer_in_r2'] = (settings_file_args['primer_in_r2'].lower() == 'true')
        settings_file_args.pop('primer_in_r2')

    settings['min_primer_aln_score'] = cmd_args.min_primer_aln_score
    if 'min_primer_aln_score' in settings_file_args:
        settings['min_primer_aln_score'] = int(settings_file_args['min_primer_aln_score'])
        settings_file_args.pop('min_primer_aln_score')

    settings['min_primer_length'] = cmd_args.min_primer_length
    if 'min_primer_length' in settings_file_args:
        settings['min_primer_length'] = int(settings_file_args['min_primer_length'])
        settings_file_args.pop('min_primer_length')

    settings['min_read_length'] = cmd_args.min_read_length
    if 'min_read_length' in settings_file_args:
        settings['min_read_length'] = int(settings_file_args['min_read_length'])
        settings_file_args.pop('min_read_length')

    settings['transposase_adapter_seq'] = cmd_args.transposase_adapter_seq
    if 'transposase_adapter_seq' in settings_file_args:
        settings['transposase_adapter_seq'] = settings_file_args['transposase_adapter_seq']
        settings_file_args.pop('transposase_adapter_seq')

    settings['arm_min_matched_start_bases'] = cmd_args.arm_min_matched_start_bases
    if 'arm_min_matched_start_bases' in settings_file_args:
        settings['arm_min_matched_start_bases'] = int(settings_file_args['arm_min_matched_start_bases'])
        settings_file_args.pop('arm_min_matched_start_bases')

    settings['arm_max_clipped_bases'] = cmd_args.arm_max_clipped_bases
    if 'arm_max_clipped_bases' in settings_file_args:
        settings['arm_max_clipped_bases'] = int(settings_file_args['arm_max_clipped_bases'])
        settings_file_args.pop('arm_max_clipped_bases')

    settings['ignore_n'] = cmd_args.ignore_n
    if 'ignore_n' in settings_file_args:
        settings['ignore_n'] = (settings_file_args['ignore_n'].lower() == 'true')
        settings_file_args.pop('ignore_n')

    settings['suppress_poor_alignment_filter'] = cmd_args.suppress_poor_alignment_filter
    if 'suppress_poor_alignment_filter' in settings_file_args:
        settings['suppress_poor_alignment_filter'] = (settings_file_args['suppress_poor_alignment_filter'].lower() == 'true')
        settings_file_args.pop('suppress_poor_alignment_filter')

    settings['run_crispresso_on_novel_sites'] = cmd_args.run_crispresso_on_novel_sites
    if 'run_crispresso_on_novel_sites' in settings_file_args:
        settings['run_crispresso_on_novel_sites'] = (settings_file_args['run_crispresso_on_novel_sites'].lower() == 'true')
        settings_file_args.pop('run_crispresso_on_novel_sites')

    settings['crispresso_min_count'] = cmd_args.crispresso_min_count
    if 'crispresso_min_count' in settings_file_args:
        settings['crispresso_min_count'] = int(settings_file_args['crispresso_min_count'])
        settings_file_args.pop('crispresso_min_count')

    settings['crispresso_max_indel_size'] = cmd_args.crispresso_max_indel_size
    if 'crispresso_max_indel_size' in settings_file_args:
        settings['crispresso_max_indel_size'] = int(settings_file_args['crispresso_max_indel_size'])
        settings_file_args.pop('crispresso_max_indel_size')

    settings['crispresso_min_aln_score'] = cmd_args.crispresso_min_aln_score
    if 'crispresso_min_aln_score' in settings_file_args:
        settings['crispresso_min_aln_score'] = int(settings_file_args['crispresso_min_aln_score'])
        settings_file_args.pop('crispresso_min_aln_score')

    settings['crispresso_quant_window_size'] = cmd_args.crispresso_quant_window_size
    if 'crispresso_quant_window_size' in settings_file_args:
        settings['crispresso_quant_window_size'] = int(settings_file_args['crispresso_quant_window_size'])
        settings_file_args.pop('crispresso_quant_window_size')

    settings['cutadapt_command'] = cmd_args.cutadapt_command
    if 'cutadapt_command' in settings_file_args:
        settings['cutadapt_command'] = settings_file_args['cutadapt_command']
        settings_file_args.pop('cutadapt_command')

    settings['samtools_command'] = cmd_args.samtools_command
    if 'samtools_command' in settings_file_args:
        settings['samtools_command'] = settings_file_args['samtools_command']
        settings_file_args.pop('samtools_command')

    settings['bowtie2_command'] = cmd_args.bowtie2_command
    if 'bowtie2_command' in settings_file_args:
        settings['bowtie2_command'] = settings_file_args['bowtie2_command']
        settings_file_args.pop('bowtie2_command')

    settings['crispresso_command'] = cmd_args.crispresso_command
    if 'crispresso_command' in settings_file_args:
        settings['crispresso_command'] = settings_file_args['crispresso_command']
        settings_file_args.pop('crispresso_command')

    settings['casoffinder_command'] = cmd_args.casoffinder_command
    if 'casoffinder_command' in settings_file_args:
        settings['casoffinder_command'] = settings_file_args['casoffinder_command']
        settings_file_args.pop('casoffinder_command')

    settings['n_processes'] = cmd_args.n_processes
    if 'n_processes' in settings_file_args:
        settings['n_processes'] = settings_file_args['n_processes']
        settings_file_args.pop('n_processes')

    if settings['n_processes'] == 'max':
        settings['n_processes'] = mp.cpu_count()
    else:
        try:
            settings['n_processes'] = int(settings['n_processes'])
        except ValueError:
            parser.print_usage()
            raise Exception('Error: n_processes must be specified as and integer or "max" (current setting: %s)' % settings['n_processes'])

    settings['novel_cut_merge_distance'] = cmd_args.novel_cut_merge_distance
    if 'novel_cut_merge_distance' in settings_file_args:
        settings['novel_cut_merge_distance'] = int(settings_file_args['novel_cut_merge_distance'])
        settings_file_args.pop('novel_cut_merge_distance')

    settings['known_cut_merge_distance'] = cmd_args.known_cut_merge_distance
    if 'known_cut_merge_distance' in settings_file_args:
        settings['known_cut_merge_distance'] = int(settings_file_args['known_cut_merge_distance'])
        settings_file_args.pop('known_cut_merge_distance')

    settings['origin_cut_merge_distance'] = cmd_args.origin_cut_merge_distance
    if 'origin_cut_merge_distance' in settings_file_args:
        settings['origin_cut_merge_distance'] = int(settings_file_args['origin_cut_merge_distance'])
        settings_file_args.pop('origin_cut_merge_distance')

    settings['short_indel_length_cutoff'] = cmd_args.short_indel_length_cutoff
    if 'short_indel_length_cutoff' in settings_file_args:
        settings['short_indel_length_cutoff'] = int(settings_file_args['short_indel_length_cutoff'])
        settings_file_args.pop('short_indel_length_cutoff')

    settings['dedup_input_on_UMI'] = cmd_args.dedup_input_on_UMI
    if 'dedup_input_on_UMI' in settings_file_args:
        settings['dedup_input_on_UMI'] = (settings_file_args['dedup_input_on_UMI'].lower() == 'true')
        settings_file_args.pop('dedup_input_on_UMI')

    settings['suppress_dedup_on_aln_pos_and_UMI_filter'] = cmd_args.suppress_dedup_on_aln_pos_and_UMI_filter
    if 'suppress_dedup_on_aln_pos_and_UMI_filter' in settings_file_args:
        settings['suppress_dedup_on_aln_pos_and_UMI_filter'] = (settings_file_args['suppress_dedup_on_aln_pos_and_UMI_filter'].lower() == 'true')
        settings_file_args.pop('suppress_dedup_on_aln_pos_and_UMI_filter')

    settings['umi_regex'] = cmd_args.umi_regex
    if 'umi_regex' in settings_file_args:
        settings['umi_regex'] = settings_file_args['umi_regex']
        settings_file_args.pop('umi_regex')

    settings['r1_r2_support_max_distance'] = cmd_args.r1_r2_support_max_distance
    if 'r1_r2_support_max_distance' in settings_file_args:
        settings['r1_r2_support_max_distance'] = int(settings_file_args['r1_r2_support_max_distance'])
        settings_file_args.pop('r1_r2_support_max_distance')

    settings['suppress_r2_support_filter'] = cmd_args.suppress_r2_support_filter
    if 'suppress_r2_support_filter' in settings_file_args:
        settings['suppress_r2_support_filter'] = (settings_file_args['suppress_r2_support_filter'].lower() == 'true')
        settings_file_args.pop('suppress_r2_support_filter')

    settings['fastq_r1'] = cmd_args.fastq_r1
    if 'fastq_r1' in settings_file_args:
        settings['fastq_r1'] = settings_file_args['fastq_r1']
        settings_file_args.pop('fastq_r1')
    if not settings['fastq_r1']:
        parser.print_usage()
        raise Exception('Error: fastq_r1 file must be provided (--fastq_r1)')
    if not os.path.isfile(settings['fastq_r1']):
        parser.print_usage()
        raise Exception('Error: fastq_r1 file %s does not exist',settings['fastq_r1'])

    settings['fastq_r2'] = cmd_args.fastq_r2
    if 'fastq_r2' in settings_file_args:
        if 'none' not in settings_file_args['fastq_r2'].lower():
            settings['fastq_r2'] = settings_file_args['fastq_r2']
        settings_file_args.pop('fastq_r2')
    if settings['fastq_r2'] and not os.path.isfile(settings['fastq_r2']):
        raise Exception('Error: fastq_r2 file %s does not exist',settings['fastq_r2'])

    settings['fastq_umi'] = cmd_args.fastq_umi
    if 'fastq_umi' in settings_file_args:
        if 'none' not in settings_file_args['fastq_umi'].lower():
            settings['fastq_umi'] = settings_file_args['fastq_umi']
        settings_file_args.pop('fastq_umi')
    if settings['fastq_umi'] and not os.path.isfile(settings['fastq_umi']):
        raise Exception('Error: fastq_umi file %s does not exist',settings['fastq_umi'])

    settings['genome'] = cmd_args.genome
    if 'genome' in settings_file_args:
        settings['genome'] = settings_file_args['genome']
        settings_file_args.pop('genome')

    if not settings['genome']:
        parser.print_usage()
        raise Exception('Error: the genome reference file must be provided (--genome)')
    if not os.path.isfile(settings['genome']):
        if os.path.isfile(settings['genome']+".fa"):
            settings['genome'] = settings['genome']+".fa"
        else:
            parser.print_usage()
            raise Exception('Error: The genome file %s does not exist'%settings['genome'])
    genome_len_file=settings['genome']+'.fai'
    if not os.path.isfile(genome_len_file):
        raise Exception('Error: The genome length file %s does not exist'%genome_len_file)

    settings['bowtie2_genome'] = cmd_args.bowtie2_genome
    if 'bowtie2_genome' in settings_file_args:
        settings['bowtie2_genome'] = settings_file_args['bowtie2_genome']
        settings_file_args.pop('bowtie2_genome')

    if settings['bowtie2_genome'] is None:
        potential_bowtie2_path = re.sub('.fa$','',settings['genome'])
        if os.path.isfile(potential_bowtie2_path+'.1.bt2'):
            settings['bowtie2_genome']= potential_bowtie2_path
        else:
            raise Exception('Error: bowtie2_genome is required in settings file, pointing to a bowtie2 genome (minus trailing .X.bt2)\nAlternatively, set genome to the .fa file in a bowtie2 directory.')

    if not settings['primer_seq']:
        parser.print_usage()
        raise Exception('Error: the primer sequences used must be provided (--primer_seq)')
    
    for guide in settings['guide_sequences']:
        wrong_nt=CRISPRessoShared.find_wrong_nt(guide)
        if wrong_nt:
            raise Exception('The guide sequence %s contains bad characters: %s'%(guide,'"'+'", "'.join(wrong_nt)+'"'))

    primer_wrong_nt=CRISPRessoShared.find_wrong_nt(settings['primer_seq'])
    if primer_wrong_nt:
        raise Exception('The primer sequence %s contains bad characters: %s'%(settings['primer_seq'],'"'+'", "'.join(wrong_nt)+'"'))

    unused_keys = settings_file_args.keys()
    if len(unused_keys) > 0:
        unused_key_str = '"'+'", "'.join(unused_keys)+'"'
        raise Exception('Unused keys in settings file: ' + unused_key_str)

    # read previous settings (we'll compare to these later to see if we can cache command outputs)
    settings_used_output_file = settings['root']+".settingsUsed.txt"
    previous_settings = {}
    if os.path.exists(settings_used_output_file):
        with open(settings_used_output_file,'r') as fin:
            for line in fin.readlines():
                settings_els = line.strip().split("\t")
                previous_settings[settings_els[0]] = settings_els[1]


    can_use_previous_analysis = False
    if len(previous_settings) > 0:
        can_use_previous_analysis = True
        for setting in settings:
            if setting in ['debug','n_processes']:
                continue
            if setting not in previous_settings:
                can_use_previous_analysis = False
                logger.info(('Not using previous analyses - got new setting %s (%s)')%(setting,settings[setting]))
                break
            elif str(settings[setting]) != str(previous_settings[setting]):
                can_use_previous_analysis = False
                logger.info(('Not using previous analyses - setting for %s has changed (%s (new) vs %s (old))')%(setting,settings[setting],previous_settings[setting]))
                break

    if can_use_previous_analysis:
        logger.info('Repeated settings detected. Using previous analyses if completed.')

    with open (settings_used_output_file,'w') as fout:
        for setting in settings:
            fout.write("%s\t%s\n"%(str(setting),str(settings[setting])))

    settings['can_use_previous_analysis'] = can_use_previous_analysis

    return settings


def assert_dependencies(cutadapt_command='cutadapt',samtools_command='samtools',bowtie2_command='bowtie2',crispresso_command='CRISPResso',casoffinder_command='cas-offinder'):
    """
    Asserts the presence of required software (faidx, bowtie2, casoffinder)

    Args:
        samtools_command: location of samtools to run
        bowtie2_command: location of bowtie2 to run
        crispresso_command: location of crispresso to run
        casoffinder_command: location of casoffinder to run

    Raises exception if any command is not found
    """
    logger = logging.getLogger('CRISPRlungo')
    logger.info('Checking dependencies')

    # check cutadapt
    try:
        cutadapt_result = subprocess.check_output('%s --version'%cutadapt_command, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: cutadapt is required')

    # check faidx
    try:
        faidx_result = subprocess.check_output('%s faidx'%samtools_command, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: samtools faidx is required')
    if not 'Usage: samtools faidx' in str(faidx_result):
        raise Exception('Error: samtools faidx is required')

    #check bowtie2
    try:
        bowtie_result = subprocess.check_output('%s --version'%bowtie2_command, stderr=subprocess.STDOUT,shell=True).decode(sys.stdout.encoding)
        m = re.search(r'bowtie2-align-s version (\d+)\.(\d+)\S+',bowtie_result)
        if m:
            bowtie_major_version = int(m.group(1))
            bowtie_minor_version = int(m.group(2))
            if bowtie_major_version < 2 or bowtie_minor_version < 4:
                    raise Exception('Bowtie version > 2.4 is required (version %s.%s found)'%(bowtie_major_version,bowtie_minor_version))
        else:
            raise Exception('Bowtie2 version cannot be found')
    except Exception:
        raise Exception('Error: bowtie2 is required')

    #check crispresso
    try:
        crispresso_result = subprocess.check_output('%s --version'%crispresso_command, stderr=subprocess.STDOUT,shell=True)
    except Exception as e:
        raise Exception('Error: CRISPResso2 is required ' + str(e))

    #check casoffinder
    try:
        casoffinder_result = subprocess.check_output('%s --version'%casoffinder_command, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: Cas-OFFinder is required')
    if not 'Cas-OFFinder v' in str(casoffinder_result):
        raise Exception('Error: Cas-OFFinder is required')

nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a','n':'n','_':'_','-':'-'})
def reverse(seq):
        return "".join([c for c in seq[-1::-1]])

def reverse_complement(seq):
        return "".join([nt_complement[c] for c in seq[-1::-1]])

def complement(seq):
    return "".join([nt_complement[c] for c in seq[:]])

def get_av_read_len(fastq,number_reads_to_check=50):
    """
    Reads the first few reads of a file to determine read length

    Args:
        fastq: read1 file
        number_reads_to_check: the number of reads to read in

    Returns:
        av_read_len: average read length
    """
    sum_len = 0

    cmd=('z' if fastq.endswith('.gz') else '' ) +('cat < \"%s\"' % fastq)+\
               r''' | awk 'BN {n=0;s=0;} NR%4 == 2 {s+=length($0);n++;} END { printf("%d\n",s/n)}' '''
    av_read_len = float(subprocess.check_output(cmd,shell=True).decode(sys.stdout.encoding).strip())

    return(int(av_read_len))

def get_num_reads_fastq(fastq):
    """
    Counts the number of reads in the specified fastq file

    Args:
        fastq: fastq file

    Returns:
        num_reads: number of reads in the fastq file
    """

    if fastq.endswith('.gz'):
        cmd = 'gunzip -c %s | wc -l'%fastq
    else:
        cmd = 'wc -l %s'%fastq

    res = int(subprocess.check_output(cmd,shell=True).decode(sys.stdout.encoding).split(" ")[0])

    return int(res/4)

def get_cut_sites_casoffinder(root,genome,pam,guides,cleavage_offset,num_mismatches,casoffinder_command,can_use_previous_analysis):
    """
    Gets off-target locations using casoffinder

    Args:
        root: root for written files
        genome: location of genome to use
        pam: PAM sequence of guide
        guides: list of sequences of on-target guides (not including PAM)
        cleavage_offset: position where cleavage occurs (relative to end of spacer seq -- for Cas9 this is -3)
        num_mismatches: number of mismatches to find
        casoffinder_command: location of casoffinder to run
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch
    Returns:
        casoffinder_cut_sites: list of off-target cleavage positions
        casoffinder_cut_annotations: dict of off-target cleavage information
    """

    logger = logging.getLogger('CRISPRlungo')
    logger.info('Calculating cut positions from guide sequence using Cas-OFFinder')
    info_file = root + '.info'

    casoffinder_cut_sites = []
    casoffinder_cut_annotations = {}
    if os.path.isfile(info_file) and can_use_previous_analysis:
        sites_found_count = -1
        with open (info_file,'r') as fin:
            sites_found_line_els = fin.readline().rstrip('\n').split("\t")
            if len(sites_found_line_els) == 2:
                sites_found_count = int(sites_found_line_els[1])
                sites_read_count = 0
                for line in fin:
                    line_els = line.rstrip('\n').split("\t")
                    casoffinder_cut_sites.append(line_els[0])
                    casoffinder_cut_annotations[line_els[0]] = line_els[1].split(",")
                    sites_read_count += 1
                if sites_found_count == sites_read_count:
                    logger.info('Using %d previously-run Cas-OFFinder off-targets'%sites_found_count)
                    return(casoffinder_cut_sites,casoffinder_cut_annotations)
        logger.info('Could not recover previously-run Cas-OFFinder off-targets')


    casoffinder_cut_sites = []
    casoffinder_cut_annotations = {}

    #link in genome -- it throws a seg fault if it's in the bowtie directory
    linked_genome = os.path.abspath(root + '.genome.fa')
    if os.path.exists(linked_genome):
        os.remove(linked_genome)
    os.symlink(os.path.abspath(genome),linked_genome)


    casoffinder_input_file = root + '.input.txt'
    casoffinder_output_file = root + '.output.txt'
    casoffinder_log_file = root + '.log'
    guide_len = max([len(x) for x in guides])
    pam_len = len(pam)
    with open (casoffinder_input_file,'w') as cout:
        cout.write(linked_genome+"\n")
        cout.write("N"*guide_len + pam + "\n")
        for guide in guides:
            cout.write(guide + "N"*pam_len + " " + str(num_mismatches) + "\n")

    casoffinder_cmd = '(%s %s C %s) &>> %s'%(casoffinder_command,casoffinder_input_file,casoffinder_output_file,casoffinder_log_file)

    with open (casoffinder_log_file,'w') as cout:
        cout.write('Linking genome from %s to %s\n'%(genome,linked_genome))
        cout.write('Command used:\n===\n%s\n===\nOutput:\n===\n'%casoffinder_cmd)

    logger.debug(casoffinder_cmd)
    if not os.path.exists(casoffinder_output_file):
        casoffinder_result = subprocess.check_output(casoffinder_cmd,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        logger.debug('Casoffinder output:' + casoffinder_result)
    else:
        logger.debug('Using previously-calculated offtargets')

    os.remove(linked_genome)

    with open (casoffinder_output_file,'r') as cin:
        for line in cin:
            line_els = line.strip().split("\t")
            seq = line_els[0][0:-pam_len]
            chrom = line_els[1]
            loc = line_els[2]
            strand = line_els[4]
            mismatch_count = line_els[5]

            if strand == "+":
                cut_loc = int(loc) + guide_len + 1 + cleavage_offset
                cut_rc_str = "FW"
            else: #neg strand
                cut_loc = int(loc) - cleavage_offset + 1 + pam_len
                cut_rc_str = "RC"
            this_key = chrom+":"+str(cut_loc)
            casoffinder_cut_sites.append(this_key)
            casoffinder_cut_annotations[this_key] = ['Cas-OFFinder OB ' + mismatch_count,'Not-origin',cut_rc_str,seq]
    logger.info('Found %d off-target sites using Cas-OFFinder'%(len(casoffinder_cut_sites)))

    with open (info_file,'w') as fout:
        fout.write('Cas-OFFinder sites\t%d\n'%len(casoffinder_cut_sites))
        for site in casoffinder_cut_sites:
            fout.write(site+"\t"+",".join(casoffinder_cut_annotations[site])+"\n")

    return casoffinder_cut_sites, casoffinder_cut_annotations

def prep_input(root, primer_seq, min_primer_length, guide_seqs, cleavage_offset, fastq_r1, samtools_command, genome, bowtie2_command, bowtie2_genome, can_use_previous_analysis=False, suppress_file_output=False):
    """
    Prepares primer info by identifying genomic location if possible and calculates statistics by analyzing input fastq_r1 file
    Prepares cut info by aligning guides to the genome
    Generally, primer_seq is genomic, and is used to amplify a region near the on-target cut site
    In cases where an exogenous sequence is amplified (e.g. GUIDE-seq) primer_seq may not be genomic
    The origin sequence is determined, and trimmed from input reads

    Args:
        root: root for written files
        primer_seq: sequence of primer sequence.
        min_primer_length: paramter for how many bases must match the primer/origin sequence. This function raises an Exception if the primer (if exogenous) or origin sequence (if genomic) shorter than this cutoff
        guide_seqs: sequences of guides used in experiment
        cleavage_offset: offset for guide where cut occurs
        fastq_r1: input fastq
        samtools_command: location of samtools to run
        genome: location of genome to extract sequence from
        bowtie2_command: location of bowtie2 to run
        bowtie2_genome: bowtie2-indexed genome to align to
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch
        suppress_file_output: suppress writing info file (useful for debugging/testing)

    Returns:
        origin_seq: common amplified region at primer to be removed from input sequences
        cut_sites: list of cut locations for input guides
        cut_annotations: dict of cut_chr:cut_site->annotation for description of cut [type, origin_status, guide_direction] 
            (type=Programmed, Off-target, Known, Casoffinder, etc) (origin_status=Not-origin or Origin:left or Origin:right) (guide_direction='FW' or 'RC' for which direction the guide binds)
        primer_chr: chromosome location of primer
        primer_loc: position of primer
        primer_is_genomic: boolean, true if primer is genomic
        av_read_length: float, average length of reads in fastq_r1
        num_reads_input: int, number of reads in input fastq_r1
    """
    logger = logging.getLogger('CRISPRlungo')
    info_file = root + '.info'

    if os.path.isfile(info_file) and can_use_previous_analysis:
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 8:
                (origin_seq,primer_chr_str,primer_loc_str,primer_seq_str,primer_is_genomic_str,cut_site_count_str,av_read_length_str,num_reads_input_str) = line_els
                primer_chr = None if primer_chr_str == "None" else primer_chr_str
                primer_loc = None if primer_loc_str == "None" else int(primer_loc_str)
                primer_is_genomic = primer_is_genomic_str == 'True'
                cut_site_count = int(cut_site_count_str)
                av_read_length = float(av_read_length_str)
                num_reads_input = int(num_reads_input_str)
                cut_sites = []
                cut_annotations = {}
                read_cut_site_count = 0
                for line in fin:
                    line_els = line.rstrip('\n').split("\t")
                    cut_sites.append(line_els[0])
                    cut_annotations[line_els[0]] = line_els[1].split(",")
                    read_cut_site_count += 1

                if read_cut_site_count == cut_site_count:
                    logger.info('Using previously-prepared primer, cut-site, and sequencing information')
                    return (origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, primer_is_genomic, av_read_length, num_reads_input)
        logger.info('Could not recover previously-created primer, cut-site, and sequencing information. Reanalyzing.')

    logger.info('Preparing primer and cut-site information')

    primer_is_genomic = False # Whether the primer is genomic (e.g. if priming off an exogenous sequence (GUIDE-seq) this would be false) This value is set below after aligning to the genome
    primer_chr = 'exogenous'
    primer_loc = 'None'
    primer_is_rc = 0
    #attempt to align primer to genome
    logger.info('Finding genomic coordinates of primer %s'%(primer_seq))
    primer_aln_cmd = '%s --no-sq --seed 2248 --end-to-end -x %s -c %s' %(bowtie2_command,bowtie2_genome,primer_seq)
    logger.debug(primer_aln_cmd)
    primer_aln_results = subprocess.check_output(primer_aln_cmd,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding).split("\n")
    logger.debug(primer_aln_results)
    if len(primer_aln_results) > 5 and primer_aln_results[3].strip() == '1 (100.00%) aligned exactly 1 time' and len(primer_aln_results[-2].split("\t")) > 5:
        sam_line_els = primer_aln_results[-2].split("\t")
        primer_chr = sam_line_els[2]
        primer_loc = int(sam_line_els[3])
        primer_is_rc = int(sam_line_els[1]) & 0x10
        if primer_is_rc:
            primer_loc = primer_loc + len(primer_seq)
        logger.info('Setting genomic coordinates for primer aligned to %s:%s'%(primer_chr,primer_loc))
        primer_is_genomic = True
    #sometimes alignments are reported by bowtie as ambiguous (e.g. TTGCAATGAAAATAAATGTTT in hg38)
    elif len(primer_aln_results) > 5 and primer_aln_results[4].strip() == '1 (100.00%) aligned >1 times' and len(primer_aln_results[-2].split("\t")) > 5:
        sam_line_els = primer_aln_results[-2].split("\t")
        primer_chr = sam_line_els[2]
        primer_loc = int(sam_line_els[3])
        primer_is_rc = int(sam_line_els[1]) & 0x10
        if primer_is_rc:
            primer_loc = primer_loc + len(primer_seq)
        logger.warning('Primer alignment resulted in an ambiguous alignment')
        logger.info('Setting genomic coordinates for primer aligned to %s:%s'%(primer_chr,primer_loc))
        primer_is_genomic = True
    else:
        logger.info('Primer is not aligned to genome. Assuming exogenous primer sequence')

    cut_sites = []
    cut_annotations = {}
    closest_cut_site_loc = -1
    closest_cut_site_dist = 99**99
    for guide_seq in guide_seqs:
        logger.info('Finding genomic coordinates of guide %s'%(guide_seq))
        guide_aln_cmd = '%s --no-sq --seed 2248 --end-to-end -x %s -c %s' %(bowtie2_command,bowtie2_genome,guide_seq)
        logger.debug(guide_aln_cmd)
        guide_aln_results = subprocess.check_output(guide_aln_cmd,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding).split("\n")
        logger.debug(guide_aln_results)
        if len(guide_aln_results) > 5 and guide_aln_results[3].strip() == '1 (100.00%) aligned exactly 1 time' and len(guide_aln_results[-2].split("\t")) > 5:
            sam_line_els = guide_aln_results[-2].split("\t")
            guide_chr = sam_line_els[2]
            guide_loc = int(sam_line_els[3])
            guide_is_rc = int(sam_line_els[1]) & 0x10
            if not guide_is_rc:
                cut_loc = guide_loc + len(guide_seq) + cleavage_offset
                guide_is_rc_str = "FW"
            else: #neg strand
                cut_loc = guide_loc - cleavage_offset
                guide_is_rc_str = "RC"
            logger.info('Setting cut site for guide %s to %s:%s'%(guide_seq,guide_chr,cut_loc))
            key = guide_chr + ":" + str(cut_loc)
            cut_sites.append(key)
            cut_annotations[key] = ['Programmed','Not-origin',guide_is_rc_str,guide_seq]
            if primer_is_genomic and primer_chr == guide_chr:
                if closest_cut_site_loc == -1:
                    closest_cut_site_dist = abs(primer_loc - cut_loc)
                    closest_cut_site_loc = cut_loc
                elif abs(primer_loc - cut_loc) < closest_cut_site_dist:
                    closest_cut_site_dist = abs(primer_loc - cut_loc)
                    closest_cut_site_loc = cut_loc
        elif len(guide_aln_results) > 5 and guide_aln_results[4].strip() == '1 (100.00%) aligned >1 times' and len(guide_aln_results[-2].split("\t")) > 5:
            sam_line_els = guide_aln_results[-2].split("\t")
            guide_chr = sam_line_els[2]
            guide_loc = int(sam_line_els[3])
            guide_is_rc = int(sam_line_els[1]) & 0x10
            if not guide_is_rc:
                cut_loc = guide_loc + len(guide_seq) + cleavage_offset
                guide_is_rc_str = "FW"
            else: #neg strand
                cut_loc = guide_loc - cleavage_offset
                guide_is_rc_str = "RC"
            logger.warning('Guide alignment resulted in an ambiguous alignment')
            logger.info('Setting cut site for guide %s to %s:%s'%(guide_seq,guide_chr,cut_loc))
            key = guide_chr + ":" + str(cut_loc)
            cut_sites.append(key)
            cut_annotations[key] = ['Programmed','Not-origin',guide_is_rc_str,guide_seq]
            if primer_is_genomic and primer_chr == guide_chr:
                if closest_cut_site_loc == -1:
                    closest_cut_site_dist = abs(primer_loc - cut_loc)
                    closest_cut_site_loc = cut_loc
                elif abs(primer_loc - cut_loc) < closest_cut_site_dist:
                    closest_cut_site_dist = abs(primer_loc - cut_loc)
                    closest_cut_site_loc = cut_loc
        else:
            logger.warning('Bowtie alignment results:')
            logger.warning(guide_aln_results)
            raise Exception('Could not find unique genomic coordinates for guide %s'%guide_seq)

    origin_seq = primer_seq
    if primer_is_genomic and closest_cut_site_loc == -1:
        primer_is_genomic = False
        logger.warning('Primer was found to be gnomic (%s:%s) but no cut sites were found on chromosome %s. Assuming primer is the origin sequence'%(primer_chr,primer_loc,primer_chr))
    elif primer_is_genomic:
        logger.debug('Closest cut site to primer (%s:%s) is %s:%s'%(primer_chr, primer_loc, primer_chr, closest_cut_site_loc))
        origin_start = primer_loc
        origin_end = closest_cut_site_loc
        cut_annotations[primer_chr + ":" + str(closest_cut_site_loc)][1] = 'Origin:left' #add annotation noting that this cut is the origin, and the primer is to the left
        if primer_is_rc:
            origin_start = closest_cut_site_loc
            origin_end = primer_loc
            cut_annotations[primer_chr + ":" + str(closest_cut_site_loc)][1] = 'Origin:right'
        if origin_start > origin_end:
            raise Exception('Primer does not point toward the closest cut site. Please check the orientation of the primer sequence parameter')

        origin_seq = subprocess.check_output(
                '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,primer_chr,origin_start,origin_end-1),shell=True).decode(sys.stdout.encoding).strip()
        if origin_seq == '':
            raise Exception('Cannot retrieve sequence information from fasta file ' + genome)
        if primer_is_rc:
            origin_seq = reverse_complement(origin_seq)
        logger.debug('Got origin sequence: ' + origin_seq)

    if len(origin_seq) < min_primer_length:
        raise Exception('The min_primer_length parameter must be less than the length of the primer/origin sequence (Origin: %s (length %d), min_primer_length: %d)'%(origin_seq,len(origin_seq),min_primer_length))

    logger.debug('Getting read length and number of total reads')
    av_read_length = get_av_read_len(fastq_r1)
    num_reads_input = get_num_reads_fastq(fastq_r1)

    cut_site_count = len(cut_sites)
    if not suppress_file_output:
        with open(info_file,'w') as fout:
            fout.write("\t".join(['origin_seq','primer_chr','primer_loc','primer_seq','primer_is_genomic','cut_site_count','av_read_length','num_reads_input'])+"\n")
            fout.write("\t".join([str(x) for x in [origin_seq,primer_chr,primer_loc,primer_seq,primer_is_genomic,cut_site_count,av_read_length,num_reads_input]])+"\n")
            for cut_site in cut_sites:
                fout.write(cut_site + "\t" + ",".join(cut_annotations[cut_site])+"\n")

    logger.info('Finished preparing primer, cut-site, and sequencing information')

    return (origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, primer_is_genomic, av_read_length, num_reads_input)

def dedup_input_file(root,fastq_r1,fastq_r2,umi_regex,min_umi_seen_to_keep_read=0,write_UMI_counts=False,can_use_previous_analysis=False):
    """
    Deduplicates fastq files based on UMI (UMI is assumed to be the last part of the read ID after the ':'

    Args:
        root: root for written files
        fastq_r1: R1 reads to dedup
        fastq_r2: R2 reads to dedup
        umi_regex: string specifying regex that UMI must match
        min_umi_seen_to_keep_read: min number of times a umi must be seen to keep the umi and read (e.g. if set to 2, a read-UMI pair that is only seen once will be discarded)
        write_UMI_counts: if True, writes a file with the UMI counts
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch
    Returns:
        fastq_r1_dedup: fastq_r1 file of deduplicated reads
        fastq_r2_dedup: fastq_r2 file of deduplicated reads
        tot_read_count: number of total reads read
        count_with_regex: number of reads with umi matching regex
        post_dedup_count: number of reads post deduplication
        post_dedup_read_count: number of reads that contributed to the deduplicated reads (post_dedup_read_count reads were seen that were deduplicated to post_dedup_count reads. If post_dedup_count == post_dedup_read_count then every UMI-read pair was seen once.)
    """

    logger = logging.getLogger('CRISPRlungo')

    dedup_stats_file = root + ".info"
    fastq_r1_dedup = None
    fastq_r2_dedup = None
    tot_read_count = -1
    count_with_regex = -1
    post_dedup_count = -1
    post_dedup_read_count = -1
    #if already finished, attempt to read in stats
    if os.path.isfile(dedup_stats_file) and can_use_previous_analysis:
        with open(dedup_stats_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) > 3:
                (fastq_r1_dedup_str,fastq_r2_dedup_str,tot_read_count_str,count_with_regex_str,post_dedup_count_str,post_dedup_read_count_str) = line_els
                fastq_r1_dedup = None if fastq_r1_dedup_str == "None" else fastq_r1_dedup_str
                fastq_r2_dedup = None if fastq_r2_dedup_str == "None" else fastq_r2_dedup_str
                tot_read_count = int(tot_read_count_str)
                count_with_regex = int(count_with_regex_str)
                post_dedup_count = int(post_dedup_count_str)
                post_dedup_read_count = int(post_dedup_read_count_str)
    if tot_read_count != -1:
        logger.info('Using previously-deduplicated fastqs with %d reads'%post_dedup_count)
        return(fastq_r1_dedup,fastq_r2_dedup, tot_read_count, count_with_regex, post_dedup_count, post_dedup_read_count)

    #otherwise perform deduplication (check if we were able to read in stats as well -- if we couldn't read them in, tot_read_count will be -1
    logger.info('Deduplicating input fastqs')

    #first, create the regex that matches UMIs
    umi_regex_string = umi_regex

    umi_regex_string = umi_regex_string.replace('I','([ATCG])')
    umi_regex_string = umi_regex_string.replace('N','([ATCG])')
    umi_regex_string = umi_regex_string.replace('R','([AG])')
    umi_regex_string = umi_regex_string.replace('Y','([CT])')
    umi_regex_string = umi_regex_string.replace('S','([GC])')
    umi_regex_string = umi_regex_string.replace('W','([AT])')
    umi_regex_string = umi_regex_string.replace('K','([GT])')
    umi_regex_string = umi_regex_string.replace('M','([AC])')
    umi_regex_string = umi_regex_string.replace('B','([CGT])')
    umi_regex_string = umi_regex_string.replace('D','([AGT])')
    umi_regex_string = umi_regex_string.replace('H','([ACT])')
    umi_regex_string = umi_regex_string.replace('V','([ACG])')

    umi_regex = re.compile(umi_regex_string)

    tot_read_count = 0
    count_with_regex = 0

    umi_keys_with_most_counts = {} #umi->key (key has most counts)
    umi_key_counts = {} #umi->count (count of fastqs key has seen)

    umi_seq_counts = {} #umi_seq->count of that pair
    umi_seq_best_qual_fastqs = {} #umi_seq->fastq to be printed
    umi_seq_best_qual_sum = {} #best qual sum for the best fastq

    if fastq_r1.endswith('.gz'):
        f1_in = gzip.open(fastq_r1,'rt')
    else:
        f1_in = open(fastq_r1,'rt')

    if fastq_r2.endswith('.gz'):
        f2_in = gzip.open(fastq_r2,'rt')
    else:
        f2_in = open(fastq_r2,'rt')

    #now iterate through f1/f2/umi files
    while (1):
        f1_id_line   = f1_in.readline().strip()
        f1_seq_line  = f1_in.readline().strip()
        f1_plus_line = f1_in.readline()
        f1_qual_line = f1_in.readline().strip()

        if not f1_qual_line : break
        if not f1_plus_line.startswith("+"):
            raise Exception("Fastq %s cannot be parsed (%s%s%s%s) "%(fastq_r1,f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))
        tot_read_count += 1

        f2_id_line   = f2_in.readline().strip()
        f2_seq_line  = f2_in.readline().strip()
        f2_plus_line = f2_in.readline()
        f2_qual_line = f2_in.readline().strip()

        this_UMI = f1_id_line.split(":")[-1].upper()

        if umi_regex.match(this_UMI):
            count_with_regex += 1

            #group 1 is the whole string
            #this_key = this_UMI + " # " + f1_seq_line + f2_seq_line
            this_key = this_UMI


            qual_sum = np.sum(np.fromstring(f1_qual_line + f2_qual_line,dtype=np.uint8))
            if this_key not in umi_seq_counts:
                umi_seq_counts[this_key] = 1
                umi_seq_best_qual_sum[this_key] = qual_sum
                umi_seq_best_qual_fastqs[this_key] = (
                        f1_id_line + "\n" + f1_seq_line + "\n" + f1_plus_line + f1_qual_line,
                        f2_id_line + "\n" + f2_seq_line + "\n" + f2_plus_line + f2_qual_line
                        )
            else:
                umi_seq_counts[this_key] += 1
                #if this sequence has the highest quality, store it
                if umi_seq_best_qual_sum[this_key] < qual_sum:
                    umi_seq_best_qual_sum[this_key] = qual_sum
                    umi_seq_best_qual_fastqs[this_key] = (
                            f1_id_line + "\n" + f1_seq_line + "\n" + f1_plus_line + f1_qual_line,
                            f2_id_line + "\n" + f2_seq_line + "\n" + f2_plus_line + f2_qual_line
                            )

            if this_UMI not in umi_key_counts:
                umi_key_counts[this_UMI] = 1
                umi_keys_with_most_counts[this_UMI] = this_key
            else:
                umi_key_counts[this_UMI] += 1
                #if this sequence is the most seen for this UMI, store it
                if umi_seq_counts[this_key] > umi_key_counts[this_UMI]:
                    umi_keys_with_most_counts[this_UMI] = this_key
    #finished iterating through fastq file
    f1_in.close()
    f2_in.close()


    if tot_read_count == 0:
        raise Exception("UMI dedup failed. Got no reads from " + fastq_r1 + " and " + fastq_r2 )

    umi_list = sorted(umi_key_counts, key=lambda k: umi_key_counts[k])
    umi_count = len(umi_list)

    logger.info('Read %d reads'%tot_read_count)
    logger.info("Processed " + str(umi_count) + " UMIs")
    if umi_count == 1 and umi_list[0] == '':
        raise Exception("Error: only the empty barcode '' was found.")

    fastq_r1_dedup = root + '.r1.gz'
    f1_out = gzip.open(fastq_r1_dedup, 'wt')
    fastq_r2_dedup = root + '.r2.gz'
    f2_out = gzip.open(fastq_r2_dedup, 'wt')

    collision_count = 0
    collision_count_reads = 0
    too_few_reads_count = 0
    too_few_reads_count_reads = 0
    post_dedup_count = 0
    post_dedup_read_count = 0
    collided_umi_reads = []
    for umi_seq in umi_seq_counts:
        if umi_seq_counts[umi_seq] < min_umi_seen_to_keep_read:
            too_few_reads_count += 1
            too_few_reads_count_reads += umi_seq_counts[umi_seq]
        else:
            (seq1,seq2) = umi_seq_best_qual_fastqs[umi_seq]
            f1_out.write(seq1+"\n")
            f2_out.write(seq2+"\n")
            post_dedup_count += 1
            post_dedup_read_count += umi_seq_counts[umi_seq]

    f1_out.close()
    f2_out.close()

    if write_UMI_counts:
        fout = open(root+".umiCounts.txt","w")
        fout.write('UMI\tCount\n')
        for umi in sorted(umi_key_counts):
            fout.write(umi + "\t" + str(umi_key_counts[umi]) + "\n")
        logger.info('Wrote UMI counts to ' + root+".umiCounts.txt")


    logger.info('Wrote %d deduplicated reads'%post_dedup_count)
    with open(dedup_stats_file,'w') as fout:
        fout.write("\t".join(["fastq_r1_dedup","fastq_r2_dedup","tot_read_count","count_with_regex","post_dedup_count","post_dedup_read_count"])+"\n")
        fout.write("\t".join([str(x) for x in [fastq_r1_dedup,fastq_r2_dedup,tot_read_count,count_with_regex,post_dedup_count,post_dedup_read_count]])+"\n")

    return(fastq_r1_dedup,fastq_r2_dedup, tot_read_count, count_with_regex, post_dedup_count, post_dedup_read_count)

def add_umi_from_umi_file(root,fastq_r1,fastq_r2,fastq_umi,can_use_previous_analysis=False):
    """
    Adds the UMI to the read ID from a UMI file (a third file with the UMI sequences per read)

    Args:
        root: root for written files
        fastq_r1: R1 reads to dedup
        fastq_r2: R2 reads to dedup
        fastq_umi: UMI fastq to dedup on
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch
    Returns:
        fastq_r1_umi: fastq_r1 with umi added
        fastq_r2_umi: fastq_r2 with umi added
        tot_read_count: number of total reads read
    """

    logger = logging.getLogger('CRISPRlungo')
    umi_stats_file = root + ".log"

    fastq_r1_dedup = "NA"
    fastq_r2_dedup = "NA"
    tot_read_count = -1
    #if already finished, attempt to read in stats
    if os.path.isfile(umi_stats_file) and can_use_previous_analysis:
        with open(umi_stats_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 3:
                logger.info('Using previously-generated fastqs with UMIs')
                (fastq_r1_dedup,fastq_r2_dedup,tot_read_count_str) = line_els
                tot_read_count = int(tot_read_count_str)

    #otherwise perform umi adding (check if we were able to read in stats as well -- if we couldn't read them in, tot_read_count will be -1
    if tot_read_count == -1:
        logger.info('Adding UMIs to input fastqs')

        tot_read_count = 0

        if fastq_r1.endswith('.gz'):
            f1_in = gzip.open(fastq_r1,'rt')
        else:
            f1_in = open(fastq_r1,'r')

        if fastq_r2.endswith('.gz'):
            f2_in = gzip.open(fastq_r2,'rt')
        else:
            f2_in = open(fastq_r2,'r')

        if fastq_umi.endswith('.gz'):
            umi_in = gzip.open(fastq_umi,'rt')
        else:
            umi_in = open(fastq_umi,'r')

        fastq_r1_dedup = root + '.r1.gz'
        f1_out = gzip.open(fastq_r1_dedup, 'wt')
        fastq_r2_dedup = root + '.r2.gz'
        f2_out = gzip.open(fastq_r2_dedup, 'wt')

        #now iterate through f1/f2/umi files
        while (1):
            f1_id_line   = f1_in.readline().strip()
            f1_seq_line  = f1_in.readline()
            f1_plus_line = f1_in.readline()
            f1_qual_line = f1_in.readline()

            if not f1_qual_line : break
            if not f1_plus_line.startswith("+"):
                raise Exception("Fastq %s cannot be parsed (%s%s%s%s) "%(fastq_r1,f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))
            tot_read_count += 1

            f2_id_line   = f2_in.readline().strip()
            f2_seq_line  = f2_in.readline()
            f2_plus_line = f2_in.readline()
            f2_qual_line = f2_in.readline()

            umi_id_line   = umi_in.readline()
            umi_seq_line  = umi_in.readline().strip()
            umi_plus_line = umi_in.readline()
            umi_qual_line = umi_in.readline()

            this_UMI = umi_seq_line

            f1_out.write(f1_id_line + ":" + this_UMI + "\n" + f1_seq_line + f1_plus_line + f1_qual_line)
            f2_out.write(f2_id_line + ":" + this_UMI + "\n" + f2_seq_line + f2_plus_line + f2_qual_line)


        #finished iterating through fastq file
        f1_in.close()
        f2_in.close()
        umi_in.close()
        f1_out.close()
        f2_out.close()


        if tot_read_count == 0:
            raise Exception("UMI command failed. Got no reads from " + fastq_r1 + " and " + fastq_r2 )

        logger.info('Added UMIs fo %d reads'%tot_read_count)
        with open(umi_stats_file,'w') as fout:
            fout.write("\t".join(["fastq_r1_dedup","fastq_r2_dedup","tot_read_count"])+"\n")
            fout.write("\t".join([str(x) for x in [fastq_r1_dedup,fastq_r2_dedup,tot_read_count]])+"\n")

    #done processing just plotting now
    return(fastq_r1_dedup,fastq_r2_dedup)

def filter_on_primer(root,fastq_r1,fastq_r2,origin_seq,min_primer_aln_score,min_primer_length=10,min_read_length=30,transposase_adapter_seq='CTGTCTCTTATACACATCTGACGCTGCCGACGA',n_processes=1,cutadapt_command='cutadapt',keep_intermediate=False,suppress_plots=False,can_use_previous_analysis=False):
    """
    Trims the primer from the input reads, only keeps reads with the primer present in R1
    Also trims the transposase adapter sequence from reads and filters reads that are too short

    Args:
        root: root for written files
        fastq_r1: R1 reads to trim
        fastq_r2: R2 reads to trim (or None)
        origin_seq: primer sequence to trim, if genomic, this includes the entire genomic sequence from the primer to the first cut site
        min_primer_aln_score: minimum score for alignment between primer/origin sequence and read sequence
        min_primer_length: minimum length of sequence that matches between the primer/origin sequence and the read sequence
        min_read_len: minimum r1 read length to keep (shorter tagmented reads are filtered)
        transposase_adapter_seq: transposase adapter seq to trim from R1
        n_processes: Number of processes to run on
        cutadapt_command: command to run cutadapt
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        suppress_plots: if true, plotting will be suppressed
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    Returns:
        filtered_on_primer_fastq_r1: fastq_r1 containing reads with primer
        filtered_on_primer_fastq_r2: fastq_r2 containing reasd with primer in r1 (or None)
        post_filter_on_primer_read_count: number of reads kept with primer
        filter_on_primer_plot_obj: plot object showing number of reads filtered
    """

    logger = logging.getLogger('CRISPRlungo')

    trim_primer_stats_file = root + ".info"
    filtered_on_primer_fastq_r1 = None
    filtered_on_primer_fastq_r2 = None
    post_trim_read_count = -1
    #if already finished, attempt to read in stats
    if os.path.isfile(trim_primer_stats_file) and can_use_previous_analysis:
        with open(trim_primer_stats_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 3:
                (filtered_on_primer_fastq_r1_str,filtered_on_primer_fastq_r2_str,post_trim_read_count_str) = line_els
                filtered_on_primer_fastq_r1 = None if filtered_on_primer_fastq_r1_str == "None" else filtered_on_primer_fastq_r1_str
                filtered_on_primer_fastq_r2 = None if filtered_on_primer_fastq_r2_str == "None" else filtered_on_primer_fastq_r2_str
                post_trim_read_count = int(post_trim_read_count_str)

                filter_on_primer_plot_obj_str = fin.readline().rstrip('\n')
                filter_on_primer_plot_obj = None
                if filter_on_primer_plot_obj_str != "" and filter_on_primer_plot_obj_str != "None":
                    filter_on_primer_plot_obj = PlotObject.from_json(filter_on_primer_plot_obj_str)

                if os.path.exists(filtered_on_primer_fastq_r1):
                    logger.info('Using %d previously-filtered fastq sequences trimmed for primers'%post_trim_read_count)

                    return(filtered_on_primer_fastq_r1,filtered_on_primer_fastq_r2,post_trim_read_count,filter_on_primer_plot_obj)
                else:
                    logger.debug('File %s does not exist.'%filtered_on_primer_fastq_r1)


        logger.info('Could not recover previously-filtered results. Reprocessing.')

    logger.info('Filtering and trimming primers from reads')

    post_trim_read_count = 'NA'
    too_short_read_count = 'NA'
    untrimmed_read_count = 'NA'
    #prep for alignments using ssw

    sLibPath = os.path.dirname(os.path.abspath(__file__))+"/lib"
    ssw_primer = ssw_lib.CSsw(sLibPath)
    ssw_align_primer = SSWPrimerAlign(ssw_primer,origin_seq)

    if fastq_r2 is None:
        filtered_on_primer_fastq_r1 = root + ".has_primer.fq.gz"
        filtered_on_primer_fastq_r2 = None
        post_trim_read_count, too_short_read_count, untrimmed_read_count = trimPrimersSingle(fastq_r1,filtered_on_primer_fastq_r1,min_primer_aln_score,min_primer_length,ssw_align_primer)
    else:
        filtered_on_primer_fastq_r1 = root + ".r1.has_primer.fq.gz"
        filtered_on_primer_fastq_r2 = root + ".r2.has_primer.fq.gz"
        rc_origin_seq = reverse_complement(origin_seq)
        ssw_primer_rc = ssw_lib.CSsw(sLibPath)
        ssw_align_primer_rc = SSWPrimerAlign(ssw_primer_rc,rc_origin_seq)
        post_trim_read_count, too_short_read_count, untrimmed_read_count = trimPrimersPair(fastq_r1,fastq_r2,filtered_on_primer_fastq_r1,filtered_on_primer_fastq_r2,min_primer_aln_score,min_primer_length,ssw_align_primer,ssw_align_primer_rc)

    contain_adapter_count = 0
    filtered_too_short_adapter_count = 0

    if post_trim_read_count > 0 and transposase_adapter_seq is not None:
        logger.info('Trimming transposase sequence from reads')
        cutadapt_transposase_log = root + ".cutadapt.transposase_adapter.log"
        if fastq_r2 is None:
            filtered_for_adapters_fastq_r1 = root + ".trimmed.fq.gz"
            filtered_for_adapters_fastq_r2 = None
            trim_command = "%s -a %s -o %s --minimum-length %d --cores %s %s > %s"%(cutadapt_command,transposase_adapter_seq,filtered_for_adapters_fastq_r1,min_read_length,n_processes,filtered_on_primer_fastq_r1,cutadapt_transposase_log)
        else:
            filtered_for_adapters_fastq_r1 = root + ".r1.trimmed.fq.gz"
            filtered_for_adapters_fastq_r2 = root + ".r2.trimmed.fq.gz"
            trim_command = "%s -a %s -o %s -p %s --minimum-length %d --pair-filter=first --cores %s %s %s > %s"%(cutadapt_command,transposase_adapter_seq,filtered_for_adapters_fastq_r1, filtered_for_adapters_fastq_r2,min_read_length,n_processes,filtered_on_primer_fastq_r1,filtered_on_primer_fastq_r2,cutadapt_transposase_log)

        logger.debug(trim_command)
        trim_result = subprocess.check_output(trim_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)

        if not os.path.exists(cutadapt_transposase_log):
            logger.error('Error found while running command:\n'+trim_command+"\nOutput: "+trim_result)

        adapter_too_short_read_count = 'NA'
        post_trim_read_count = 'NA' #overwrites post_trim_read_count above - ok because this is the final number if adapter seqs are given
        with open(cutadapt_transposase_log,'r') as fin:
            for line in fin:
                m = re.search(r'written \(passing filters\):\s+([\d,]+) \(\d.*%\)',line)
                if m:
                    post_trim_read_count = int(m.group(1).replace(',',''))

                m = re.search(r'that were too short:\s+([\d,]+) \(\d.*%\)',line)
                if m:
                    adapter_too_short_read_count = int(m.group(1).replace(',',''))

        if post_trim_read_count == 'NA' or adapter_too_short_read_count == 'NA':
            logger.error('Could not parse trim read count from file ' + cutadapt_transposase_log)
        too_short_read_count += adapter_too_short_read_count

        if not keep_intermediate:
            logger.debug('Deleting intermediate fastq files trimmed for primers')
            os.remove(filtered_on_primer_fastq_r1)
            if fastq_r2 is not None:
                os.remove(filtered_on_primer_fastq_r2)

        filtered_on_primer_fastq_r1 = filtered_for_adapters_fastq_r1
        filtered_on_primer_fastq_r2 = filtered_for_adapters_fastq_r2

    #plot summary
    filter_labels = ['Filtered: too short','Filtered: no primer','Retained']
    values = [too_short_read_count,untrimmed_read_count,post_trim_read_count]
    filter_on_primer_plot_obj_root = root + ".counts"
    with open(filter_on_primer_plot_obj_root+".txt",'w') as summary:
        summary.write("\t".join(filter_labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
    if suppress_plots:
        filter_on_primer_plot_obj = None
    else:
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(filter_labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(filter_labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Read filtering based on primer/origin presence')
        plt.savefig(filter_on_primer_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(filter_on_primer_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(filter_labels,values)])
        plot_label = "Reads were filtered to contain primer origin sequence <p class='text-break text-monospace'>" + origin_seq +"</p>"
        plot_label += "Reads must contain at least " + str(min_primer_length) + " bases that match the primer/origin sequence<br>"
        plot_label += "Reads containing the primer/origin sequence but shorter than " + str(min_read_length) + " were filtered as too short."
        filter_on_primer_plot_obj = PlotObject(
                plot_name = filter_on_primer_plot_obj_root,
                plot_title = 'Read filtering based on primer/origin presence',
                plot_label = plot_label + '<br>'+plot_count_str,
                plot_datas = [
                    ('Primer/origin information',trim_primer_stats_file),
                    ('Primer filtering statistics',filter_on_primer_plot_obj_root + ".txt")
                    ]
                )


    with open(trim_primer_stats_file,'w') as fout:
        fout.write("\t".join(["filtered_on_primer_fastq_r1","filtered_on_primer_fastq_r2","post_trim_read_count"])+"\n")
        fout.write("\t".join([str(x) for x in [filtered_on_primer_fastq_r1,filtered_on_primer_fastq_r2,post_trim_read_count]])+"\n")

        filter_on_primer_plot_obj_str = "None"
        if filter_on_primer_plot_obj is not None:
            filter_on_primer_plot_obj_str = filter_on_primer_plot_obj.to_json()
        fout.write(filter_on_primer_plot_obj_str+"\n")


    logger.info('Keeping %d reads with primer sequence'%post_trim_read_count)


    return(filtered_on_primer_fastq_r1,filtered_on_primer_fastq_r2,post_trim_read_count,filter_on_primer_plot_obj)

def align_reads(root,fastq_r1,fastq_r2,bowtie2_reference,bowtie2_command='bowtie2',bowtie2_threads=1,use_old_bowtie=False,samtools_command='samtools',keep_intermediate=False,can_use_previous_analysis=False):
    """
    Aligns reads to the provided reference

    Args:
        root: root for written files
        fastq_r1: fastq_r1 to align
        fastq_r2: fastq_r2 to align
        bowtie2_reference: bowtie2 reference to align to (either artificial targets or reference)
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        use_old_bowtie: boolean for whether to use old bowtie2 version. Versions before v2.3.4.2 didn't implement the '--soft-clipped-unmapped-tlen' option
        samtools_command: location of samtools to run
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    Returns:
        mapped_bam_file: aligned reads in bam format
    """

    logger = logging.getLogger('CRISPRlungo')

    mapped_bam_file = root + ".bam"

    info_file = root + '.info'
    if os.path.isfile(info_file) and os.path.isfile(mapped_bam_file) and can_use_previous_analysis:
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 1:
                mapped_bam_file = line_els[0]
                logger.info('Using previously-performed alignment')
                return mapped_bam_file
        logger.info('Could not recover previously-performed alignment. Reanalyzing.')

    tlen_option = "--soft-clipped-unmapped-tlen"
    if use_old_bowtie:
        tlen_option = ""


    bowtie_log = root + '.bowtie2Log'
    if fastq_r2 is not None: #paired-end reads
        logger.info('Aligning paired reads using %s'%(bowtie2_command))
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --very-sensitive-local {tlen_option} --threads {bowtie2_threads} -x {bowtie2_reference} -1 {fastq_r1} -2 {fastq_r2} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                tlen_option=tlen_option,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logger.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        if 'error' in aln_result.lower():
            logger.error('Error found while running command:\n'+aln_command+"\nOutput: "+aln_result)
            raise Exception('Alignment error: ' + aln_result)

        logger.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning paired reads\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(aln_command,aln_result))
    #unpaired reads
    else:
        logger.info('Aligning single-end reads using %s'%(bowtie2_command))
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --very-sensitive-local {tlen_option} --threads {bowtie2_threads} -x {bowtie2_reference} -U {fastq_r1} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                tlen_option=tlen_option,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logger.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        if 'error' in aln_result.lower():
            logger.error('Error found while running command:\n'+aln_command+"\nOutput: "+aln_result)
            raise Exception('Alignment error: ' + aln_result)

        logger.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning single-end reads\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(aln_command,aln_result))



    with open(info_file,'w') as fout:
        fout.write("\t".join([str(x) for x in ["mapped_bam_file"]])+"\n")
        fout.write("\t".join([str(x) for x in [mapped_bam_file]])+"\n")
    logger.info('Finished genome alignment')
    return(mapped_bam_file)

def make_final_read_assignments(root,genome_mapped_bam,origin_seq,
    cut_sites,cut_annotations,cut_classification_annotations,guide_seqs,cleavage_offset,min_primer_length,genome,r1_r2_support_max_distance=100000,
    novel_cut_merge_distance=50,known_cut_merge_distance=50,origin_cut_merge_distance=10000,short_indel_length_cutoff=50,guide_homology_max_gaps=2,guide_homology_max_mismatches=5,
    arm_min_matched_start_bases=10,arm_max_clipped_bases=0,genome_map_resolution=1000000,crispresso_max_indel_size=50,suppress_dedup_on_aln_pos_and_UMI_filter=False,
    suppress_r2_support_filter=False,suppress_poor_alignment_filter=False,write_discarded_read_info=False,
    samtools_command='samtools',keep_intermediate=False,suppress_plots=False,can_use_previous_analysis=False):
    """
    Makes final read assignments (after deduplicating based on UMI and alignment location)

    Args:
        root: root for written files
        genome_mapped_bam: bam of reads aligned to genome
        origin_seq: common amplified region at primer (to first cut site)
        cut_sites: array of known/input cut sites
        cut_annotations: dict of cut_chr:cut_site->annotation for description of cut [type, origin_status, guide_direction] 
            (type=Programmed, Off-target, Known, Casoffinder, etc) (origin_status=Not-origin or Origin:left or Origin:right) (guide_direction='FW' or 'RC' for which direction the guide binds)
        cut_classification_annotations: list of cut_chr:cut_site:direction:annotation (as provided by users as parameters)
        guide_seqs: sequences of guides used in experiment (required here for finding homology at novel cut sites)
        cleavage_offset: position where cleavage occurs (relative to end of spacer seq -- for Cas9 this is -3)
        min_primer_length: minimum length of sequence that matches between the primer/origin sequence and the read sequence
        genome: path to genome fa file (required here for finding genome sequences for guide homology)
        r1_r2_support_max_distance: max distance between r1 and r2 for the read pair to be classified as 'supported' by r2
        novel_cut_merge_distance: cuts within this distance (bp) will be merged
        known_cut_merge_distance: novel cuts within this distance to a site with homology to a given guide. Novel cuts are first merged to the nearest site with homology (within known_cut_merge_distance), and remaining cuts are merged with eachother (specified by the 'novel_cut_merge_distance' parameter).
        origin_cut_merge_distance: reads aligned with this distance to the origin will be merged to the origin
        short_indel_length_cutoff: indels this length or shorter will be classified 'short indels' while those longer will be classified 'long indels'
        guide_homology_max_gaps: when searching a sequence for homology to a guide sequence, it will be homologous if his has this number of gaps or fewer in the alignment of the guide to a sequence
        guide_homology_max_mismatches: when searching a sequence for homology to a guide sequence, it will be homologous if his has this number of mismatches or fewer in the alignment of the guide to a sequence
        arm_min_matched_start_bases: number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each 'side' a read.
        arm_max_clipped_bases: the maximum number of bases that can be clipped on both sides for a read to be accepted.
        genome_map_resolution: window size (bp) for reporting number of reads aligned
        crispresso_max_indel_size: maximum length of indel (as determined by genome alignment position) for a read to be analyzed by CRISPResso. Reads with indels longer than this length will not be analyzed by CRISPResso.
        dedup_based_on_aln_pos: if true, reads with the same UMI and alignment positions will be removed from analysis
        suppress_r2_support_filter: if true, reads without r2 support will be included in final analysis and counts. If false, reads without r2 support will be filtered from the final analysis.
        suppress_poor_alignment_filter: if true, reads with poor alignment (fewer than --arm_min_matched_start_bases matches at the alignment ends or more than --arm_max_clipped_bases on both sides of the read) are included in the final analysis and counts. If false, reads with poor alignment are filtered from the final analysis.
        write_discarded_read_info: if true, files are written containing info for reads that are discarded from final analysis
        samtools_command: location of samtools to run
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        suppress_plots: if true, plotting will be suppressed
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    Returns:
        final_assignment_filename: filename of final assignments for each read
        final_read_ids_for_crispresso: dict of readID=>cut
        final_cut_counts: dict of cutID=>count how many times each cut was seen -- when we iterate through the reads (fastq) in the next step, we check to see that the cut was seen above a threshold before printing those reads. The count is stored in this dict.
        cut_classification_lookup: dict of cutID=> type (e.g. Linear, Translocation, etc)
        final_read_count: number of r1 reads used in final analysis
        discarded_read_counts: dict of category -> count for reasons reads were discarded
        classification_read_counts: dict of classification -> count for final assignment categories
        classification_indel_read_counts: dict of classification -> count for final assignment categories (including 'short indels' category)
        chr_aln_plot_obj: plot showing alignment locations of reads
        tlen_plot_object: plot showing distribution of template lengths
        deduplication_plot_obj: plot showing how many reads were deuplicated using UMI + location
        classification_plot_obj: plot showing assignments (linear, translocation, etc)
        classification_indel_plot_obj: plot showing assignments (linear, translocation, etc) with indel classification (short indels)
        tx_order_plot_obj: plot showing ordered translocations and counts
        tx_count_plot_obj: plot showing translocations and counts
        tx_circos_plot_obj: circos plot showing translocations and counts
        origin_indel_hist_plot_obj: plot of indel lenghts for origin
        origin_inversion_hist_plot_obj: plot of inversion lengths at origin (reads pointing opposite directions)
        origin_deletion_hist_plot_obj: plot of deletion lengths for ontarget (reads pointing same direction/away from the primer)
        origin_indel_depth_plot_obj: plot of read depth around origin
        r1_r2_support_plot_obj: plot of how many reads for which r1 and r2 support each other
        r1_r2_support_dist_plot_obj: plot of how far apart the r1/r2 supporting reads were aligned from each other (showing reads supported and not supported by r2 separately)
        discarded_reads_plot_obj: plot of why reads were discarded from final analysis

    """
    logger = logging.getLogger('CRISPRlungo')

    final_file = root + '.final_assignments.txt'
    info_file = root + '.info'
    #check to see if this analysis has been completed previously
    if os.path.isfile(info_file) and os.path.isfile(final_file) and can_use_previous_analysis:
        read_total_read_count = 0 #how many lines we read in this time by reading the final assignments file
        previous_total_read_count = -1 #how many lines were reported to be read the previous time (if the analysis was completed previously)
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 9:
                (total_reads_processed,total_r1_processed,discarded_reads_unaligned,discarded_reads_duplicates,discarded_reads_not_supported_by_R2,discarded_reads_poor_alignment,final_total_count,found_cut_point_total,final_cut_point_total) = [int(x) for x in line_els]
                discarded_read_counts = [('Unaligned',discarded_reads_unaligned),('Duplicate read',discarded_reads_duplicates),('Not supported by R2',discarded_reads_not_supported_by_R2),('Poor alignment',discarded_reads_poor_alignment)]
                classification_categories = fin.readline().rstrip('\n').split("\t")
                classification_category_counts_str_els = fin.readline().rstrip('\n').split("\t")
                classification_category_counts = [int(x) for x in classification_category_counts_str_els]
                classification_read_counts = list(zip(classification_categories, classification_category_counts))

                classification_indel_categories = fin.readline().rstrip('\n').split("\t")
                classification_indel_category_counts_str_els = fin.readline().rstrip('\n').split("\t")
                classification_indel_category_counts = [int(x) for x in classification_indel_category_counts_str_els]
                classification_indel_read_counts = list(zip(classification_indel_categories, classification_indel_category_counts))

                previous_total_read_count = final_total_count
                cut_classification_lookup_str = fin.readline().rstrip('\n')
                cut_classification_lookup = json.loads(cut_classification_lookup_str)

                #load plots
                chr_aln_plot_obj_str = fin.readline().rstrip('\n')
                chr_aln_plot_obj = None
                if chr_aln_plot_obj_str != "" and chr_aln_plot_obj_str != "None":
                    chr_aln_plot_obj = PlotObject.from_json(chr_aln_plot_obj_str)

                tlen_plot_obj_str = fin.readline().rstrip('\n')
                tlen_plot_obj = None
                if tlen_plot_obj_str != "" and tlen_plot_obj_str != "None":
                    tlen_plot_obj = PlotObject.from_json(tlen_plot_obj_str)

                deduplication_plot_obj_str = fin.readline().rstrip('\n')
                deduplication_plot_obj = None
                if deduplication_plot_obj_str != "" and deduplication_plot_obj_str != "None":
                    deduplication_plot_obj = PlotObject.from_json(deduplication_plot_obj_str)

                tx_order_plot_obj_str = fin.readline().rstrip('\n')
                tx_order_plot_obj = None
                if tx_order_plot_obj_str != "" and tx_order_plot_obj_str != "None":
                    tx_order_plot_obj = PlotObject.from_json(tx_order_plot_obj_str)

                tx_count_plot_obj_str = fin.readline().rstrip('\n')
                tx_count_plot_obj = None
                if tx_count_plot_obj_str != "" and tx_count_plot_obj_str != "None":
                    tx_count_plot_obj = PlotObject.from_json(tx_count_plot_obj_str)

                tx_circos_plot_obj_str = fin.readline().rstrip('\n')
                tx_circos_plot_obj = None
                if tx_circos_plot_obj_str != "" and tx_circos_plot_obj_str != "None":
                    tx_circos_plot_obj = PlotObject.from_json(tx_circos_plot_obj_str)

                classification_plot_obj_str = fin.readline().rstrip('\n')
                classification_plot_obj = None
                if classification_plot_obj_str != "" and classification_plot_obj_str != "None":
                    classification_plot_obj = PlotObject.from_json(classification_plot_obj_str)

                classification_indel_plot_obj_str = fin.readline().rstrip('\n')
                classification_indel_plot_obj = None
                if classification_indel_plot_obj_str != "" and classification_indel_plot_obj_str != "None":
                    classification_indel_plot_obj = PlotObject.from_json(classification_indel_plot_obj_str)

                origin_indel_hist_plot_obj_str = fin.readline().rstrip('\n')
                origin_indel_hist_plot_obj = None
                if origin_indel_hist_plot_obj_str != "" and origin_indel_hist_plot_obj_str != "None":
                    origin_indel_hist_plot_obj = PlotObject.from_json(origin_indel_hist_plot_obj_str)

                origin_inversion_hist_plot_obj_str = fin.readline().rstrip('\n')
                origin_inversion_hist_plot_obj = None
                if origin_inversion_hist_plot_obj_str != "" and origin_inversion_hist_plot_obj_str != "None":
                    origin_inversion_hist_plot_obj = PlotObject.from_json(origin_inversion_hist_plot_obj_str)

                origin_deletion_hist_plot_obj_str = fin.readline().rstrip('\n')
                origin_deletion_hist_plot_obj = None
                if origin_deletion_hist_plot_obj_str != "" and origin_deletion_hist_plot_obj_str != "None":
                    origin_deletion_hist_plot_obj = PlotObject.from_json(origin_deletion_hist_plot_obj_str)

                origin_indel_depth_plot_obj_str = fin.readline().rstrip('\n')
                origin_indel_depth_plot_obj = None
                if origin_indel_depth_plot_obj_str != "" and origin_indel_depth_plot_obj_str != "None":
                    origin_indel_depth_plot_obj = PlotObject.from_json(origin_indel_depth_plot_obj_str)

                r1_r2_support_plot_obj_str = fin.readline().rstrip('\n')
                r1_r2_support_plot_obj = None
                if r1_r2_support_plot_obj_str != "" and r1_r2_support_plot_obj_str != "None":
                    r1_r2_support_plot_obj = PlotObject.from_json(r1_r2_support_plot_obj_str)
                r1_r2_support_dist_plot_obj_str = fin.readline().rstrip('\n')
                r1_r2_support_dist_plot_obj = None
                if r1_r2_support_dist_plot_obj_str != "" and r1_r2_support_dist_plot_obj_str != "None":
                    r1_r2_support_dist_plot_obj = PlotObject.from_json(r1_r2_support_dist_plot_obj_str)

                discarded_reads_plot_obj_str = fin.readline().rstrip('\n')
                discarded_reads_plot_obj = None
                if discarded_reads_plot_obj_str != "" and discarded_reads_plot_obj_str != "None":
                    discarded_reads_plot_obj = PlotObject.from_json(discarded_reads_plot_obj_str)

                if previous_total_read_count > -1:
                    #read final assignments from file
                    final_read_ids_for_crispresso = defaultdict(int)
                    final_cut_counts = defaultdict(int)
                    with open (final_file,'r') as fin:
                        head_line = fin.readline().rstrip('\n')
                        head_line_els = head_line.split("\t")
                        read_total_read_count = 0
                        final_cut_ind = 10
                        final_cut_indel_ind = 11
                        if head_line_els[final_cut_ind] == "final_cut_pos" and head_line_els[final_cut_indel_ind] == 'final_cut_indel':
                            logger.info('Reading previously-processed assignments')
                            #if header matches, read file
                            for line in fin:
                                read_total_read_count += 1
                                line_els = line.rstrip("\r\n").split("\t")
                                final_cut_indel = int(line_els[final_cut_indel_ind]) #how long this indel was 
                                final_cut_counts[line_els[final_cut_ind]] += 1
                                if abs(final_cut_indel) <= crispresso_max_indel_size:
                                    final_read_ids_for_crispresso[line_els[0]] = line_els[final_cut_ind]

                        if final_total_count == read_total_read_count:
                            logger.info('Using previously-processed assignments for ' + str(read_total_read_count) + ' total reads')
                            return (final_file,final_read_ids_for_crispresso,final_cut_counts,cut_classification_lookup,final_total_count,discarded_read_counts,classification_read_counts,classification_indel_read_counts,
                                chr_aln_plot_obj,tlen_plot_obj,deduplication_plot_obj,tx_order_plot_obj,tx_count_plot_obj,tx_circos_plot_obj,classification_plot_obj,classification_indel_plot_obj,
                                origin_indel_hist_plot_obj,origin_inversion_hist_plot_obj,origin_deletion_hist_plot_obj,origin_indel_depth_plot_obj,
                                r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,discarded_reads_plot_obj)
                        else:
                            logger.info('In attempting to recover previously-processed assignments, expecting ' + str(previous_total_read_count) + ' reads, but only read ' + str(read_total_read_count) + ' from ' + final_file + '. Reprocessing.')

    logger.info('Analyzing reads aligned to genome')

    # First pass through alignment file:
    # -track cut sites that will be later aggregated into final cut sites
    # -deduplicate based on UMI and r1/r2 alignment
    # -classify reads based on r1/r2 support

    #strings for r1/r2 support (the text is changed in just one place: here)
    r1_r2_support_str_supported = "Supported by R2"
    r1_r2_support_str_not_supported_orientation = "Not supported by R2 - different R1/R2 orientations"
    r1_r2_support_str_not_supported_distance = "Not supported by R2 - alignment separation more than maximum"
    r1_r2_support_str_not_supported_diff_chrs = "Not supported by R2 - different chrs"
    r1_r2_support_str_not_supported_not_aln = "Not supported by R2 - unaligned"

    #analyze alignments
    seen_reads_for_dedup = {} # put read id into this dict for deduplicating
    seen_read_counts_for_dedup = {} # how many times each umi/loc was seen
    total_reads_processed = 0 # includes pairs
    total_r1_processed = 0
    discarded_reads_unaligned = 0 #count of pairs
    discarded_reads_duplicates = 0
    discarded_reads_poor_alignment = 0
    discarded_reads_not_supported_by_R2 = 0
    cut_points_by_chr = defaultdict(lambda: defaultdict(int)) # dict of cut points for each chr contianing observed cuts cut_points_by_chr[chr1][1559988] = count seen
    aligned_tlens = defaultdict(int)
    aligned_chr_counts = defaultdict(int)
    aln_pos_by_chr = defaultdict(lambda: defaultdict(int))

    #keep track of distance between R1/R2
    r1_r2_support_distances = defaultdict(int) #for reads that are separated by less than r1_r2_support_max_distance
    r1_r2_not_support_distances = defaultdict(int) #for reads that are separated by more than r1_r2_support_max_distance

    discarded_read_file = root+'.discarded_read_info.txt'

    if write_discarded_read_info:
        fout_discarded = open(discarded_read_file,'w')

    #counts for 'support' status
    r1_r2_support_status_counts = defaultdict(int)
    final_file_tmp = root + '.final_assignments.tmp'
    with open(final_file_tmp,'w') as af1:
        af1.write("\t".join([str(x) for x in ['read_id','curr_position','cut_pos','cut_direction','del_primer','ins_target','insert_size','r1_r2_support_status','r1_r2_support_dist','r1_r2_orientation_str']])+"\n")

        if not os.path.exists(genome_mapped_bam):
            raise Exception('Cannot find bam file at ' + genome_mapped_bam)
        for line in read_command_output('%s view %s'%(samtools_command,genome_mapped_bam)):
            if line.strip() == "": break

            line_els = line.split("\t")

            total_reads_processed += 1
            read_has_multiple_segments = int(line_els[1]) & 0x1
            read_is_paired_read_1 = int(line_els[1]) & 0x40 # only set with bowtie if read is from paired reads

            if read_has_multiple_segments and not read_is_paired_read_1:
                continue
            total_r1_processed += 1

            read_is_not_primary_alignment = int(line_els[1]) & 0x100
            if read_is_not_primary_alignment:
                discarded_reads_unaligned += 1
                if write_discarded_read_info:
                    fout_discarded.write(line_els[0]+'\tUnaligned\t'+line_els[1]+"\n")
                continue

            line_unmapped = int(line_els[1]) & 0x4
            if line_unmapped:
                discarded_reads_unaligned += 1
                if write_discarded_read_info:
                    fout_discarded.write(line_els[0]+'\tUnmapped\t'+line_els[1]+"\n")
                continue

            (line_info, primer_trimmed_len_str) = line_els[0].split(" primer_len=")
            primer_trimmed_len = int(primer_trimmed_len_str)

            left_matches = 0
            right_matches = 0
            mez_val = ""
            for line_el in line_els:
                if line_el.startswith('MD:Z'):
                    mez_val = line_el[5:]
                    left_matches,right_matches = getLeftRightMismatchesMDZ(mez_val)
                    break

            cigar_str = line_els[5]
            start_clipped = 0
            m = re.match(r"^(\d+)S\d",cigar_str)
            if m:
                start_clipped = int(m.group(1))
            end_clipped = 0
            m = re.match(r"(\d+)S$",cigar_str)
            if m:
                end_clipped = int(m.group(1))

            seq = line_els[9]
            is_rc = int(line_els[1]) & 0x10

            if not suppress_poor_alignment_filter and \
                    ((start_clipped > arm_max_clipped_bases and end_clipped > arm_max_clipped_bases) or \
                    left_matches < arm_min_matched_start_bases or \
                    right_matches < arm_min_matched_start_bases):

                discarded_reads_poor_alignment += 1
                if write_discarded_read_info:
                    fout_discarded.write(line_els[0]+'\tPoor alignment\tclipped: ['+str(start_clipped)+','+str(end_clipped)+'] matches: ['+str(left_matches)+','+str(right_matches)+']\n')
                continue

            line_chr = line_els[2]
            line_mapq = line_els[5]
            line_start = int(line_els[3])
            seq_len = len(seq) - (start_clipped + end_clipped)
            line_end = line_start+(seq_len-1)
            cut_direction = "right"
            if is_rc:
                tmp = line_end
                line_end = line_start
                line_start = tmp
                tmp = start_clipped
                start_clipped = end_clipped
                end_clipped = tmp
                cut_direction = "left"

            curr_position = "%s:%s-%s"%(line_chr,line_start,line_end)
            cut_pos = "%s:%s"%(line_chr,line_start)

            del_primer = len(origin_seq) - primer_trimmed_len #how many bases are deleted from the full primer
            ins_target = start_clipped # how many bases were clipped from alignment (insertions at target)

#            print(f'{line_els[0]=}')
#            print(f'{del_primer=}')
#            print(f'{len(origin_seq)=}')
#            print(f'{primer_trimmed_len=}')
#            print(f'{ins_target=}')
#
            # assign r1/r2 support for this read
            r1_r2_support_status = None
            r1_r2_support_dist = None # distance between R1 and R2
            r1_r2_orientation_str = None

            #these variables are for assignment of this read
            r1_aln_chr_for_support = None #values to see if r1 and r2 agree
            r1_orientation_for_support = None
            r2_aln_chr_for_support = None
            r2_orientation_for_support = None

            if read_is_paired_read_1:
                r1_aln_chr_for_support = line_chr
                r1_orientation_for_support = "-" if is_rc else "+"
                r2_aln_chr_for_support = line_els[6]
                r2_orientation_for_support = "-" if int(line_els[1]) & 0x20 else "+"
                r1_r2_support_dist = abs(int(line_els[8]))

                #determine R1/R2 support
                if r1_aln_chr_for_support == "*" or r2_aln_chr_for_support == "*":
                    r1_r2_support_status = r1_r2_support_str_not_supported_not_aln
                elif r2_aln_chr_for_support == "=" or r1_aln_chr_for_support == r2_aln_chr_for_support:
                    #distance must be within cutoff
                    if r1_r2_support_dist <= r1_r2_support_max_distance:
                        #must have different orientations (R1 is pos, R2 is rc/egative)
                        #for standard illumina sequencing, R1 is forward and R2 is reverse-complement
                        if r1_orientation_for_support != r2_orientation_for_support:
                            r1_r2_support_status = r1_r2_support_str_supported
                        else:
                            r1_r2_support_status = r1_r2_support_str_not_supported_orientation
                    else:
                        r1_r2_support_status = r1_r2_support_str_not_supported_distance
                else:
                    r1_r2_support_status = r1_r2_support_str_not_supported_diff_chrs

                r1_r2_orientation_str = '%s/%s'%(r1_orientation_for_support,r2_orientation_for_support)
                r1_r2_support_status_counts[r1_r2_support_status] += 1
                if r1_r2_support_status == r1_r2_support_str_supported:
                    r1_r2_support_distances[int(r1_r2_support_dist)] += 1
                elif r1_r2_support_status == r1_r2_support_str_not_supported_distance:
                    r1_r2_not_support_distances[int(r1_r2_support_dist)] += 1
                if not suppress_r2_support_filter and r1_r2_support_status != r1_r2_support_str_supported:
                    discarded_reads_not_supported_by_R2 += 1
                    if write_discarded_read_info:
                        fout_discarded.write(line_els[0]+'\tNot supported by R2\t'+r1_r2_support_status+'\n')
                    continue
            #done setting r1/r2 support status for paired reads

            #as a last step, deduplicate
            barcode = line_info.split(":")[-1]
            aln1 = line_els[2] + ":" + line_els[3]
            aln2 = line_els[6] + ":" + line_els[7] # if unpaired, this will just be *:0, so we can still use it in our key
            tlen = abs(int(line_els[8]))
            key = "%s %s %s %s"%(barcode,aln1,aln2,tlen)
            if key in seen_reads_for_dedup:
                seen_read_counts_for_dedup[key] += 1
                if not suppress_dedup_on_aln_pos_and_UMI_filter:
                    discarded_reads_duplicates += 1
                    if write_discarded_read_info:
                        fout_discarded.write(line_els[0]+'\tDuplicate\tDuplicate of '+seen_reads_for_dedup[key]+'('+key+')\n')
                    continue
            else:
                seen_reads_for_dedup[key] = line_els[0]
                seen_read_counts_for_dedup[key] = 1
            #finished deduplication by barcode and aln position

            cut_points_by_chr[line_chr][line_start] += 1
            insert_size = 'None'
            if read_is_paired_read_1: #only set if paired
                insert_size = abs(int(line_els[8]))
                aligned_tlens[insert_size] += 1

            aligned_chr_counts[line_chr] += 1
            aln_pos_window = int(line_start)/genome_map_resolution
            aln_pos_by_chr[line_chr][aln_pos_window] += 1

            af1.write("\t".join([str(x) for x in [line_info,curr_position,cut_pos,cut_direction,del_primer,ins_target,insert_size, r1_r2_support_status, r1_r2_support_dist, r1_r2_orientation_str]]) + "\n")

        #done iterating through bam file
    #close assignments file
    logger.debug('Finished first pass through alignments')
    if write_discarded_read_info:
        fout_discarded.close()

    chr_aln_plot_root = root + ".chr_alignments"
    keys = sorted(aligned_chr_counts.keys())
    vals = [aligned_chr_counts[key] for key in keys]
    with open(chr_aln_plot_root+".txt","w") as chrs:
        chrs.write('chr\tnumReads\n')
        for key in sorted(aligned_chr_counts.keys()):
            chrs.write(key + '\t' + str(aligned_chr_counts[key]) + '\n')

    if suppress_plots:
        chr_aln_plot_obj = None
    else:
        fig = plt.figure(figsize=(12,6))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(range(len(keys)),vals,tick_label=keys)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Location of read alignments')
        plt.savefig(chr_aln_plot_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(chr_aln_plot_root+".png",pad_inches=1,bbox_inches='tight')

        plot_label = 'Bar plot showing alignment location of read alignments'
        if len(keys) == 0:
            plot_label = '(No reads aligned to the genome)'

        chr_aln_plot_obj = PlotObject(
                plot_name = chr_aln_plot_root,
                plot_title = 'Genome alignment summary',
                plot_label = plot_label,
                plot_datas = [('Genome alignment summary',chr_aln_plot_root + ".txt")]
                )

    tlen_plot_obj = None
    if len(aligned_tlens) > 0: #paired-end reads
        tlen_plot_root = root + ".insertSizes"
        keys = sorted(aligned_tlens.keys())
        vals = [aligned_tlens[key] for key in keys]
        with open(tlen_plot_root+".txt","w") as fout:
            fout.write('insertSize\tnumReads\n')
            for key in keys:
                fout.write(str(key) + '\t' + str(aligned_tlens[key]) + '\n')

        if not suppress_plots:
            fig = plt.figure(figsize=(12,6))
            ax = plt.subplot(111)
            if len(vals) > 0:
                ax.bar(keys,vals)
                ax.set_ymargin(0.05)
            else:
                ax.bar(0,0)
            ax.set_ylabel('Number of Reads')
            ax.set_title('Insert size')
            plt.savefig(tlen_plot_root+".pdf",pad_inches=1,bbox_inches='tight')
            plt.savefig(tlen_plot_root+".png",pad_inches=1,bbox_inches='tight')

            tlen_plot_obj = PlotObject(
                    plot_name = tlen_plot_root,
                    plot_title = 'Genome Alignment Insert Size Summary',
                    plot_label = 'Bar plot showing insert size of reads aligned to genome',
                    plot_datas = [('Genome alignment insert size summary',tlen_plot_root + ".txt")]
                    )

    #deduplication plot

    dup_counts = defaultdict(int)
    dup_keys = seen_read_counts_for_dedup.keys()
    for dup_key in dup_keys:
        dup_counts[seen_read_counts_for_dedup[dup_key]] += 1
    dup_count_keys = sorted(dup_counts.keys())
    dup_count_vals = [dup_counts[key] for key in dup_count_keys]

    deduplication_plot_obj_root = root + ".deduplication_by_UMI_and_aln_pos"
    with open(deduplication_plot_obj_root+".txt","w") as summary:
        summary.write('Reads per UMI\tNumber of UMIs\n')
        for key in dup_count_keys:
            summary.write(str(key) + '\t' + str(dup_counts[key]) + '\n')

    deduplication_plot_obj = None
    if total_reads_processed > 0 and not suppress_plots:
        labels = ['Not duplicate','Duplicate']
        values = [total_r1_processed-discarded_reads_duplicates,discarded_reads_duplicates]

        fig = plt.figure(figsize=(12,6))
        ax = plt.subplot(121)
        pie_values = []
        pie_labels = []
        for i in range(len(labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Duplicate read counts')

        ax2 = plt.subplot(122)
        if len(dup_count_vals) > 0:
            ax2.bar(dup_count_keys,dup_count_vals)
            ax2.set_ymargin(0.05)
        else:
            ax2.bar(0,0)
        ax2.set_ylabel('Number of UMIs')
        ax2.set_title('Reads per UMI')

        plt.savefig(deduplication_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(deduplication_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        dedup_note = ". Note that deduplication of reads by alignment position and UMI is performed by default. To disable deduplication, set the parameter --suppress_dedup_on_aln_pos_and_UMI_filter."
        if suppress_dedup_on_aln_pos_and_UMI_filter:
            dedup_note = ". Deduplication based on alignment and UMI position is not performed because the flag --suppress_dedup_on_aln_pos_and_UMI_filter is set."

        deduplication_plot_obj = PlotObject(
                plot_name = deduplication_plot_obj_root,
                plot_title = 'Duplicate read counts',
                plot_label = 'Number of reads that were duplicates based on alignment and UMI' + dedup_note,
                plot_datas = [('Read assignments',final_file)]
                )

    #now that we've aggregated all cut sites, find likely cut positions and record assignments to those positions in final_cut_point_lookup


    cut_point_homology_info = {} #dict containing cut points with sufficient homology
    with open(root+".cut_homology.txt",'w') as fout:
        fout.write("\t".join([str(x) for x in ['guide_seq','cut_chr','cut_point','is_valid_homology_site','potential_guide','best_score','match_direction','n_matches','n_mismatches','n_gaps','aln_guide','aln_ref']])+"\n")
        aln_match_score = 2
        aln_mismatch_score = -1
        aln_gap_extend_score = -2
        aln_gap_open_score = -5
        aln_matrix = CRISPResso2Align.make_matrix(match_score=aln_match_score,mismatch_score=aln_mismatch_score)
        padded_seq_len = 30 #how many bp to extend around cut to search for guide
        guide_homology_max_gaps = 2
        guide_homology_max_mismatches = 5
        for guide_seq in guide_seqs:
            for cut_chr in cut_points_by_chr:
                these_cut_points = sorted(cut_points_by_chr[cut_chr])
                for cut_point in these_cut_points:
                    faidx_cmd = '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,cut_chr,cut_point-padded_seq_len,cut_point+padded_seq_len)
                    logger.debug('Searching for homology. Getting sequence around cut site ' + cut_chr + ':' + str(cut_point) + ': ' + str(faidx_cmd))
                    cut_padded_seq = subprocess.check_output(faidx_cmd,shell=True).decode(sys.stdout.encoding).strip()
                    potential_guide = None
                    amp_incentive = np.zeros(len(cut_padded_seq)+1,dtype=int)
                    f1,f2,fw_score=CRISPResso2Align.global_align(guide_seq.upper(),cut_padded_seq.upper(),matrix=aln_matrix,gap_incentive=amp_incentive,gap_open=aln_gap_open_score,gap_extend=aln_gap_extend_score)
                    r1,r2,rv_score=CRISPResso2Align.global_align(reverse_complement(guide_seq.upper()),cut_padded_seq.upper(),matrix=aln_matrix,gap_incentive=amp_incentive,gap_open=aln_gap_open_score,gap_extend=aln_gap_extend_score)

                    is_rc = False
                    best_score = fw_score
                    best_1 = f1
                    best_2 = f2
                    if fw_score > rv_score:
                        n_matches,n_mismatches,n_gaps,cleavage_ind,potential_guide = get_guide_match_from_aln(f1,f2,len(guide_seq)+cleavage_offset-1) #cleavage ind returns
                    else:
                        is_rc = True
                        best_score = rv_score
                        best_1 = r1
                        best_2 = r2
                        n_matches,n_mismatches,n_gaps,cleavage_ind,potential_guide = get_guide_match_from_aln(r1,r2,-1*cleavage_offset-1) #cleavage ind returns 
                        potential_guide = reverse_complement(potential_guide)

                    cut_key = cut_chr+":"+str(cut_point-padded_seq_len+cleavage_ind+1)
                    is_rc_string = 'Reverse' if is_rc else 'Forward'
                    is_valid_homology = False

                    if n_gaps <= guide_homology_max_gaps and n_mismatches <= guide_homology_max_mismatches:
                        is_valid_homology = True
                        cut_key = cut_chr+":"+str(cut_point-padded_seq_len+cleavage_ind+1)
                        if cut_key not in cut_point_homology_info:
                            cut_point_homology_info[cut_key] = (best_score,potential_guide,n_matches,n_mismatches,n_gaps,is_rc)
                        elif best_score > cut_point_homology_info[cut_key][0]: #sometimes a suboptimal alignment (guide is aligned partially beyond the amplicon) could already identify this cut, so take the one with the better alignment score
                            cut_point_homology_info[cut_key] = (best_score,potential_guide,n_matches,n_mismatches,n_gaps,is_rc)

                    fout.write("\t".join([str(x) for x in [guide_seq,cut_chr,cut_point,is_valid_homology,potential_guide,best_score,is_rc_string,n_matches,n_mismatches,n_gaps,best_1,best_2]])+"\n")

    #the counts in cut_points include only reads that weren't marked as duplicates
    #final_cut_point_lookup = dict from old(fuzzy/imprecise) position to new
    #final_cut_points_by_chr = dict of chr > cut points
    origin_chr, origin_cut_pos, origin_direction, final_cut_point_lookup, final_cut_points_by_chr, found_cut_point_total = collapse_cut_points(
        novel_cut_merge_distance = novel_cut_merge_distance,
        known_cut_merge_distance = known_cut_merge_distance,
        origin_cut_merge_distance = origin_cut_merge_distance,
        cut_points_by_chr = cut_points_by_chr,
        cut_annotations = cut_annotations,
        cut_point_homology_info = cut_point_homology_info
        )

    final_cut_point_total = 0
    cut_classification_lookup = {}
    for chrom in sorted(final_cut_points_by_chr.keys()):
        final_cut_point_total += len(final_cut_points_by_chr[chrom])
        for pos in sorted(final_cut_points_by_chr[chrom]):
            this_anno=['Novel','Not-origin','NA']
            this_str = chrom+":"+str(pos)
            if this_str in cut_point_homology_info:
                this_anno = ['Homologous','Not-origin','NA']
            if this_str in cut_annotations: #overwrite homology annotation if present in cut_annotations
                this_anno = cut_annotations[this_str]
            if origin_direction is None: # guide seq primer
                cut_classification_lookup['%s:%s:%s'%(chrom,pos,'left')] = this_anno[0]
                cut_classification_lookup['%s:%s:%s'%(chrom,pos,'right')] = this_anno[0] 
            else: #origin is genomic
                if chrom != origin_chr:
                    cut_classification_lookup['%s:%s:%s'%(chrom,pos,'left')] = this_anno[0]+' translocation'
                    cut_classification_lookup['%s:%s:%s'%(chrom,pos,'right')] = this_anno[0] + ' translocation'
                else:
                    if origin_direction == 'left':
                        if pos == origin_cut_pos:
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'right')] = 'Linear'
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'left')] = 'Chimera'
                        else:
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'left')] = this_anno[0] + ' large deletion'
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'right')] = this_anno[0] + ' large inversion'

                    else: # origin direction is right
                        if pos == origin_cut_pos:
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'left')] = 'Linear'
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'right')] = 'Chimera'
                        else:
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'left')] = this_anno[0] + ' large deletion'
                            cut_classification_lookup['%s:%s:%s'%(chrom,pos,'right')] = this_anno[0] + ' large inversion'
    logger.info(str(final_cut_point_total) + ' final cut points after collapsing')


    #add cut classifications from parameters
    for cut_classification_param in cut_classification_annotations:
        cut_classification_els = cut_classification_param.split(":")
        cut_key = ":".join(cut_classification_els[0:3])
        cut_val = cut_classification_els[3]
        cut_classification_lookup[cut_key] = cut_val

    with open(root+".cut_classification.txt",'w') as fout:
        for key in sorted(cut_classification_lookup.keys()):
            fout.write(key+'\t'+cut_classification_lookup[key]+"\n")

    with open(root + ".final_cut_points.txt",'w') as fout:
        fout.write('chr\tpos\tcount\tannotation\n')
        for chrom in sorted(final_cut_points_by_chr.keys()):
            for pos in sorted(final_cut_points_by_chr[chrom]):
                this_anno=['Novel','Not-origin','NA']
                this_str = chrom+":"+str(pos)
                if this_str in cut_point_homology_info:
                    this_anno = ['Homologous','Not-origin','NA']
                if this_str in cut_annotations: #overwrite homology annotation if present in cut_annotations
                    this_anno = cut_annotations[this_str]
                fout.write("%s\t%d\t%d\t%s\n"%(chrom,pos,final_cut_points_by_chr[chrom][pos],",".join(this_anno)))

    with open(root+".final_cut_points_lookup.txt",'w') as fout:
        fout.write('discovered cut point\tfinal mapped cut point\tcount\tannotation\n');
        for key in sorted(final_cut_point_lookup.keys(),key=lambda k: (k.split(":")[0],int(k.split(":")[1]))):
            (chrom,pos)=key.split(":")
            this_anno=['Novel','Not-origin','NA']
            this_final = final_cut_point_lookup[key]
            if this_final in cut_point_homology_info:
                this_anno = ['Homologous','Not-origin','NA']
            if this_final in cut_annotations: # overwrite homology annotation if present in cut_annotations
                this_anno = cut_annotations[this_final]
            fout.write("%s\t%s\t%d\t%s\n"%(key,final_cut_point_lookup[key],cut_points_by_chr[chrom][int(pos)],",".join(this_anno)))

    final_total_count = 0
    final_classification_counts = defaultdict(int) # classification is Linear, Large Deletion, Translocation etc.
    final_classification_indel_counts = defaultdict(int) # classification is Linear, Large Deletion, Translocation etc., this adds '_indels' for possible classes

    #keep track of read ids aligning to cuts
    #reads with these ids will be pulled out of the fastq for analysis by CRISPResso
    #cut is chr:pos:direction
    final_read_ids_for_crispresso = defaultdict(int) # dict readID->cut
    final_cut_counts = defaultdict(int) # dict cut->count of reads

    if origin_direction == 'left':
        origin_direction_opposite = "right"
    else:
        origin_direction_opposite = "left"
    origin_indel_lengths = defaultdict(int) # dict indel_len -> count
    origin_inversion_lengths = defaultdict(int) # dict inversion_len -> count of reads on same chr as primer and same orientation with specified distance between origin and read (primer is to the genomic left of the cut site, this read is oriented to the genomic left side)
    origin_deletion_lengths = defaultdict(int) # dict deletion_len -> count of reads on same chr as primer and opposite direction with specified distance between origin and read (primer is to the genomic left of the cut site, this read is oriented to the genomic right side)
    origin_depth_counts_100bp = np.zeros(101) #100bp count the number of reads with deletions covering each position
    origin_depth_counts_100bp_primer = np.zeros(len(origin_seq)+1) #count deletions covering each base of primer (for those with deletions within 100bp)
    origin_depth_counts_100bp_total = 0 #count of total reads

    long_distance = 100000
    long_distance_str = '100kb'
    origin_depth_counts_long = np.zeros(long_distance + 1) #long count of the number of reads with deletions covering that position
    origin_depth_counts_long_primer = np.zeros(len(origin_seq)+1) #count reads covering each base of primer (for those with deletions within {long_distance})
    origin_depth_counts_long_total = 0 #count of total reads

    read_id_ind = 0
    aln_ind = 1
    cut_pos_ind = 2
    cut_direction_ind = 3
    del_primer_ind = 4
    ins_target_ind = 5
    insert_size_ind = 6
    r1_r2_support_status_ind = 7
    r1_r2_support_dist_ind = 8
    r1_r2_orientation_str_ind = 9

    with open(final_file_tmp,'r') as fin, open(final_file,'w') as fout:
        old_head = fin.readline().strip()
        fout.write(old_head + "\t" + "\t".join([str(x) for x in ['final_cut_pos','final_cut_indel','classification']]) +"\n")
        for line in fin:
            line = line.rstrip('\n')
            final_total_count += 1
            line_els = line.split("\t")

            line_id = line_els[read_id_ind]
            aln_loc = line_els[aln_ind]
            cut_point = line_els[cut_pos_ind]
            cut_point_direction = line_els[cut_direction_ind]
            final_cut_point = final_cut_point_lookup[cut_point]
            final_cut_point_and_direction = final_cut_point + ":" + cut_point_direction

            (final_cut_point_chr,final_cut_point_start) = final_cut_point.split(":")
            final_cut_point_start = int(final_cut_point_start)
            aln_start = int(aln_loc.split(":")[1].split("-")[0])

            #read extends to the right of the cut point
            if cut_point_direction == "right":
                aln_indel = -1*(final_cut_point_start - aln_start)
            else: #read extends to the left of the cut point
                aln_indel = final_cut_point_start - aln_start - 1

            del_primer = int(line_els[del_primer_ind])
            ins_target = int(line_els[ins_target_ind])

            final_cut_indel = ins_target - del_primer - aln_indel

            final_cut_counts[final_cut_point_and_direction] += 1
            if abs(final_cut_indel) <= crispresso_max_indel_size:
                final_read_ids_for_crispresso[line_id] = final_cut_point_and_direction

            final_classification = cut_classification_lookup[final_cut_point_and_direction]
            final_classification_counts[final_classification] += 1

            indel_str = ' no indels'
            if abs(final_cut_indel) > 0:
                if abs(final_cut_indel) > short_indel_length_cutoff:
                    indel_str = ' long indels'
                else:
                    indel_str = ' short indels'

            final_classification_indel_counts[final_classification+indel_str] += 1

            fout.write(line + "\t" + "\t".join([str(x) for x in [final_cut_point_and_direction,final_cut_indel,final_classification]]) +"\n")
#            fout.write(line + "\t" + "\t".join([str(x) for x in [final_cut_point_and_direction,'it:'+str(ins_target)+';dp:'+str(del_primer)+';ai:'+str(aln_indel)+';fci:'+str(final_cut_indel),final_classification]]) +"\n")

            if final_cut_point_chr == origin_chr:
                if final_cut_point_start == origin_cut_pos:
                    if cut_point_direction == origin_direction_opposite:
                        origin_indel_lengths[final_cut_indel] += 1
                this_dist_to_origin = aln_start - origin_cut_pos
                if cut_point_direction == origin_direction:
                    origin_inversion_lengths[this_dist_to_origin] += 1
                else:
                    origin_deletion_lengths[this_dist_to_origin] += 1

                if origin_direction == 'left' and cut_point_direction == 'right' and aln_start >= origin_cut_pos: #  primer--origin_seq--->cut..aln_start (origin extends to the left from cut)
                    this_del_len = aln_start - origin_cut_pos
                    if this_del_len <= 100:
                        origin_depth_counts_100bp_total += 1
                        origin_depth_counts_100bp[0:this_del_len+1] += 1
                        origin_depth_counts_100bp_primer[0:del_primer+1] += 1
                    if this_del_len <= long_distance:
                        origin_depth_counts_long_total += 1
                        origin_depth_counts_long[0:this_del_len+1] += 1
                        origin_depth_counts_long_primer[0:del_primer+1] += 1
                elif origin_direction == 'right' and cut_point_direction == 'left' and aln_start <= origin_cut_pos: #  aln_start..cut<--origin_seq--primer (origin extends to the right from cut)
                    this_del_len = (aln_start - origin_cut_pos) * -1
                    if this_del_len <= 100:
                        origin_depth_counts_100bp_total += 1
                        origin_depth_counts_100bp[0:this_del_len+1] += 1
                        origin_depth_counts_100bp_primer[0:del_primer+1] += 1
                    if this_del_len <= long_distance:
                        origin_depth_counts_long_total += 1
                        origin_depth_counts_long[0:this_del_len+1] += 1
                        origin_depth_counts_long_primer[0:del_primer+1] += 1

            #finish iterating through read ids
        #close final file
    logger.info('Processed %d reads.'%(final_total_count))

    tx_keys = sorted(final_cut_counts.keys())
    tx_list_report = root + ".translocation_list.txt"
    written_tx_keys = {}
    with open (tx_list_report,'w') as fout:
        fout.write('cut\tcut_annotation\tcount_left\tcount_right\tcount_total\n')
        for tx_key in tx_keys:
            this_cut_annotation = 'Novel'
            tx_cut = ':'.join(tx_key.split(":")[0:2])
            if tx_cut in written_tx_keys:
                continue
            if tx_cut in cut_annotations:
                this_cut_annotation = cut_annotations[tx_cut][0]
            this_left_count = final_cut_counts[tx_cut+":left"]
            this_right_count = final_cut_counts[tx_cut+":right"]
            this_total_count = this_left_count + this_right_count
            fout.write("%s\t%s\t%s\t%s\t%s\n"%(tx_cut,this_cut_annotation,this_left_count,this_right_count,this_total_count))
            written_tx_keys[tx_cut] = 1
    logger.debug('Wrote translocation list ' + tx_list_report)

    # fragment translocations
    tx_list_report_root = root + '.fragment_translocation_list'
    tx_list_report = tx_list_report_root + '.txt'
    with open (tx_list_report,'w') as fout:
        fout.write('fragment\tcut_annotation\tfragment_annotation\tcount\n')
        for tx_key in tx_keys:
            cut_annotation = 'Novel'
            tx_cut = ':'.join(tx_key.split(":")[0:2])
            if tx_cut in cut_annotations:
                cut_annotation = cut_annotations[tx_cut][0]
            tx_classification = 'Unknown'
            if tx_key in cut_classification_lookup:
                tx_classification = cut_classification_lookup[tx_key]
            fout.write("%s\t%s\t%s\t%s\n"%(tx_key,cut_annotation,tx_classification,final_cut_counts[tx_key]))
        logger.debug('Wrote fragment translocation list ' + tx_list_report)

    sorted_tx_list = sorted(final_cut_counts, key=final_cut_counts.get,reverse=True)
    top_number_to_plot = 20
    top_sorted_tx_list = sorted_tx_list[:min(top_number_to_plot,len(sorted_tx_list))]
    other_tx_count = 0 #number of locations not plotted
    other_tx_read_count = 0 # number of reads to other locations not plotted
    if len(sorted_tx_list) > top_number_to_plot:
        for i in range(top_number_to_plot,len(sorted_tx_list)):
            other_tx_count += 1
            other_tx_read_count += final_cut_counts[sorted_tx_list[i]]

    tx_order_plot_obj = None
    tx_count_plot_obj = None
    tx_circos_plot_obj = None
    classification_plot_obj = None
    classification_indel_plot_obj = None

    classification_labels = []
    if origin_chr is not None:
        classification_labels = ['Linear'] #linear goes first in order
    for key in sorted(final_classification_counts.keys()):
        if key not in classification_labels and final_classification_counts[key] > 0:
            classification_labels.append(key)

    classification_values = [final_classification_counts[x] for x in classification_labels]
    classification_read_counts_str = "\t".join(classification_labels)+"\n"+"\t".join([str(x) for x in classification_values])+"\n"
    classification_read_counts = list(zip(classification_labels,classification_values))
    classification_plot_obj_root = root + ".classifications"
    with open(classification_plot_obj_root+".txt",'w') as summary:
        summary.write(classification_read_counts_str)

    pie_values = []
    pie_labels = []
    for i in range(len(classification_labels)):
        pie_values.append(classification_values[i])
        pie_labels.append(classification_labels[i]+"\n("+str(classification_values[i])+")")
    noindel_pie_values = pie_values #for use below

    short_indel_categories = [' no indels',' short indels',' long indels']
    classification_indel_labels = []
    classification_indel_counts = []
    inner_pie_values = []
    inner_pie_labels = []
    outer_pie_values = []
    inner_other_count = 0
    outer_other_counts = [0]*len(short_indel_categories)
    sum_inner = sum(noindel_pie_values)
    cutoff_pct = 0.05 #counts < this percent are shown as 'other'
    for label in classification_labels:
        this_counts = []
        for short_indel_category in short_indel_categories:
            new_label = label + short_indel_category
            classification_indel_labels.append(new_label)
            classification_indel_counts.append(final_classification_indel_counts[new_label])
            this_counts.append(final_classification_indel_counts[new_label])

        this_sum = sum(this_counts)
        if sum_inner > 0:
            if this_sum/sum_inner > cutoff_pct:
                inner_pie_labels.append(label+"\n("+str(final_classification_counts[label])+")")
                inner_pie_values.append(final_classification_counts[label])
                outer_pie_values.extend(this_counts)
            else:
                inner_other_count += this_sum

                for i in range(len(short_indel_categories)):
                    outer_other_counts[i] += this_counts[i]

    if inner_other_count > 0:
        inner_pie_labels.append('Other\n('+str(inner_other_count) + ')')
        inner_pie_values.append(inner_other_count)
        outer_pie_values.extend(outer_other_counts)

    classification_indel_read_counts_str = "\t".join(classification_indel_labels)+"\n"+"\t".join([str(x) for x in classification_indel_counts])+"\n"
    classification_indel_read_counts = list(zip(classification_indel_labels,classification_indel_counts))
    classification_indel_plot_obj_root = root + ".classifications_with_indels"
    with open(classification_indel_plot_obj_root+".txt",'w') as summary:
        summary.write("\t".join(classification_indel_labels)+"\n")
        summary.write("\t".join([str(x) for x in classification_indel_counts])+"\n")
    if len(top_sorted_tx_list) > 0 and not suppress_plots:
        logger.debug('Plotting translocations')
        # make tx order plot
        pos_list = [':'.join(x.split(':')[0:2]) for x in top_sorted_tx_list]
        pos_list = sorted(list(set(pos_list)))

        cut_categories = ['Programmed', 'Off-target', 'Known','Cas-OFFinder','Novel']
        cut_colors = plt.get_cmap('Set1',len(cut_categories))
        cut_category_lookup = {}
        for idx,cat in enumerate(cut_categories):
            cut_category_lookup[cat] = cut_colors(idx/len(cut_categories))

        color_grad = plt.get_cmap('viridis',len(pos_list))
        color_lookup = {}
        for idx,pos in enumerate(pos_list):
            color_lookup[pos] = color_grad(idx/len(pos_list))

        left_labs = []
        right_labs = []
        counts = []
        fill_cols = []
        outline_cols = []
        for tx_order_obj in top_sorted_tx_list:
            left_labs.append(cut_classification_lookup[tx_order_obj])
            right_labs.append(tx_order_obj)
            counts.append(final_cut_counts[tx_order_obj])

            this_chr_pos = ':'.join(tx_order_obj.split(':')[0:2])
            this_fill_col = color_lookup[this_chr_pos]
            this_cut_anno='Novel'
            if this_chr_pos in cut_annotations:
                this_cut_anno = cut_annotations[this_chr_pos][0]

            fill_cols.append(this_fill_col)
            if 'Cas-OFFinder' in this_cut_anno: #casoffinder categories look like "Cas-OFFinder OB 3"
                this_cut_anno = 'Cas-OFFinder'
            this_outline_col = cut_category_lookup[this_cut_anno]
            outline_cols.append(this_outline_col)

        if other_tx_count > 0:
            left_labs.append('Other')
            right_labs.append('('+str(other_tx_count)+' locations)')
            counts.append(other_tx_read_count)

            fill_cols.append('0.8') #light gray
            outline_cols.append('None')

        tx_plt = makeTxCountPlot(left_labs=left_labs, right_labs=right_labs, counts=counts, fill_cols=fill_cols, outline_cols=outline_cols,legend_outline_cols=cut_colors.colors, legend_outline_labs=cut_categories)
        plot_name = tx_list_report_root
        tx_plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
        tx_plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')

        truncated_str = ''
        if len(sorted_tx_list) < len(final_cut_counts):
            truncated_str = '. Note that only the top ' + str(len(sorted_tx_list)) + '/' + str(len(final_cut_counts)) + ' events are shown. All events and counts can be found in the Fragment translocation list.'
        tx_order_plot_obj = PlotObject(plot_name = plot_name,
                plot_title = 'Translocation Summary',
                plot_label = 'Final translocation counts' + truncated_str,
                plot_datas = [('Fragment translocation list',tx_list_report)]
                )

        #make count plot obj
        def get_chr_sort_key(chrom):
            chrom = chrom.replace("chr","")
            if chrom == 'X': return 23
            elif chrom == 'Y': return 24
            elif chrom == 'M': return 25
            elif '_' in chrom : return 26

            else: chrom = int(chrom)
            return(chrom)

        chrom_max_count = 0
        height_at_locs = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        #here only take the top items (not all of them or the plot gets wonky)
        #but we'll take the left and right sides
        for loc_key in top_sorted_tx_list:
            loc_key_els = loc_key.split(":")
            cut_key = ":".join(loc_key_els[0:2])
            height_at_locs[loc_key_els[0]][int(loc_key_els[1])]['left'] = final_cut_counts[cut_key+":"+'left']
            height_at_locs[loc_key_els[0]][int(loc_key_els[1])]['right'] = final_cut_counts[cut_key+":"+'right']
            #loc_key will always be bigger than the opposite (left vs right) cut because it's in the top_sorted_tx_list
            if final_cut_counts[loc_key] > chrom_max_count: chrom_max_count = final_cut_counts[loc_key]

        sorted_chroms = sorted(height_at_locs.keys(), key=lambda x: get_chr_sort_key(x))

        chrom_max_cuts = 0
        for idx,chrom in enumerate(sorted_chroms):
            chrom_locs = sorted(height_at_locs[chrom].keys())
            if len(chrom_locs) > chrom_max_cuts: chrom_max_cuts = len(chrom_locs)

        width = 1/float(chrom_max_cuts)
        width=0.4

        origin_idx = None
        origin_mapped_pos = None
        if len(sorted_chroms) > 1:
            fig, axs = plt.subplots(len(sorted_chroms), 1, figsize=(12,max(12,3*len(sorted_chroms))))
            for idx,chrom in enumerate(sorted_chroms):
                chrom_locs = sorted(height_at_locs[chrom].keys())
                left_vals = []
                right_vals = []
                xs = np.arange(len(chrom_locs))

                for idx2,chrom_loc in enumerate(chrom_locs):
                    left_vals.append(height_at_locs[chrom][chrom_loc]['left'])
                    right_vals.append(height_at_locs[chrom][chrom_loc]['right'])
                    if origin_chr == chrom and origin_cut_pos == chrom_loc:
                        origin_idx = idx
                        origin_mapped_pos = idx2

                axs[idx].bar(xs,left_vals,width=width,label='left')
                axs[idx].bar(xs+width,right_vals,width=width,label='right')
                axs[idx].set_ylabel(chrom)
                axs[idx].set_xticks(xs + width/2)
                #axs[idx].set_xticklabels([chrom+":"+str(x) for x in chrom_locs], rotation = 90)
                axs[idx].set_xticklabels([str(x) for x in chrom_locs], rotation = 90)
                axs[idx].set_xlim(0-width,2*width+chrom_max_cuts-1)
                axs[idx].set_ylim(1,chrom_max_count*2)
                axs[idx].tick_params(axis='both', labelsize=6)
                axs[idx].set_yscale('log')

            fig.subplots_adjust(hspace=1)
            axs[0].set_title('Number of reads aligned at left and right sides of cut sites')

            if origin_idx is not None and origin_mapped_pos is not None:
                if origin_direction == 'left':
                    axs[origin_idx].annotate("Primer", xy=(origin_mapped_pos+width/2, chrom_max_count/10), xytext=(origin_mapped_pos-width/2, chrom_max_count/10),
                        arrowprops=dict(arrowstyle="->"),va='center')
                else:
                    axs[origin_idx].annotate("Primer", xy=(origin_mapped_pos+width/2, chrom_max_count/10), xytext=(origin_mapped_pos+width, chrom_max_count/10),
                        arrowprops=dict(arrowstyle="->"),va='center')

        else: #only one chr to plot
            fig, ax = plt.subplots(figsize=(12,max(12,3*len(sorted_chroms))))
            for idx,chrom in enumerate(sorted_chroms):
                chrom_locs = sorted(height_at_locs[chrom].keys())
                left_vals = []
                right_vals = []
                xs = np.arange(len(chrom_locs))

                for idx2,chrom_loc in enumerate(chrom_locs):
                    left_vals.append(height_at_locs[chrom][chrom_loc]['left'])
                    right_vals.append(height_at_locs[chrom][chrom_loc]['right'])
                    if origin_chr == chrom and origin_cut_pos == chrom_loc:
                        origin_idx = idx
                        origin_mapped_pos = idx2

                ax.bar(xs,left_vals,width=width,label='left')
                ax.bar(xs+width,right_vals,width=width,label='right')
                ax.set_ylabel(chrom)
                ax.set_xticks(xs + width/2)
                #ax.set_xticklabels([chrom+":"+str(x) for x in chrom_locs], rotation = 90)
                ax.set_xticklabels([str(x) for x in chrom_locs], rotation = 90)
                ax.set_xlim(0-width,2*width+chrom_max_cuts-1)
                ax.set_ylim(1,chrom_max_count*2)
                ax.tick_params(axis='both', labelsize=6)
                ax.set_yscale('log')
                ax.set_title('Number of reads aligned at left and right sides of cut sites')

                if origin_idx is not None and origin_mapped_pos is not None:
                    if origin_direction == 'left':
                        ax.annotate("Primer", xy=(origin_mapped_pos+width/2, chrom_max_count/10), xytext=(origin_mapped_pos-width/2, chrom_max_count/10),
                            arrowprops=dict(arrowstyle="->"),va='center')
                    else:
                        ax.annotate("Primer", xy=(origin_mapped_pos+width/2, chrom_max_count/10), xytext=(origin_mapped_pos+width, chrom_max_count/10),
                            arrowprops=dict(arrowstyle="->"),va='center')

        plot_name = tx_list_report_root + ".counts"
        plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')

        truncated_str = ''
        if len(top_sorted_tx_list) < len(final_cut_counts):
            truncated_str = '. Note that only the top ' + str(len(top_sorted_tx_list)) + '/' + str(len(final_cut_counts)) + ' events are shown. All events and counts can be found in the Fragment translocation list.'

        tx_count_plot_obj = PlotObject(plot_name = plot_name,
                plot_title = 'Translocation Counts',
                plot_label = 'Final translocation positions and counts' + truncated_str,
                plot_datas = [('Fragment translocation list',tx_list_report)]
                )

        #plot circos plot
        plot_name = tx_list_report_root + ".circos"
        plot_tx_circos(genome,final_cut_counts,origin_chr,origin_cut_pos,plot_name)
        tx_circos_plot_obj = PlotObject(plot_name = plot_name,
                plot_title = 'Translocation Circos Plot',
                plot_label = 'Final translocation positions and counts',
                plot_datas = [('Fragment translocation list',tx_list_report)]
                )

        #assignment plot
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Read classifications')
        plt.savefig(classification_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(classification_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(classification_labels,classification_values)])
        classification_plot_obj = PlotObject(
                plot_name = classification_plot_obj_root,
                plot_title = 'Read Classification',
                plot_label = 'Read classification<br>'+plot_count_str,
                plot_datas = [
                    ('Alignment classifications',classification_plot_obj_root + ".txt"),
                    ('Read assignments',final_file)
                    ]
                )

        #assignment (with indels) plot
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        width = 0.1
        inner_wedge_properties = {"edgecolor":"w",'linewidth': 2}
        wedge_properties = {"width":width, "edgecolor":"w",'linewidth': 2}

        ax.pie(inner_pie_values, labels=inner_pie_labels, # plot categories without indels
            wedgeprops=inner_wedge_properties,autopct="%1.2f%%",radius=1-width, labeldistance=0.50,pctdistance=0.30)
        ax.pie(outer_pie_values, labels=None,
            wedgeprops=wedge_properties,autopct="%1.2f%%",colors=['silver','pink','crimson'],pctdistance=1.08)
        patch1 = Patch(color='silver', label='No indels')
        patch2 = Patch(color='pink', label='Short indels')
        patch3 = Patch(color='crimson', label='Long indels')

        plt.legend(handles=[patch1, patch2, patch3])

        ax.set_title('Read classifications')
        plt.savefig(classification_indel_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(classification_indel_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(classification_indel_labels,classification_indel_counts)])
        classification_indel_plot_obj = PlotObject(
                plot_name = classification_indel_plot_obj_root,
                plot_title = 'Read Classification including indels',
                plot_label = 'Read classification including annotations for presence of indels<br>(Short indels: <=' +str(short_indel_length_cutoff)+'bp, Long indels: >'+str(short_indel_length_cutoff)+'bp)<br>'+plot_count_str,
                plot_datas = [
                    ('Alignment classifications with indels',classification_indel_plot_obj_root + ".txt"),
                    ('Alignment classifications',classification_plot_obj_root + ".txt"),
                    ('Read assignments',final_file)
                    ]
                )
    #end if len(top_sorted_tx_list) > 0

    origin_indel_hist_plot_obj = None # plot of indel lengths at origin
    origin_indel_hist_plot_obj_root = root + ".origin_indel_histogram"
    keys = sorted(origin_indel_lengths.keys())
    vals = [origin_indel_lengths[key] for key in keys]
    with open(origin_indel_hist_plot_obj_root+".txt","w") as summary:
        summary.write('Indel_length\tnumReads\n')
        for key in keys:
            summary.write(str(key) + '\t' + str(origin_indel_lengths[key]) + '\n')
    if len(origin_indel_lengths) > 0 and not suppress_plots:
        fig = plt.figure(figsize=(12,6))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(keys,vals)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Indel length at origin')
        plt.savefig(origin_indel_hist_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(origin_indel_hist_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_label = 'Bar plot showing lengths of deletions (negative) or insertions (positive) for linear reads supporting the origin cut site at ' + origin_chr + ': ' + str(origin_cut_pos) 
        if len(keys) == 0:
            plot_label = '(No indel lengths)'

        origin_indel_hist_plot_obj = PlotObject(
                plot_name = origin_indel_hist_plot_obj_root,
                plot_title = 'Indel length at origin',
                plot_label = plot_label,
                plot_datas = [('Indel length at origin histogram',origin_indel_hist_plot_obj_root + ".txt")]
                )

    origin_dist_cutoff = 100000 # only plot items within 100kb of origin cut

    origin_inversion_hist_plot_obj = None # plot of inversion lengths at origin (where read points in opposite genomic direction as primer)
    origin_inversion_hist_plot_obj_root = root + ".origin_inversion_histogram"
    keys = sorted(origin_inversion_lengths.keys())
    vals = [origin_inversion_lengths[key] for key in keys]
    with open(origin_inversion_hist_plot_obj_root+".txt","w") as summary:
        summary.write('Inversion_distance\tnumReads\n')
        for key in keys:
            summary.write(str(key) + '\t' + str(origin_inversion_lengths[key]) + '\n')
    if len(origin_inversion_lengths) > 0 and not suppress_plots:

        plot_keys = [x for x in keys if abs(x) < origin_dist_cutoff]
        plot_vals = [origin_inversion_lengths[key] for key in plot_keys]
        fig = plt.figure(figsize=(12,6))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(plot_keys,plot_vals)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Inversion distance from origin')
        plt.savefig(origin_inversion_hist_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(origin_inversion_hist_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_label = 'Bar plot showing distance from origin cut of reads supporting inversions (up to 100kb). The origin cut is at ' + origin_chr + ':' + str(origin_cut_pos) +' and the origin extends ' + origin_direction + ' from the cut to the primer. ' + \
            'Reads supporting inversions extend ' + origin_direction + ' after genomic alignment. Distances greater than ' + str(origin_dist_cutoff) + 'bp are not shown.'
        if len(keys) == 0:
            plot_label = '(No reads supporting inversions found)'

        origin_inversion_hist_plot_obj = PlotObject(
                plot_name = origin_inversion_hist_plot_obj_root,
                plot_title = 'Inversion distance relative to origin cut',
                plot_label = plot_label,
                plot_datas = [('Inversion distance relative to origin cut (all values)',origin_inversion_hist_plot_obj_root + ".txt")]
                )

    origin_deletion_hist_plot_obj = None # plot of deletion lengths at origin
    origin_deletion_hist_plot_obj_root = root + ".origin_deletion_histogram"
    keys = sorted(origin_deletion_lengths.keys())
    vals = [origin_deletion_lengths[key] for key in keys]
    with open(origin_deletion_hist_plot_obj_root+".txt","w") as summary:
        summary.write('Deletion_distance\tnumReads\n')
        for key in keys:
            summary.write(str(key) + '\t' + str(origin_deletion_lengths[key]) + '\n')

    if len(origin_deletion_lengths) > 0 and not suppress_plots:
        plot_keys = [x for x in keys if abs(x) < origin_dist_cutoff]
        plot_vals = [origin_deletion_lengths[key] for key in plot_keys]
        fig = plt.figure(figsize=(12,6))
        ax = plt.subplot(111)
        if len(plot_vals) > 0:
            ax.bar(plot_keys,plot_vals)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Deletion distance from origin')
        plt.savefig(origin_deletion_hist_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(origin_deletion_hist_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_label = 'Bar plot showing distance from origin cut of reads supporting deletions (up to 100kb). The origin cut is at ' + origin_chr + ':' + str(origin_cut_pos) +' and the origin extends ' + origin_direction + ' from the cut to the primer. ' + \
            'Reads supporting deletions extend ' + origin_direction_opposite + ' after genomic alignment. Distances greater than ' + str(origin_dist_cutoff) + 'bp are not shown.'
        if len(keys) == 0:
            plot_label = '(No reads supporting deletions found)'

        origin_deletion_hist_plot_obj = PlotObject(
                plot_name = origin_deletion_hist_plot_obj_root,
                plot_title = 'Deletion distance relative to origin cut',
                plot_label = plot_label,
                plot_datas = [('Deletion distance relative to origin cut (all values)',origin_deletion_hist_plot_obj_root + ".txt")]
                )

    # plot of coverage of 100bp after cut
    origin_indel_depth_plot_obj_root = root + ".origin_depth_counts"
    origin_depth_counts_100bp = origin_depth_counts_100bp_total - origin_depth_counts_100bp
    origin_depth_counts_100bp_primer = origin_depth_counts_100bp_total - origin_depth_counts_100bp_primer

    xs_primer = list(np.arange(-1*(len(origin_depth_counts_100bp_primer)-1),0,1))
    ys_primer = list(np.flip(origin_depth_counts_100bp_primer[1:]))
    xs_post = list(np.arange(1,len(origin_depth_counts_100bp)))
    ys_post = list(origin_depth_counts_100bp[1:])

    xs = xs_primer + xs_post
    ys = ys_primer + ys_post

    with open(origin_indel_depth_plot_obj_root+"_100bp.txt","w") as summary:
        summary.write('bpFromOriginCut\treadDepth\n')
        for (xval, yval) in zip(xs,ys):
            summary.write(str(xval)+"\t"+str(yval)+"\n")

    origin_depth_counts_long_plot_obj_root = root + ".origin_depth_counts_"+long_distance_str
    origin_depth_counts_long = origin_depth_counts_long_total - origin_depth_counts_long
    origin_depth_counts_long_primer = origin_depth_counts_long_total - origin_depth_counts_long_primer

    xs_long_primer = list(np.arange(-1*(len(origin_depth_counts_long_primer)-1),0,1))
    ys_long_primer = list(np.flip(origin_depth_counts_long_primer[1:]))
    xs_long_post = list(np.arange(1,len(origin_depth_counts_long)))
    ys_long_post = list(origin_depth_counts_long[1:])

    xs_long = xs_long_primer + xs_long_post
    ys_long = ys_long_primer + ys_long_post

    with open(origin_indel_depth_plot_obj_root+"_"+long_distance_str+".txt","w") as summary:
        summary.write('bpFromOriginCut\treadDepth\n')
        for (xval, yval) in zip(xs_long,ys_long):
            summary.write(str(xval)+"\t"+str(yval)+"\n")

    if suppress_plots:
        origin_indel_depth_plot_obj = None
    else:
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(211)
        if len(origin_depth_counts_100bp) > 0:
            ax.bar(xs_primer+xs_post,ys_primer+ys_post)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Read depth (Number of reads)')
        ax.set_title('Read depth of 100bp after origin')
        line1 = ax.axvline(xs_primer[min_primer_length]+0.5,color='k',ls='dotted')
        line2 = ax.axvline(0,color='r',ls='dotted')
        ax.legend([line1,line2],['min primer length','origin location'],loc='lower right',bbox_to_anchor=(1,0))

        # plot of coverage of {long_distance} after cut
        ax = plt.subplot(212)
        all_long_xs = xs_long_primer+xs_long_post
        all_long_ys = ys_long_primer+ys_long_post
        inds_to_keep = list(range(0,len(all_long_xs),int(len(all_long_xs)/100)))
        subset_long_xs = [all_long_xs[x] for x in inds_to_keep]
        subset_long_ys = [all_long_ys[x] for x in inds_to_keep]
        if len(origin_depth_counts_long) > 0:
            ax.bar(subset_long_xs,subset_long_ys)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Read depth (Number of reads)')
        ax.set_title('Read depth of '+long_distance_str+' after origin')
        line1 = ax.axvline(xs_long_primer[min_primer_length]+0.5,color='k',ls='dotted')
        line2 = ax.axvline(0,color='r',ls='dotted')
        ax.legend([line1,line2],['min primer length','origin location'],loc='lower right',bbox_to_anchor=(1,0))

        plt.savefig(origin_indel_depth_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(origin_indel_depth_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        primer_info = 'The vertical red dotted line shows the end of the primer sequence. '
        if origin_chr is not None:
            primer_info = 'Origin cut is at ' + origin_chr + ':' + str(origin_cut_pos) + ' and the origin extends ' + origin_direction + ' from the cut to the primer. The vertical red dotted line shows the origin location. '
        plot_label = 'Read depth surrounding origin. The top plot shows reads with deletions smaller than 100bp (N='+str(origin_depth_counts_100bp_total)+'). '+ \
            'The bottom plot shows reads with deletions smaller than '+long_distance_str+' (N='+str(origin_depth_counts_long_total)+'). ' + primer_info + \
            'The vertical black dotted line shows the minimum requred primer length (' + str(min_primer_length)+'bp) set by the parameter --min_primer_length.'
        if sum(ys_long) < 0:
            plot_label = '(No reads found at origin)'

        origin_indel_depth_plot_obj = PlotObject(
                plot_name = origin_indel_depth_plot_obj_root,
                plot_title = 'Read depth around origin',
                plot_label = plot_label,
                plot_datas = [('Depth around origin (100bp)',origin_indel_depth_plot_obj_root + "_100bp.txt"),
                            ('Depth around origin ('+long_distance_str+')',origin_indel_depth_plot_obj_root + "_"+long_distance_str+".txt")]
                )

    #make r1/r2/support plots
    r1_r2_support_plot_obj = None # plot of how many reads for which r1 and r2 support each other
    r1_r2_support_dist_plot_obj = None # plot of how far apart the r1/r2 supporting reads were aligned from each other
    r1_r2_support_labels = [r1_r2_support_str_supported,r1_r2_support_str_not_supported_orientation,r1_r2_support_str_not_supported_distance,r1_r2_support_str_not_supported_diff_chrs,r1_r2_support_str_not_supported_not_aln]

    #check to make sure labels match
    for k in r1_r2_support_status_counts.keys():
        if k not in r1_r2_support_labels:
            raise Exception('R1/R2 support label %s not found in %s'%(k,str(r1_r2_support_labels)))
    values = [r1_r2_support_status_counts[x] for x in r1_r2_support_labels]
    r1_r2_support_plot_obj_root = root + ".r1_r2_support"
    with open(r1_r2_support_plot_obj_root+".txt",'w') as summary:
        summary.write("\t".join(r1_r2_support_labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")

    if total_reads_processed > 0 and not suppress_plots:
        fig = plt.figure(figsize=(8,8))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(r1_r2_support_labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(r1_r2_support_labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('R1/R2 support status')

        discard_str = '(Reads without R1/R2 support are currently excluded from the final analyses. Set the --suppress_r2_support_filter flag to include these reads.)'
        if suppress_r2_support_filter:
            discard_str = '(Reads without R1/R2 support are currently included in the final analyses because the --suppress_r2_support_filter flag is set.)'
        plt.savefig(r1_r2_support_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(r1_r2_support_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(r1_r2_support_labels,values)])
        r1_r2_support_plot_obj = PlotObject(
                plot_name = r1_r2_support_plot_obj_root,
                plot_title = 'R1/R2 Support Classification',
                plot_label = 'R1/R2 support classification<br>'+plot_count_str+'<br>'+discard_str,
                plot_datas = [
                    ('R1/R2 support classifications',r1_r2_support_plot_obj_root + ".txt")
                    ]
                )

    #plot of distance between R1 and R2 for R1/R2 supportive read pairs
    r1_r2_support_dist_plot_obj_root = root + ".r1_r2_distances"
    keys = sorted(r1_r2_support_distances.keys())
    vals = [r1_r2_support_distances[key] for key in keys]
    with open(r1_r2_support_dist_plot_obj_root+"_supported.txt","w") as summary:
        summary.write('R1_R2_distance\tnumReads\n')
        for key in keys:
            summary.write(str(key) + '\t' + str(r1_r2_support_distances[key]) + '\n')

    #plot of distance between R1 and R2 for R1/R2 NOT supportive read pairs
    keys_2 = sorted(r1_r2_not_support_distances.keys())
    vals_2 = [r1_r2_not_support_distances[key] for key in keys_2]
    with open(r1_r2_support_dist_plot_obj_root+"_not_supported.txt","w") as summary:
        summary.write('R1_R2_distance\tnumReads\n')
        for key in keys_2:
            summary.write(str(key) + '\t' + str(r1_r2_not_support_distances[key]) + '\n')

    if total_reads_processed > 0 and not suppress_plots:
        #plot of distance between R1 and R2 for R1/R2 supportive read pairs
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(211)
        if len(vals) > 0:
            ax.bar(keys,vals)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Distance between R1/R2 for read pairs in which R2 supports the R1 alignment')

        #plot of distance between R1 and R2 for R1/R2 NOT supportive read pairs
        ax2 = plt.subplot(212)
        if len(vals) > 0:
            ax2.bar(keys,vals)
            ax2.set_ymargin(0.05)
        else:
            ax2.bar(0,0)
        ax2.set_ylabel('Number of Reads')
        ax2.set_title('Distance between R1/R2 for read pairs in which R2 does not support the R1 alignment (but R1 and R2 are on the same chromosome)')

        plt.savefig(r1_r2_support_dist_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(r1_r2_support_dist_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_label = '(top) Bar plot showing distance between R1 and R2 reads for reads with R2 supporting the R1 alignment.<br>'
        plot_label += '(bottom) Bar plot showing distance between R1 and R2 reads for reads with R2 NOT supporting the R1 alignment due to incorrect read orientation or exceeding maximum alignment distance (max distance for support is ' + str(r1_r2_support_max_distance) + 'bp)'
        if len(keys) == 0 and len(keys_2) == 0:
            plot_label = '(No R1/R2 reads)'

        r1_r2_support_dist_plot_obj = PlotObject(
                plot_name = r1_r2_support_dist_plot_obj_root,
                plot_title = 'Distance between R1/R2 Reads',
                plot_label = plot_label,
                plot_datas = [('R1/R2 supported distance histogram',r1_r2_support_dist_plot_obj_root + "_supported.txt"),
                    ('R1/R2 not supported distance histogram',r1_r2_support_dist_plot_obj_root + "_not_supported.txt")

                    ]
                )

    discarded_reads_plot_obj_root = root+".discarded_reads"
    labels = ['Unaligned','Duplicate','Poor alignment','Not supported by R2']
    values = [discarded_reads_unaligned,discarded_reads_duplicates,discarded_reads_poor_alignment,discarded_reads_not_supported_by_R2]
    with open(discarded_reads_plot_obj_root+".txt",'w') as summary:
        summary.write("\t".join(labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
    if suppress_plots:
        discarded_reads_plot_obj = None
    else:
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        plot_name = discarded_reads_plot_obj_root
        plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')
        plot_count_str = "<br>".join(["%s N=%d"%x for x in zip(labels,values)])

        discarded_count = discarded_reads_unaligned+discarded_reads_duplicates+discarded_reads_poor_alignment+discarded_reads_not_supported_by_R2

        discard_r1r2_str = '(Reads without R1/R2 support are currently excluded from the final analyses. Set the --suppress_r2_support_filter flag to include these reads.)'
        if suppress_r2_support_filter:
            discard_r1r2_str = '(Reads without R1/R2 support are currently included in the final analyses because the --suppress_r2_support_filter flag is set.)'
        discard_dup_str = '(Reads that may be duplicates based on UMI/alignment position are currently excluded from the final analyses. Set the --suppress_dedup_on_aln_pos_and_UMI_filter flag to include these reads.)'
        if suppress_dedup_on_aln_pos_and_UMI_filter:
            discard_dup_str = '(Duplicate reads based on UMI/alignment position are currently included in the final analyses because the --suppress_dedup_on_aln_pos_and_UMI_filter flag is set.)'

        discard_poor_aln_str = '(Reads with poor alignment are currently excluded from the final analyses. Set the --suppress_poor_alignment_filter flag to include these reads.)'
        if suppress_poor_alignment_filter:
            discard_poor_aln_str = '(Reads with poor alignment are currently included in the final analyses because the --suppress_poor_alignment_filter flag is set.)'


        discarded_reads_plot_obj = PlotObject(
                plot_name = discarded_reads_plot_obj_root,
                plot_title = 'Discarded read count summary',
                plot_label = str(discarded_count) + '/'+str(total_r1_processed) + ' reads were discarded.<br>'+str(final_total_count) + ' R1 reads were used for the final analysis. <br>' + plot_count_str + '<br>' +
                    discard_dup_str + '<br>' + discard_poor_aln_str + '<br>' + discard_r1r2_str + '<br>Poor alignment: Reads are aligned such that the front AND back of the read have greater than ' +
                    str(arm_max_clipped_bases) +'bp clipped, or the start or the end of the read have greater than ' + str(arm_min_matched_start_bases) + 'bp matching the reference. These cutoffs can be adjusted using the --arm_max_clipped_bases and --arm_min_matched_start_bases parameters.',
                plot_datas = [
                    ('Discarded reads summary',discarded_reads_plot_obj_root + ".txt")
                    ]
                )
    if write_discarded_read_info:
        discarded_reads_plot_obj.datas.append(('Discarded reads',discarded_read_file))

    with open(info_file,'w') as fout:
        fout.write("\t".join(['total_reads_processed','total_r1_processed','discarded_reads_unaligned','discarded_reads_duplicates','discarded_reads_not_supported_by_R2','discarded_reads_poor_alignment','final_read_count','discovered_cut_point_count','final_cut_point_count'])+"\n")
        fout.write("\t".join(str(x) for x in [total_reads_processed,total_r1_processed,discarded_reads_unaligned,discarded_reads_duplicates,discarded_reads_not_supported_by_R2,discarded_reads_poor_alignment,final_total_count,found_cut_point_total,final_cut_point_total])+"\n")
        fout.write(classification_read_counts_str)
        fout.write(classification_indel_read_counts_str)
        classification_json_str = json.dumps(cut_classification_lookup,separators=(',',':'))
        fout.write(classification_json_str+"\n")

        chr_aln_plot_obj_str = "None"
        if chr_aln_plot_obj is not None:
            chr_aln_plot_obj_str = chr_aln_plot_obj.to_json()
        fout.write(chr_aln_plot_obj_str+"\n")

        tlen_plot_obj_str = "None"
        if tlen_plot_obj is not None:
            tlen_plot_obj_str = tlen_plot_obj.to_json()
        fout.write(tlen_plot_obj_str+"\n")

        deduplication_plot_obj_str = "None"
        if deduplication_plot_obj is not None:
            deduplication_plot_obj_str = deduplication_plot_obj.to_json()
        fout.write(deduplication_plot_obj_str+"\n")

        tx_order_plot_obj_str = "None"
        if tx_order_plot_obj is not None:
            tx_order_plot_obj_str = tx_order_plot_obj.to_json()
        fout.write(tx_order_plot_obj_str+"\n")

        tx_count_plot_obj_str = "None"
        if tx_count_plot_obj is not None:
            tx_count_plot_obj_str = tx_count_plot_obj.to_json()
        fout.write(tx_count_plot_obj_str+"\n")

        classification_plot_obj_str = "None"
        if classification_plot_obj is not None:
            classification_plot_obj_str = classification_plot_obj.to_json()
        fout.write(classification_plot_obj_str+"\n")

        classification_indel_plot_obj_str = "None"
        if classification_indel_plot_obj is not None:
            classification_indel_plot_obj_str = classification_indel_plot_obj.to_json()
        fout.write(classification_indel_plot_obj_str+"\n")

        origin_indel_hist_plot_obj_str = "None"
        if origin_indel_hist_plot_obj is not None:
            origin_indel_hist_plot_obj_str = origin_indel_hist_plot_obj.to_json()
        fout.write(origin_indel_hist_plot_obj_str+"\n")

        origin_inversion_hist_plot_obj_str = "None"
        if origin_inversion_hist_plot_obj is not None:
            origin_inversion_hist_plot_obj_str = origin_inversion_hist_plot_obj.to_json()
        fout.write(origin_inversion_hist_plot_obj_str+"\n")

        origin_deletion_hist_plot_obj_str = "None"
        if origin_deletion_hist_plot_obj is not None:
            origin_deletion_hist_plot_obj_str = origin_deletion_hist_plot_obj.to_json()
        fout.write(origin_deletion_hist_plot_obj_str+"\n")

        origin_indel_depth_plot_obj_str = "None"
        if origin_indel_depth_plot_obj is not None:
            origin_indel_depth_plot_obj_str = origin_indel_depth_plot_obj.to_json()
        fout.write(origin_indel_depth_plot_obj_str+"\n")

        r1_r2_support_plot_obj_str = "None"
        if r1_r2_support_plot_obj is not None:
            r1_r2_support_plot_obj_str = r1_r2_support_plot_obj.to_json()
        fout.write(r1_r2_support_plot_obj_str+"\n")

        r1_r2_support_dist_plot_obj_str = "None"
        if r1_r2_support_dist_plot_obj is not None:
            r1_r2_support_dist_plot_obj_str = r1_r2_support_dist_plot_obj.to_json()
        fout.write(r1_r2_support_dist_plot_obj_str+"\n")

        discarded_reads_plot_obj_str = "None"
        if discarded_reads_plot_obj is not None:
            discarded_reads_plot_obj_str = discarded_reads_plot_obj.to_json()
        fout.write(discarded_reads_plot_obj_str+"\n")

    if not keep_intermediate:
        if os.path.exists(final_file_tmp):
            logger.debug('Removing tmp assignments file: ' + final_file_tmp)
            os.remove(final_file_tmp)
        if os.path.exists(genome_mapped_bam):
            logger.debug('Removing genome alignment bam file: ' + genome_mapped_bam)
            os.remove(genome_mapped_bam)
        if os.path.exists(genome_mapped_bam + ".bai"):
            logger.debug('Removing genome alignment bam index file: ' + genome_mapped_bam + ".bai")
            os.remove(genome_mapped_bam + ".bai")

    discarded_read_counts = [('Unaligned',discarded_reads_unaligned),('Duplicate read',discarded_reads_duplicates),('Not supported by R2',discarded_reads_not_supported_by_R2),('Bad alignment',discarded_reads_poor_alignment)]
    return (final_file,final_read_ids_for_crispresso,final_cut_counts,cut_classification_lookup, final_total_count, discarded_read_counts, classification_read_counts, classification_indel_read_counts,
        chr_aln_plot_obj,tlen_plot_obj,deduplication_plot_obj,tx_order_plot_obj,tx_count_plot_obj,tx_circos_plot_obj,classification_plot_obj,classification_indel_plot_obj,
        origin_indel_hist_plot_obj,origin_inversion_hist_plot_obj,origin_deletion_hist_plot_obj,origin_indel_depth_plot_obj,
        r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,discarded_reads_plot_obj)

def collapse_cut_points(novel_cut_merge_distance, known_cut_merge_distance, origin_cut_merge_distance,
        cut_points_by_chr, cut_annotations, cut_point_homology_info):
    """
    Collapse observed cut points to final cut points based on merging distances

    Args:
        novel_cut_merge_distance (int): distance to merge novel cut points
        known_cut_merge_distance (int): distance to merge known cut points
        origin_cut_merge_distance (int): distance to merge origin cut points
        cut_points_by_chr (dict): dictionary of known cut points by chromosome
        cut_annotations (dict): dict of cut_chr:cut_site->annotation for description of cut [type, origin_status, guide_direction, guide_sequence] for known/given cut sites
        cut_point_homology_info (dict): dictionary of homology info by cut point

    Returns:
        origin_chr (str): name of origin chromosome
        origin_cut (int): position of origin cut
        origin_direction (str): direction of primer from origin cut
        final_cut_point_lookup (dict): dictionary of final cut point lookups old_cut_point > new_cut_point
        final_cut_points_by_chr (dict): dictionary of final cut points chr->cut point
        found_cut_point_total (int): number of final cut points

    """

    logger = logging.getLogger('CRISPRlungo')
    logger.info('Merging cut points')
    origin_chr = None
    origin_cut_pos = None
    origin_direction = None
    known_cut_points_by_chr = {} #cut sites are input by user (or found by casoffinder)
    final_cut_points_by_chr = defaultdict(lambda: defaultdict(int)) #dict of final cut points by chromosome, contains count of reads with that cut point
    final_cut_point_lookup = {} #dict from old(fuzzy/imprecise) position to new
    for cut_site in cut_annotations:
        cut_chr,cut_pos_str = cut_site.split(":")
        cut_pos = int(cut_pos_str)
        if cut_site in cut_annotations:
            origin_status = cut_annotations[cut_site][1]
            if origin_status == 'Origin:right':
                origin_chr = cut_chr
                origin_cut_pos = cut_pos
                origin_direction = 'right'
            elif origin_status == 'Origin:left':
                origin_chr = cut_chr
                origin_cut_pos = cut_pos
                origin_direction = 'left'

        if cut_chr not in known_cut_points_by_chr:
            known_cut_points_by_chr[cut_chr] = []
        known_cut_points_by_chr[cut_chr].append(cut_pos)
        #add known cut points to final cut point lookup as well
        if cut_chr not in final_cut_points_by_chr:
            final_cut_points_by_chr[cut_chr] = defaultdict(int)
        final_cut_points_by_chr[cut_chr][cut_pos] = 0

    #create dict of cut_point > closest cut point that is origin, known or with homology
    closest_origin_known_homology_lookup = {}

    #first, assign applicable points to origin - sites known or homology will overwrite this
    origin_range_min = -1 # for non-genomic primers
    origin_range_max = -1
    if origin_cut_pos is not None: #for genomic primers
        origin_range_min = origin_cut_pos - origin_cut_merge_distance
        origin_range_max = origin_cut_pos + origin_cut_merge_distance

    for cut_point in cut_points_by_chr[origin_chr]:
        if cut_point >= origin_range_min and cut_point <= origin_range_max:
            closest_origin_known_homology_lookup[origin_chr+":"+str(cut_point)] = origin_chr+":"+str(origin_cut_pos)

    #next create dict with known sites and sites with homology
    cut_points_known_or_homology_by_chr = defaultdict(list)
    #add known sites
    for cut_chr in known_cut_points_by_chr:
        cut_points_known_or_homology_by_chr[cut_chr] = known_cut_points_by_chr[cut_chr][:]
    #split cut points with homology by chr to allow sorting
    for cut_point in cut_point_homology_info:
        cut_chr,cut_pos = cut_point.split(":")
        if int(cut_pos) not in cut_points_known_or_homology_by_chr[cut_chr]:
            cut_points_known_or_homology_by_chr[cut_chr].append(int(cut_pos))


    #assign each cut point to the closest origin, known site or site with homology if possible
    for cut_chr in cut_points_by_chr:
        these_cut_points = sorted(cut_points_by_chr[cut_chr])

        these_cut_points_known_or_homology = sorted(cut_points_known_or_homology_by_chr[cut_chr])
        #if no cut points on this chr
        if len(these_cut_points_known_or_homology) == 0:
            for cut_point in these_cut_points:
                closest_origin_known_homology_lookup[cut_chr+":"+str(cut_point)] = None
        else:
            curr_ind = 0
            for cut_point in these_cut_points:
                dist_to_curr = abs(cut_point - these_cut_points_known_or_homology[curr_ind])
                #if curr_ind is at last item, choose it
                if curr_ind + 1 >= len(these_cut_points_known_or_homology):
                    if dist_to_curr <= known_cut_merge_distance:
                        closest_origin_known_homology_lookup[cut_chr+":"+str(cut_point)] = cut_chr+":"+str(these_cut_points_known_or_homology[curr_ind])
                    else:
                        closest_origin_known_homology_lookup[cut_chr+":"+str(cut_point)] = None
                #otherwise choose closest of curr or next
                else:
                    dist_to_next = abs(cut_point - these_cut_points_known_or_homology[curr_ind+1])
                    if dist_to_curr <= dist_to_next:
                        if dist_to_curr <= known_cut_merge_distance:
                            closest_origin_known_homology_lookup[cut_chr+":"+str(cut_point)] = cut_chr+":"+str(these_cut_points_known_or_homology[curr_ind])
                        else:
                            closest_origin_known_homology_lookup[cut_chr+":"+str(cut_point)] = None
                    else:
                        curr_ind += 1
                        if dist_to_next <= known_cut_merge_distance:
                            closest_origin_known_homology_lookup[cut_chr+":"+str(cut_point)] = cut_chr+":"+str(these_cut_points_known_or_homology[curr_ind])
                        else:
                            closest_origin_known_homology_lookup[cut_chr+":"+str(cut_point)] = None

    #cluster remaining novel points
    found_cut_point_total = 0
    closest_cut_point_lookup = {}
    for cut_chr in cut_points_by_chr:
        found_cut_point_total += len(cut_points_by_chr[cut_chr])
        these_cut_points = sorted(cut_points_by_chr[cut_chr])
        for cut_point in these_cut_points:
            this_cut_str = cut_chr+":"+str(cut_point)
            if closest_origin_known_homology_lookup[this_cut_str] is not None:
                these_cut_points.remove(cut_point)
                closest_cut_point_lookup[this_cut_str] = closest_origin_known_homology_lookup[this_cut_str]

        if len(these_cut_points) > 0:
            last_seen_point = these_cut_points[0]
            curr_points = [last_seen_point]
            for i in range(1,len(these_cut_points)+1):
                #if next point is beyond novel_cut_merge_dist, merge and reset!
                #do weighted sum
                this_sum = 0
                this_tot = 0
                for curr_point in curr_points:
                    this_sum += curr_point*cut_points_by_chr[cut_chr][curr_point]
                    this_tot += cut_points_by_chr[cut_chr][curr_point]
                this_weighted_mean = 0
                if this_tot > 0:
                    this_weighted_mean = this_sum/float(this_tot)
                #this_mean = int(sum(curr_points)/float(len(curr_points)))
                if i == len(these_cut_points) or abs(these_cut_points[i] - this_weighted_mean) > novel_cut_merge_distance:
                    this_pos = int(this_weighted_mean)
                    #add this assignment to the lookup
                    for curr_point in curr_points:
                        closest_cut_point_lookup[cut_chr+":"+str(curr_point)] = cut_chr+":"+str(this_pos)

                    #reset current points
                    curr_points = []

                if i < len(these_cut_points):
                    curr_points.append(these_cut_points[i])

    #make final assignments - either to the next homologous site (if available), to the origin (if nearby), or to the next closest site
    all_cut_points_count = 0
    all_cut_points_assigned_to_homology_count = 0 #how many cut sites were assigned to their nearest homologous/known site
    all_cut_points_assigned_to_origin_count = 0 #how many cut sites were assigned to the origin site
    all_cut_points_collapsed_count = 0 #how many cut sites were collapsed (not assigned to origin or homologous sites)
    for cut_chr in cut_points_by_chr:
        for cut_point in cut_points_by_chr[cut_chr]:
            all_cut_points_count += 1
            cut_key = cut_chr+":"+str(cut_point)
            if closest_origin_known_homology_lookup[cut_key] is not None:
                cut_assignment = closest_origin_known_homology_lookup[cut_key]
                all_cut_points_assigned_to_homology_count += 1
            elif cut_chr == origin_chr and cut_point >= origin_range_min and cut_point <= origin_range_max:
                cut_assignment = origin_chr+":"+str(origin_cut_pos)
                all_cut_points_assigned_to_origin_count += 1
            else:
                cut_assignment = closest_cut_point_lookup[cut_key]
                all_cut_points_collapsed_count += 1

            cut_assignment_pos = int(cut_assignment.split(":")[1])
            final_cut_point_lookup[cut_key] = cut_assignment
            final_cut_points_by_chr[cut_chr][cut_assignment_pos] += cut_points_by_chr[cut_chr][cut_point]

    logger.info('Processed ' + str(all_cut_points_count) + ' observed cut points')
    logger.info(str(all_cut_points_assigned_to_homology_count) + ' cut points were collapsed to known/given/homologous cut sites')
    logger.info(str(all_cut_points_assigned_to_origin_count) + ' cut points were collapsed to the origin site')
    logger.info(str(all_cut_points_collapsed_count) + ' cut points were collapsed to nearby cut sites')

    return origin_chr, origin_cut_pos, origin_direction, final_cut_point_lookup, final_cut_points_by_chr, found_cut_point_total

def plot_tx_circos(genome,final_cut_counts,origin_chr,origin_cut_pos,figure_root,max_translocations_to_plot=100,include_all_chrs=False,plot_black_white=True):
    """Plot a circos plot showing translocations.

    Args:
        genome (str): Genome file location (specifically, this function uses the .fai index file to find chrom lengths)
        final_cut_counts: dict of cutID=>count how many times each cut was seen -- when we iterate through the reads (fastq) in the next step, we check to see that the cut was seen above a threshold before printing those reads. The count is stored in this dict.
        origin_chr (str): name of chromosome where origin is (or None if primer is not genomic)
        origin_cut_pos (int): genomic location of origin cut
        figure_root (str): root of files to plot - '.pdf' and '.png' will be added
        max_translocations_to_plot (int): maximum number of translocations to show in circos plot
        include_all_chrs (bool, optional): Whether to plot all chromosomes. Defaults to False and only plots chromosomes that don't include an '_' character.
        plot_black_white (bool, optional): If true, arcs will be plotted grayscale. Otherwise, they will be colored the color of the target chromosome
    """
    genome_index_file = genome + ".fai" #we've checked that this file exists before
    all_chrs = []
    chr_lens = {}
    with open(genome_index_file,'r') as fin:
        for line in fin:
            line_els = line.split("\t")
            all_chrs.append(line_els[0])
            chr_lens[line_els[0]] = int(line_els[1])

    tx_data = defaultdict(int)
    tx_chrs = {}
    max_tx_count = 0
    min_tx_count = None
    tx_keys = final_cut_counts.keys()
    for tx_key in tx_keys:
        cut_els = tx_key.split(":")
        tx_chrs[cut_els[0]] = 1
        this_count = final_cut_counts[tx_key]
        tx_data[cut_els[0]+":"+cut_els[1]] += this_count #could have multiple counts for a single site (both strands)
        this_total_count = tx_data[cut_els[0]+":"+cut_els[1]]

        if this_total_count > max_tx_count:
            max_tx_count = this_total_count
        if min_tx_count is None:
            min_tx_count = this_total_count
        if this_total_count < min_tx_count:
            min_tx_count = this_total_count

    circle = pc.Gcircle(figsize=(8,8))

    chrom_raxis_range=(600,700)
    chrom_label_raxis_range=(700,710)
    chrom_counts_raxis_range=(500,600)
    arc_raxis = 490

    arcdata_dict = defaultdict(dict)

    min_chr_size = chr_lens[all_chrs[0]]
    max_chr_size = min_chr_size
    good_chrs = []
    for chrom in all_chrs:
        is_good_chr = True
        if '_' in chrom and not include_all_chrs:
            is_good_chr = False
        if chrom in tx_chrs:
            is_good_chr = True

        if is_good_chr:
            good_chrs.append(chrom)
            if chr_lens[chrom] < min_chr_size:
                min_chr_size = chr_lens[chrom]
            if chr_lens[chrom] > max_chr_size:
                max_chr_size = chr_lens[chrom]

    #add origin chr
    if origin_chr is None:
        chr_lens['primer'] = int(float(min_chr_size/2))
        good_chrs.append('primer')
        origin_multiplier = 1
        origin_chr = 'primer'
        origin_cut_pos = 40
    else:
        origin_multiplier = int((max_chr_size*4)/chr_lens[origin_chr])
        chr_lens[origin_chr] *= origin_multiplier

    for chrom in good_chrs:
        arc = pc.Garc(arc_id=chrom, size=chr_lens[chrom], interspace=2, raxis_range=chrom_raxis_range, labelposition=80, label_visible=True)
        circle.add_garc(arc)

        #for scatter plot
        if chrom not in arcdata_dict:
            arcdata_dict[chrom]["positions"] = []
            arcdata_dict[chrom]["values"] = []


    circle.set_garcs()

    for arc_id in circle.garc_dict:
        this_tick_interval = 20000000
        if arc_id == origin_chr:
            this_tick_interval *= origin_multiplier

        circle.tickplot(arc_id, raxis_range=chrom_label_raxis_range, tickinterval=this_tick_interval, ticklabels=None)




    log_max_tx_count = math.log(max_tx_count)
    if min_tx_count == max_tx_count:
        min_tx_count = 0.001
    log_min_tx_count = math.log(min_tx_count)
    #scatter plot

    sorted_tx_keys = sorted(tx_data.items(), key=lambda item: item[1], reverse=True)
    sorted_tx_keys_to_plot = sorted_tx_keys[0:min(len(sorted_tx_keys), max_translocations_to_plot)]
    values_all   = []
    for tx_key,tx_count in sorted_tx_keys_to_plot:
        tx_chr,tx_loc = tx_key.split(":")
        tx_loc = int(tx_loc)
        if tx_chr == origin_chr:
            tx_loc = tx_loc * origin_multiplier
        this_log = math.log(tx_data[tx_key])
        this_transform = (log_max_tx_count - this_log) + log_min_tx_count
        values_all.append(this_transform)
        arcdata_dict[tx_chr]["positions"].append(tx_loc)
        arcdata_dict[tx_chr]["values"].append(this_transform)

    for key in arcdata_dict:
        circle.scatterplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
                        rlim=[log_min_tx_count-0.1*abs(log_min_tx_count), log_max_tx_count+0.1*abs(log_max_tx_count)], raxis_range=chrom_counts_raxis_range, facecolor='k', spine=True)


    #add txs
    if plot_black_white:
        cmap = plt.cm.Greys
    arcdata_dict = defaultdict(dict)
    # reverse so the highest counts get plotted last/on top
    for tx_key,tx_count in sorted_tx_keys_to_plot[::-1]:
        tx_chr,tx_loc = tx_key.split(":")
        tx_loc = int(tx_loc)
        name1  = origin_chr
        start1 = origin_cut_pos * origin_multiplier
        end1   = start1 + 20
        name2  = tx_chr
        start2 = tx_loc
        if name2 == origin_chr:
            start2 = start2 * origin_multiplier
        end2   = start2 + 20
        this_count = tx_data[tx_key]
        source = (name1, start1, end1, arc_raxis)
        destination = (name2, start2, end2, arc_raxis)
        color_pct = this_count/max_tx_count
        if color_pct < 0.1:
            color_pct = 0.1
        if plot_black_white:
            circle.chord_plot(source, destination,facecolor=cmap(color_pct),linewidth=color_pct,edgecolor=cmap(color_pct))
        else:
            circle.chord_plot(source, destination, edgecolor=circle.garc_dict[name2].facecolor,linewidth=color_pct*4)

    circle.figure.savefig(figure_root + ".pdf")
    circle.figure.savefig(figure_root + ".png")

def read_command_output(command):
    """
    Runs a shell command and returns an iter to read the output

    Args:
        command: shell command to run

    Returns:
        iter to read the output
    """

    p = subprocess.Popen(command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,shell=True,
#            encoding='utf-8',universal_newlines=True)
            universal_newlines=True,
            bufsize=-1) #bufsize system default
    return iter(p.stdout.readline, b'')


def prep_crispresso2(root,input_fastq_file,read_ids_for_crispresso,cut_counts_for_crispresso,cut_classification_lookup,cut_annotations,
        av_read_length,origin_seq,cleavage_offset,genome,genome_len_file,crispresso_min_count,crispresso_min_aln_score,crispresso_quant_window_size,
        run_crispresso_on_novel_sites,samtools_command,crispresso_command,n_processes=1,use_counts_from_both_sides_of_cuts=True):
    """
    Prepares reads for analysis by crispresso2
    Frequently-aligned locations with a min number of reads (crispresso_min_count) are identified
    Reads from these locations are extracted and prepared for analysis by CRISPResso2

    Args:
        root: root for written files
        input_fastq_file: input fastq with read sequences
        read_ids_for_crispresso: dict of readID=>cut assignment
        cut_counts_for_crispresso: dict of cut=>count for # reads assigned to each cut
        cut_classification_lookup: dict of cutID=> type (e.g. Linear, Translocation, etc)
        cut_annotations: dict of cut_chr:cut_site->annotation for description of cut [type, origin_status, guide_direction, guide_sequence] 
            (type=Programmed, Off-target, Known, Casoffinder, etc) (origin_status=Not-origin or Origin:left or Origin:right) (guide_direction='FW' or 'RC' for which direction the guide binds)
        av_read_length: average read length (used to calculate reference region size for CRISPResso)
        origin_seq: common amplified region at primer (to first cut site)
        cleavage_offset: position where cleavage occurs (relative to end of spacer seq -- for Cas9 this is -3)
        genome: path to genome fa file
        genome_len_file: path to tab-sep file with lengths of chrs in genome (ends with .fai)
        crispresso_min_count: min number of reads at site for crispresso2 processing
        crispresso_min_aln_score: minimum score for reads to align to amplicons
        crispresso_quant_window_size: number of bases on each side of the cut site in which to consider editing events
        run_crispresso_on_novel_sites: if false, only predicted sites (as input by the user as on or off-targets) are analyzed using CRISPResso (no novel sites)
        samtools_command: location of samtools to run
        crispresso_command: location of crispresso to run
        n_processes: number of processes to run CRISPResso commands on
        use_counts_from_both_sides_of_cuts: boolean; cuts must have at least cut_counts_for_crispresso reads aligned to them for crispresso analysis. If use_counts_from_both_sides_of_cuts is true, this count is the sum of the number of left and right reads. Otherwise if false, the count is only the left (or right) count. If use_counts_from_both_sides_of_cuts, for sites where the left or right side has a lot of reads, both sides will be analyzed by CRISPResso (if they have reads)

    Returns:
        crispresso_infos: dict containing metadata about each crispresso run
            #dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command
    """
    logger = logging.getLogger('CRISPRlungo')
    logger.info('Preparing reads for analysis by CRISPResso2')

    chrom_lens = {}
    chroms = []
    with open(genome_len_file,'r') as gfile:
        for line in gfile:
            line_els = line.split("\t")
            line_chrom = line_els[0]
            line_len = line_els[1]
            chrom_lens[line_chrom] = line_len
            chroms.append(line_chrom)

    crispresso_infos = [] #information about each crispresso_name

    data_dir = root+'_data'
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    total_read_count = 0
    printed_read_count = 0
    filehandles = {}

    printed_cuts = defaultdict(int)
    cut_names = {}

    #iterate through fastq file
    # if read id is in the read_ids_for_crispresso, print it to that file
    if input_fastq_file.endswith('.gz'):
        f_in = gzip.open(input_fastq_file,'rt')
    else:
        f_in = open(input_fastq_file,'r')
    while (1):
        id_line   = f_in.readline().strip()
        seq_line  = f_in.readline().strip()
        plus_line = f_in.readline()
        qual_line = f_in.readline().strip()

        if not qual_line : break

        total_read_count+=1

        id_line_els = id_line.split("\t")
        read_id = id_line_els[0][1:] #chop off leading @

        if read_id not in read_ids_for_crispresso:
            continue

        this_cut = read_ids_for_crispresso[read_id]

        this_count = cut_counts_for_crispresso[this_cut]
        if use_counts_from_both_sides_of_cuts:
            this_cut_els = this_cut.split(":")
            this_cut_pos = ":".join(this_cut_els[0:2])
            left_count = cut_counts_for_crispresso[this_cut_pos+":left"]
            right_count = cut_counts_for_crispresso[this_cut_pos+":right"]
            this_count = left_count + right_count

        if this_count < crispresso_min_count or cut_counts_for_crispresso[this_cut] == 0:
            continue

        if this_cut not in cut_names:
            cut_name = (cut_classification_lookup[this_cut] + ':' + this_cut).replace(" ","_").replace(":","_")
            cut_names[this_cut] = cut_name
        cut_name = cut_names[this_cut]

        printed_read_count += 1
        reads_file = os.path.join(data_dir,cut_name+".fq")
        if reads_file not in filehandles:
            fh = open(reads_file,'w')
            filehandles[reads_file] = fh
        filehandles[reads_file].write("%s\n%s\n%s%s\n"%(id_line,seq_line.upper(),plus_line,qual_line))

        printed_cuts[this_cut] += 1

    amp_half_length = av_read_length/1.5
    for i,cut in enumerate(printed_cuts):
        (cut_chr,cut_pos,cut_direction) = cut.split(":")
        cut_pos = int(cut_pos)
        cut_name = cut_names[cut]

        cut_guide_annotation = ['Novel','Not-origin','NA']
        cut_str = cut_chr+':'+str(cut_pos)
        if cut_str in cut_annotations:
            cut_guide_annotation = cut_annotations[cut_str]

        if cut_direction == "right":
            cut_amp_start = cut_pos
            cut_amp_end = cut_pos + amp_half_length
            faidx_cmd = '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,cut_chr,cut_amp_start,cut_amp_end)
            logger.debug('Getting sequence on right for cut ' + str(cut) + ': ' + str(faidx_cmd))
            cut_amp_seq = subprocess.check_output(faidx_cmd,shell=True).decode(sys.stdout.encoding).strip()

            # I had attempted to recover the actual guide on each side of the cut, but sometimes it is only 6bp long and not unique within the amplicon sequence
#            if cut_guide_annotation[2] == 'FW':
#                cut_guide_seq = cut_guide_annotation[3][cleavage_offset:]
#            elif cut_guide_annotation[2] == 'RV':
#                cut_guide_seq = cut_guide_annotation[3][0:-cleavage_offset]
#            else:
#                cut_guide_seq = cut_amp_seq[0:10]


        elif cut_direction == "left":
            cut_amp_end = cut_pos
            cut_amp_start = cut_pos - amp_half_length

            faidx_cmd = '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,cut_chr,cut_amp_start,cut_amp_end - 1)
            logger.debug('Getting sequence on left for cut ' + str(cut) + ': ' + str(faidx_cmd))
            cut_amp_seq = subprocess.check_output(faidx_cmd,shell=True).decode(sys.stdout.encoding).strip()
            cut_amp_seq = reverse_complement(cut_amp_seq)

#            if cut_guide_annotation[2] == 'FW':
#                cut_guide_seq = cut_guide_annotation[3][0:-cleavage_offset]
#            elif cut_guide_annotation[2] == 'RV':
#                cut_guide_seq = cut_guide_annotation[3][cleavage_offset:]
#            else:
#                cut_guide_seq = cut_amp_seq[0:10]

        else:
            raise Exception('Got unexpected cut direction: ' + cut_direction + ' from cut ' + cut)


        amp_seq = origin_seq + cut_amp_seq
        post_cut_seq = cut_amp_seq[0:15] #3bp + 3bp is sometimes too short to uniquely identify the location so we'll just use 15bp on each side. 
        guide_seq = origin_seq[-15:] + post_cut_seq
        quant_window_start = len(origin_seq) - crispresso_quant_window_size
        quant_window_coords = "%d-%d"%(quant_window_start, quant_window_start + 2*(crispresso_quant_window_size) - 1)

        reads_file = os.path.join(data_dir,cut_name+".fq")
        processes_str = ""
        if n_processes > 4:
            processes_str = "--n_processes 4"
        output_folder = os.path.join(data_dir,'CRISPResso_on_'+cut_names[cut])
        crispresso_cmd = "%s -o %s -n %s --default_min_aln_score %d -a %s -g %s -wc %s -w %s -r1 %s -gn 'Fragment junction' --fastq_output --quantification_window_coordinates %s %s &> %s.log"%(crispresso_command,data_dir,cut_name,crispresso_min_aln_score,
            amp_seq,guide_seq,-15,crispresso_quant_window_size,reads_file,quant_window_coords,processes_str,reads_file)
        
        run_this_one = 'True'
        if not run_crispresso_on_novel_sites and 'Novel' in cut_classification_lookup[cut]:
            run_this_one = 'False: novel location'

        if not use_counts_from_both_sides_of_cuts and printed_cuts[cut] < crispresso_min_count:
            run_this_one = 'False: too few reads (' + str(printed_cuts[cut]) + ')'
        crispresso_infos.append({
                "cut":cut,
                "name":cut_names[cut],
                "type":cut_classification_lookup[cut],
                "cut_loc": cut_chr + ":" + str(cut_pos),
                "amp_seq": amp_seq,
                "output_folder":output_folder,
                "reads_file": reads_file,
                "printed_read_count": printed_cuts[cut],
                "command": crispresso_cmd,
                "run": run_this_one
                })

    return (crispresso_infos)

def run_and_aggregate_crispresso(root,crispresso_infos,final_assignment_file,n_processes,skip_failed=True,keep_intermediate=False,suppress_plots=False,can_use_previous_analysis=False):
    """
    Runs CRISPResso2 commands and aggregates output

    Args:
        root: root for written files
        crispresso_infos: array of metadata information for CRISPResso
            dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command, run
        final_assignment_file: filename of final assignments for each read id pair
        n_processes: number of processes to run CRISPResso commands on
        skip_failed: if true, failed CRISPResso runs are skipped. Otherwise, if one fails, the program will fail.
        keep_intermediate: whether to keep intermediate files (if False, intermediate files including produced fastqs will be deleted)
        suppress_plots: if true, plotting will be suppressed
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    Returns:
        crispresso_results: dict of run_names and run_sub_htmls for display in report
        crispresso_classification_plot_obj: plot object with classifications includign crispresso analysis

    """
    logger = logging.getLogger('CRISPRlungo')
    crispresso_stats_file = root + ".info"

    crispresso_commands = []
    for crispresso_info in crispresso_infos:
        if crispresso_info['run'] == 'True':
            crispresso_commands.append(crispresso_info['command'])
    expected_completed_command_count = len(crispresso_commands)

    #check to see if runs have been completed previously
    crispresso_run_names = []
    crispresso_sub_htmls = {}
    read_completed_command_count = 0
    if os.path.isfile(crispresso_stats_file) and can_use_previous_analysis:
        classification_plot_obj = None
        with open(crispresso_stats_file,'r') as fin:
            classification_plot_obj_str = fin.readline().rstrip('\n')
            if classification_plot_obj_str != "" and classification_plot_obj_str != "None":
                classification_plot_obj = PlotObject.from_json(classification_plot_obj_str)
            spacer_line = fin.readline()

            head_line = fin.readline()
            for line in fin:
                line_els = line.strip().split()
                crispresso_run_names.append(line_els[0])
                crispresso_sub_htmls[line_els[0]] = line_els[1]
                read_completed_command_count += 1
        if read_completed_command_count == expected_completed_command_count:
            logger.info('Using %d previously-completed CRISPResso runs'%read_completed_command_count)
            crispresso_results = {}
            crispresso_results['run_names'] = crispresso_run_names
            crispresso_results['run_sub_htmls'] = crispresso_sub_htmls
            return(crispresso_results,classification_plot_obj)
        else:
            logger.info('Could not recover previously-completed CRISPResso runs. Rerunning.')


    #otherwise run crispresso runs
    logger.info('Running and analyzing ' + str(len(crispresso_infos)) + ' alignments using CRISPResso2')


    CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_commands,n_processes,'run',skip_failed)

    all_crispresso_read_modified_data = {} #dict storing read id > modifed/unmodified
    crispresso_results = {}
    crispresso_results['run_names'] = []
    crispresso_results['run_sub_htmls'] = {}
    crispresso_command_file = root + '.commands.txt'
    crispresso_info_file = root + '.summary.txt'
    with open(crispresso_info_file,'w') as crispresso_info_fh, open(crispresso_command_file,'w') as crispresso_command_fh:
        crispresso_info_fh.write('name\tcut_annotation\tcut_type\tcut_location\treference\treads_printed\tn_total\treads_aligned\treads_unmod\treads_mod\treads_discarded\t'+\
            'reads_insertion\treads_deletion\treads_substitution\treads_only_insertion\treads_only_deletion\treads_only_substitution\treads_insertion_and_deletion\t'+\
            'reads_insertion_and_substitution\treads_deletion_and_substitution\treads_insertion_and_deletion_and_substitution\tamplicon_sequence\n')
        for crispresso_info in crispresso_infos:
            crispresso_command_fh.write(crispresso_info['command'])
            name = crispresso_info['name']
            cut = crispresso_info['cut']
            cut_type = crispresso_info['type']
            cut_loc = crispresso_info['cut_loc']
            cut_amp_seq = crispresso_info['amp_seq']
            n_printed = crispresso_info['printed_read_count']
            run_status = crispresso_info['run']
            ref_name = "NA"
            n_total = "NA"
            n_aligned = "NA"
            n_unmod = "NA"
            n_mod = "NA"
            n_discarded = "NA"

            n_insertion = "NA"
            n_deletion = "NA"
            n_substitution = "NA"
            n_only_insertion = "NA"
            n_only_deletion = "NA"
            n_only_substitution = "NA"
            n_insertion_and_deletion = "NA"
            n_insertion_and_substitution = "NA"
            n_deletion_and_substitution = "NA"
            n_insertion_and_deletion_and_substitution = "NA"

            unmod_pct = "NA"
            mod_pct = "NA"

            if run_status == 'True':
                try:
                    run_data = CRISPRessoShared.load_crispresso_info(crispresso_info['output_folder'])

                except:
                    logger.debug('Could not load CRISPResso run information from ' + crispresso_info['name'] + ' at ' + crispresso_info['output_folder'])
                else:
                    #if running crispresso python 3
                    if 'running_info' in run_data:
                        report_filename = run_data['running_info']['report_filename']
                        report_file_loc = os.path.join(os.path.dirname(crispresso_info['output_folder']),report_filename)
                        if os.path.isfile(report_file_loc):
                            crispresso_results['run_names'].append(name)
                            crispresso_results['run_sub_htmls'][name] = report_file_loc

                        ref_name = run_data['results']['ref_names'][0] #only expect one amplicon sequence
                        n_total = run_data['running_info']['alignment_stats']['N_TOT_READS']
                        n_aligned = run_data['results']['alignment_stats']['counts_total'][ref_name]
                        n_unmod = run_data['results']['alignment_stats']['counts_unmodified'][ref_name]
                        n_mod = run_data['results']['alignment_stats']['counts_modified'][ref_name]
                        n_discarded = run_data['results']['alignment_stats']['counts_discarded'][ref_name]

                        n_insertion = run_data['results']['alignment_stats']['counts_insertion'][ref_name]
                        n_deletion = run_data['results']['alignment_stats']['counts_deletion'][ref_name]
                        n_substitution = run_data['results']['alignment_stats']['counts_substitution'][ref_name]
                        n_only_insertion = run_data['results']['alignment_stats']['counts_only_insertion'][ref_name]
                        n_only_deletion = run_data['results']['alignment_stats']['counts_only_deletion'][ref_name]
                        n_only_substitution = run_data['results']['alignment_stats']['counts_only_substitution'][ref_name]
                        n_insertion_and_deletion = run_data['results']['alignment_stats']['counts_insertion_and_deletion'][ref_name]
                        n_insertion_and_substitution = run_data['results']['alignment_stats']['counts_insertion_and_substitution'][ref_name]
                        n_deletion_and_substitution = run_data['results']['alignment_stats']['counts_deletion_and_substitution'][ref_name]
                        n_insertion_and_deletion_and_substitution = run_data['results']['alignment_stats']['counts_insertion_and_deletion_and_substitution'][ref_name]

                    else: #python 2 version
                        report_filename = run_data['report_filename']
                        report_file_loc = os.path.join(os.path.dirname(crispresso_info['output_folder']),report_filename)
                        if os.path.isfile(report_file_loc):
                            crispresso_results['run_names'].append(name)
                            crispresso_results['run_sub_htmls'][name] = report_file_loc

                        ref_name = run_data['ref_names'][0] #only expect one amplicon sequence
                        n_total = run_data['aln_stats']['N_TOT_READS']
                        n_aligned = run_data['counts_total'][ref_name]
                        n_unmod = run_data['counts_unmodified'][ref_name]
                        n_mod = run_data['counts_modified'][ref_name]
                        n_discarded = run_data['counts_discarded'][ref_name]

                        n_insertion = run_data['counts_insertion'][ref_name]
                        n_deletion = run_data['counts_deletion'][ref_name]
                        n_substitution = run_data['counts_substitution'][ref_name]
                        n_only_insertion = run_data['counts_only_insertion'][ref_name]
                        n_only_deletion = run_data['counts_only_deletion'][ref_name]
                        n_only_substitution = run_data['counts_only_substitution'][ref_name]
                        n_insertion_and_deletion = run_data['counts_insertion_and_deletion'][ref_name]
                        n_insertion_and_substitution = run_data['counts_insertion_and_substitution'][ref_name]
                        n_deletion_and_substitution = run_data['counts_deletion_and_substitution'][ref_name]
                        n_insertion_and_deletion_and_substitution = run_data['counts_insertion_and_deletion_and_substitution'][ref_name]

                    if n_aligned > 0:
                        unmod_pct = 100*n_unmod/float(n_aligned)
                        mod_pct = 100*n_mod/float(n_aligned)

                    if 'fastq_output' not in run_data:
                        raise Exception('Cannot read fastq output for ' + crispresso_info['name'])
                    else:
                        mod_read_count = 0
                        unmod_read_count = 0
                        fastq_output = run_data['fastq_output']
                        if fastq_output.endswith('.gz'):
                            fq_in = gzip.open(fastq_output,'rt')
                        else:
                            fq_in = open(fastq_output,'rt')

                        while (1):
                            fq_id_line   = fq_in.readline().strip()
                            fq_seq_line  = fq_in.readline().strip()
                            fq_plus_line = fq_in.readline()
                            fq_qual_line = fq_in.readline().strip()

                            if not fq_qual_line: break

                            read_id = fq_id_line[1:]
                            if 'CLASS=Reference_MODIFIED' in fq_plus_line:
                                all_crispresso_read_modified_data[read_id] = 'M'
                                mod_read_count += 1
                            else:
                                all_crispresso_read_modified_data[read_id] = 'U'
                                unmod_read_count += 1
                        fq_in.close()

                        if n_mod != mod_read_count:
                            raise Exception("Couldn't get correct number of modified reads from the fastq_output file " + fastq_output)
            new_vals = [name,cut,cut_type,cut_loc,ref_name,run_status,n_printed,n_total,n_aligned,n_unmod,n_mod,n_discarded,n_insertion,n_deletion,n_substitution,n_only_insertion,n_only_deletion,n_only_substitution,n_insertion_and_deletion,n_insertion_and_substitution,n_deletion_and_substitution,n_insertion_and_deletion_and_substitution,cut_amp_seq]
            crispresso_info_fh.write("\t".join([str(x) for x in new_vals])+"\n")

    annotated_final_assignments_file = root + ".annotated_final_assignments"
    annotated_final_class_counts = defaultdict(int)
    with open(final_assignment_file,'r') as f_assignments, open(annotated_final_assignments_file,'w') as fout:
        head = f_assignments.readline().rstrip('\n')
        head_line_els = head.split("\t")
        final_classification_ind = 12
        if head_line_els[final_classification_ind] != "classification":
            raise Exception("Couldn't parse final assignment file " + final_assignment_file)
        fout.write(head+"\tcrispresso_status\n")
        for line in f_assignments:
            line = line.rstrip('\n')
            f_assignments_line_els = line.split("\t")
            crispresso_status = '(not analyzed)'
            this_id = f_assignments_line_els[0]
            if this_id in all_crispresso_read_modified_data:
                if all_crispresso_read_modified_data[this_id] == 'M':
                    crispresso_status = "short indels"
                else:
                    crispresso_status = "no indels"
            if f_assignments_line_els[final_classification_ind] != 'Unidentified':
                annotated_final_class_counts[f_assignments_line_els[final_classification_ind] + ' ' + crispresso_status] += 1
            fout.write("%s\t%s\n"%(line, crispresso_status))

    classification_plot_obj_root = root+".classificationSummary"
    labels = sorted(list(annotated_final_class_counts.keys()))
    values = [annotated_final_class_counts[x] for x in labels]
    sum_values = sum(values)
    with open(classification_plot_obj_root+".txt",'w') as summary:
        summary.write("\t".join(labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
    if sum_values == 0 or suppress_plots:
        classification_plot_obj = None
    else:
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        plot_name = classification_plot_obj_root
        plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')
        plot_count_str = "<br>".join(["%s N=%d"%x for x in zip(labels,values)])

        classification_plot_obj = PlotObject(
                plot_name = classification_plot_obj_root,
                plot_title = 'Read classification including CRISPResso analysis',
                plot_label = 'CRISPResso classification<br>' + plot_count_str + '<br>Reads labeled "(not analyzed)" had too few reads for CRISPResso analysis or indels (based on genome alignment) that were too long. ' + \
                    'You can adjust the minimum number of reads to run CRISPResso with the parameter "crispresso_min_count". You can adjust the maximum indel length for analysis by CRISPResso with the parameter "crispresso_max_indel_size".',
                plot_datas = [
                    ('CRISPResso classification summary',classification_plot_obj_root + ".txt")
                    ]
                )


    with open(crispresso_stats_file,'w') as fout:
        classification_plot_obj_str = "None"
        if classification_plot_obj is not None:
            classification_plot_obj_str = classification_plot_obj.to_json()
        fout.write(classification_plot_obj_str+"\n\n")

        fout.write('Run Name\tRun sub-html\n')
        for run_name in crispresso_results['run_names']:
            fout.write(run_name + "\t" + crispresso_results['run_sub_htmls'][run_name]+"\n")

    if not keep_intermediate:
        for crispresso_info in crispresso_infos:
            reads_file = crispresso_info['reads_file']
            if os.path.exists(reads_file):
                logger.debug('Removing reads file: ' + reads_file)
                os.remove(reads_file)

    return crispresso_results, classification_plot_obj


def getLeftRightMismatchesMDZ(mdz_str):
    """
    Gets the number of mismates on the left and right of an alignment from the MD:Z string

    From the Samtool spec:

    The MD string consists of the following items, concatenated without additional delimiter characters:
      [0-9]+, indicating a run of reference bases that are identical to the corresponding SEQ bases;
      [A-Z], identifying a single reference base that differs from the SEQ base aligned at that position;
      \^[A-Z]+, identifying a run of reference bases that have been deleted in the alignment.

    Args:
        mdz_str: MD:Z string from alignment

    Returns:
        left_matches: number of bp that match on the left side of the alignment
        right_matches
   """
    num_set = set([str(x) for x in range(10)])
    left_matches = 0
    curr_num_str = ""
    for str_index in range(len(mdz_str)):
        if mdz_str[str_index] in num_set:
            curr_num_str += mdz_str[str_index]
        else:
            break

    if curr_num_str == '':
        left_matches = 0
    else:
        left_matches = int(curr_num_str)

    right_matches = 0
    curr_num_str = ""
    for str_index in range(len(mdz_str)-1,-1,-1):
        if mdz_str[str_index] in num_set:
            curr_num_str = mdz_str[str_index] + curr_num_str
        else:
            break

    if curr_num_str == '':
        right_matches = 0
    else:
        right_matches = int(curr_num_str)

    return left_matches,right_matches

def makeTxCountPlot(left_labs = [],
        right_labs = [],
        counts = [],
        fill_cols = None,
        outline_cols=None,
        legend_outline_cols=None,
        legend_outline_labs=None
        ):
    """
    Make a count plot for the most common types of translocations

    Args:
        left_labs: (list) labels for the left side translocations
        right_labs: (list) labels for the right side of translocations
        counts: (list) read counts for each translocation
        fill_cols: (list) colors for blocks
        outline_cols: (list) outline colors for blocks
        legend_outline_cols: (list) colors for outline color legend
        legend_outline_labs: (list) labels for outline color legend

        The length of all params should be the same - one for each translocation to plot
    Returns:
        plt: matplotlib plot object
    """

    num_boxes = len(left_labs)

    if fill_cols is None:
        fill_cols = ['b']*num_boxes

    if outline_cols is None:
        outline_cols = ['None']*num_boxes

    x_width = 2 - 0.1
    y_height = 1 - 0.1
    ys = range(0,num_boxes)[::-1]
    x_start = 2
    count_bar_ydiff = 0.3

    fig, (ax1, ax2) = plt.subplots(1,2,sharey=True, figsize=(12,8))

    boxes = [Rectangle((x_start, y), x_width, y_height)
                      for y in ys]

    ax1.add_collection(PatchCollection(boxes, facecolor=fill_cols, edgecolor=outline_cols, linewidth=2))

    # Add collection to axes
    ax1.set_ylim(0,num_boxes)

    max_right_len = max([len(lab) for lab in right_labs])

    ax1.set_xlim(-2,10)

    for ind in range(num_boxes):
        ax1.text(x_start-0.1,ys[ind]+y_height/2,left_labs[ind],ha='right',va='center')
        ax1.text(4,ys[ind]+y_height/2,right_labs[ind],ha='left',va='center')

    ax1.axis('off')

    if legend_outline_labs is not None:
        legend_patches = []
        for col,lab in zip(legend_outline_cols,legend_outline_labs):
            legend_patches.append(Patch(facecolor='None',edgecolor=col,label=lab))
        ax1.legend(handles=legend_patches,loc="lower center", bbox_to_anchor=(0.5, -0.2))

    rects = []
    for ind in range(num_boxes):
        val = max(1,counts[ind]) #min value is 1 or matplotlib flips out
        rects.append(Rectangle((1,ys[ind]+count_bar_ydiff),val,y_height-(count_bar_ydiff*2)))
        ax2.text(x=val,y=ys[ind]+y_height/2,s=" " + str(counts[ind]),ha='left',va='center')

    pc = PatchCollection(rects)
    ax2.add_collection(pc)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_frame_on(False)

    ax2.set_xscale("log")
    ax2.set_xlim(1,max(5,max(counts)))
    ax2.set_xlabel('Number of reads')

    return plt

def make_final_summary(root, num_reads_input, post_dedup_count, post_filter_on_primer_read_count, final_read_count, discarded_read_counts, classification_read_counts, classification_indel_read_counts,suppress_plots=False):
    """
    Make final summary plot object showing number of reads deduplicated, filtered, and classified, etc.

    Args:
        root: where file should be written to
        num_reads_input (int): Number of reads in input
        post_dedup_count (int): Number of reads post-deduplication (of initial UMIs)
        post_filter_on_primer_read_count (int): Number of reads post filtering on primers
        final_read_count (int): Number of reads in final analysis
        discarded_read_counts (list): List of tuples (why read was discarded, number of reads)
        classification_read_counts (list): List of tuples (read classification, number of reads)
        classification_indel_read_counts (list): List of tuples (read classification, number of reads) including 'short indels' and 'long indels' categories
        suppress_plots: if true, plotting will be suppressed

    Returns:
        final_summary_file: file containing summary of results
        final_summary_plot_obj: plot shwoing final summary
    """
    final_summary_str = 'Total input reads: ' + str(num_reads_input) + '\n'
    final_summary_head = ['total_input_reads']
    final_summary_vals = [num_reads_input]

    post_dedup_count_pct = None
    if num_reads_input > 0:
        post_dedup_count_pct=round(100*post_dedup_count/num_reads_input,2)
    final_summary_str += '\tSurvived initial UMI deduplication: %d/%d (%s%%)\n'%(post_dedup_count,num_reads_input,post_dedup_count_pct)
    final_summary_head.append('post_dedup_reads')
    final_summary_vals.append(post_dedup_count)

    post_filter_on_primer_read_count_pct = None
    if post_dedup_count > 0:
        post_filter_on_primer_read_count_pct =round(100* post_filter_on_primer_read_count/post_dedup_count,2)
    final_summary_str += '\tSurvived filtering for primer/origin presence: %d/%d (%s%%)\n'%(post_filter_on_primer_read_count,post_dedup_count,post_filter_on_primer_read_count_pct)
    final_summary_head.append('post_primer_filter_reads')
    final_summary_vals.append(post_filter_on_primer_read_count)

    discarded_total = sum([x[1] for x in discarded_read_counts])
    survived_filtering_count = post_filter_on_primer_read_count - discarded_total

    survived_filtering_count_pct = None
    if post_filter_on_primer_read_count > 0:
        survived_filtering_count_pct =round(100* survived_filtering_count/post_filter_on_primer_read_count,2)
    final_summary_str += '\tSurvived quality filtering: %d/%d (%s%%)\n'%(survived_filtering_count,post_filter_on_primer_read_count,survived_filtering_count_pct)
    final_summary_head.append('post_quality_filter_reads')
    final_summary_vals.append(survived_filtering_count)

    discarded_total_pct = None
    if post_filter_on_primer_read_count > 0:
        discarded_total_pct =round(100* discarded_total/post_filter_on_primer_read_count,2)
    final_summary_str += '\tFailed quality filtering: %d/%d (%s%%)\n'%(discarded_total,post_filter_on_primer_read_count,discarded_total_pct)
    final_summary_head.append('discarded_quality_filter_reads')
    final_summary_vals.append(discarded_total)

    for (category,category_count) in discarded_read_counts:
        category_pct = None
        if discarded_total > 0:
            category_pct =round(100* category_count/discarded_total,2)
        final_summary_str += '\t\t%s: %d/%d (%s%%)\n'%(category,category_count,discarded_total,category_pct)
        final_summary_head.append(category)
        final_summary_vals.append(category_count)

    final_read_count_pct = None
    if post_filter_on_primer_read_count > 0:
        final_read_count_pct =round(100* final_read_count/post_filter_on_primer_read_count,2)
    final_summary_str += '\tReads used for final analysis: %d/%d (%s%%)\n'%(final_read_count,post_filter_on_primer_read_count,final_read_count_pct)
    final_summary_head.append('analyzed_read_count')
    final_summary_vals.append(final_read_count)

    for (category,category_count) in classification_read_counts:
        category_pct = None
        if final_read_count > 0:
            category_pct =round(100* category_count/final_read_count,2)
        final_summary_str += '\t\t%s: %d/%d (%s%%)\n'%(category,category_count,final_read_count,category_pct)
        final_summary_head.append(category)
        final_summary_vals.append(category_count)

    final_summary_str += '\tReads used for final analysis (including breakdown of indels): %d/%d (%s%%)\n'%(final_read_count,post_filter_on_primer_read_count,final_read_count_pct)

    for (category,category_count) in classification_indel_read_counts:
        category_pct = None
        if final_read_count > 0:
            category_pct =round(100* category_count/final_read_count,2)
        final_summary_str += '\t\t%s: %d/%d (%s%%)\n'%(category,category_count,final_read_count,category_pct)
        final_summary_head.append(category)
        final_summary_vals.append(category_count)

    #plot summary
    filter_labels = ['Deduplicated by initial UMI deduplication','Did not contain origin/primer seqeunce','Failed quality filtering','Used for final analyses']
    values = [num_reads_input - post_dedup_count,post_dedup_count - post_filter_on_primer_read_count,post_filter_on_primer_read_count - final_read_count, final_read_count]
    summary_plot_obj_root = root
    with open(summary_plot_obj_root+".txt","w") as fout:
        fout.write("\t".join(final_summary_head)+"\n")
        fout.write("\t".join([str(x) for x in final_summary_vals])+"\n")
        fout.write(final_summary_str + "\n")
    if suppress_plots:
        summary_plot_obj = None
    else:
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(filter_labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(filter_labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Read Assignment Summary')
        plt.savefig(summary_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(summary_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(filter_labels,values)])
        plot_label = 'Assignment summary for N=' + str(num_reads_input) + ' input reads:'
        summary_plot_obj = PlotObject(
                plot_name = summary_plot_obj_root,
                plot_title = 'Read summary',
                plot_label = plot_label + '<br>'+plot_count_str,
                plot_datas = [
                    ('Read assignment summary',summary_plot_obj_root + ".txt")
                    ]
                )
    return summary_plot_obj_root + ".txt", summary_plot_obj

class PlotObject:
    """
    Holds information for plots for future output, namely:
        the plot name: root of plot (name.pdf and name.png should exist)
        the plot title: title to be shown to user
        the plot label: label to be shown under the plot
        the plot data: array of (tuple of display name and file name)
        the plot order: int specifying the order to display on the report (lower numbers are plotted first, followed by higher numbers)
    """
    def __init__(self,plot_name,plot_title,plot_label,plot_datas,plot_order=50):
        self.name = plot_name
        self.title = plot_title
        self.label = plot_label
        self.datas = plot_datas
        self.order = plot_order

    def to_json(self):
        obj = {
                'plot_name':self.name,
                'plot_title':self.title,
                'plot_label':self.label,
                'plot_datas':self.datas,
                'plot_order':self.order
                }
        obj_str = json.dumps(obj,separators=(',',':'))
        return obj_str

    #construct from a json string
    @classmethod
    def from_json(cls, json_str):
        obj = json.loads(json_str)
        return cls(plot_name=obj['plot_name'],
                plot_title=obj['plot_title'],
                plot_label=obj['plot_label'],
                plot_datas=obj['plot_datas'],
                plot_order=obj['plot_order'])

    def __str__(self):
        return 'Plot object with name ' + self.name

    def __repr__(self):
        return f'PlotObject(name={self.name}, title={self.title}, label={self.label}, datas={self.datas} order={self.order})'



def make_report(report_file,report_name,crisprlungo_folder,
            crispresso_run_names,crispresso_sub_html_files,
            summary_plot_objects=[]
        ):
    """
    Makes an HTML report for a CRISPRlungo run

    Parameters:
    report_file: path to the output report
    report_name: description of report type to be shown at top of report
    crisprlungo_folder (string): absolute path to the crisprlungo output

    crispresso_run_names (arr of strings): names of crispresso runs
    crispresso_sub_html_files (dict): dict of run_name->file_loc

    summary_plot_objects (list): list of PlotObjects to plot
    """

    logger = logging.getLogger('CRISPRlungo')
    ordered_plot_objects = sorted(summary_plot_objects,key=lambda x: x.order)

    html_str = """
<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <title>"""+report_name+"""</title>

    <!-- Bootstrap core CSS -->
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootswatch/4.5.2/flatly/bootstrap.min.css" integrity="sha384-qF/QmIAj5ZaYFAeQcrQ6bfVMAh4zZlrGwTPY7T/M+iTTLJqJBJjwwnsE5Y0mV7QK" crossorigin="anonymous">
  </head>

  <body>
<style>
html,
body {
  height: 100%;
}

body {
  padding-bottom: 40px;
  background-color: #f5f5f5;
}

.navbar-fixed-left {
  width: 200px;
  position: fixed;
  border-radius: 0;
  height: 100%;
  padding: 10px;
}

.navbar-fixed-left .navbar-nav > li {
  /*float: none;   Cancel default li float: left */
  width: 160px;
}

</style>

<nav class="navbar navbar-fixed-left navbar-dark bg-dark" style="overflow-y:auto">
     <a class="navbar-brand" href="#">CRISPRlungo</a>
      <ul class="nav navbar-nav me-auto">
"""
    for idx,plot_obj in enumerate(ordered_plot_objects):
        html_str += """        <li class="nav-item">
          <a class="nav-link active" href="#plot"""+str(idx)+"""">"""+plot_obj.title+"""
          </a>
        </li>"""
    if len(crispresso_run_names) > 0:
        html_str += """        <li class="nav-item">
          <a class="nav-link active" href="#crispresso_output">CRISPResso Output
          </a>
        </li>"""
    html_str += """      </ul>
</nav>
<div class='container'>
<div class='row justify-content-md-center'>
<div class='col-8'>
    <div class='text-center pb-4'>
    <h1 class='display-3 pt-5'>CRISPRlungo</h1><hr><h2>"""+report_name+"""</h2>
    </div>
"""

    data_path = ""
    for idx,plot_obj in enumerate(ordered_plot_objects):
        plot_path = plot_obj.name
        plot_path = os.path.basename(plot_path)
        plot_str = "<div class='card text-center mb-2' id='plot"+str(idx)+"'>\n\t<div class='card-header'>\n"
        plot_str += "<h5>"+plot_obj.title+"</h5>\n"
        plot_str += "</div>\n"
        plot_str += "<div class='card-body'>\n"
        plot_str += "<a href='"+data_path + plot_path+".pdf'><img src='"+data_path + plot_path + ".png' width='80%' ></a>\n"
        plot_str += "<label>"+plot_obj.label+"</label>\n"
        for (plot_data_label,plot_data_path) in plot_obj.datas:
            plot_data_path = os.path.basename(plot_data_path)
            plot_str += "<p class='m-0'><small>Data: <a href='"+data_path+plot_data_path+"'>" + plot_data_label + "</a></small></p>\n";
        plot_str += "</div></div>\n";
        html_str += plot_str

    if len(crispresso_run_names) > 0:
        run_string = """<div class='card text-center mb-2' id='crispresso_output'>
          <div class='card-header'>
            <h5>CRISPResso Output</h5>
          </div>
          <div class='card-body p-0'>
            <div class="list-group list-group-flush">
            """
        for crispresso_run_name in crispresso_run_names:
            crispresso_run_names,crispresso_sub_html_files,
            run_string += "<a href='"+data_path+crispresso_sub_html_files[crispresso_run_name]+"' class='list-group-item list-group-item-action'>"+crispresso_run_name+"</a>\n"
        run_string += "</div></div></div>"
        html_str += run_string

    html_str += """
                </div>
            </div>
        </div>
    </body>
</html>
"""
    with open(report_file,'w') as fo:
        fo.write(html_str)
    logger.info('Wrote ' + report_file)

class SSWPrimerAlign:
    """
    Class for storing objects for SSW alignment

    Attributes:
        ssw_obj: ssw object for alignment
        list_letters: list of possible letters for alignment
        dict_letter_lookup: lookup of letter>number
        matrix: alignment matrix scores for matches/mismatches
        qProfile: profile of primer sequence for ssw alignment
        nFlag: flag for ssw alignment
        match_score: match score for ssw alignment (should be positive int)
        mismatch_penalty: mismatch score for ssw alignment (should be positive int)
        gap_open_penalty: gap open score for ssw alignment (should be positive int)
        gap_extend_penalty: gap extend score for ssw alignment (should be positive int)
        primer_seq: primer seq to look for
    
    """
    def __init__(self,ssw_obj,primer_seq,list_letters=['A','C','G','T','N'],nFlag=1,match_score=2,mismatch_penalty=5,gap_open_penalty=5,gap_extend_penalty=5):
        """
        Args:
            ssw_obj: ssw object for alignment
            primer_seq: primer sequence which will be searched for in reads
            list_letters: list of possible letters for alignment
            nFlag: flag for ssw alignment - 0 returns starts only, 1 returns both starts and ends
            match_score: match score for ssw alignment (should be positive int)
            mismatch_penalty: mismatch score for ssw alignment (should be positive int)
            gap_open_penalty: gap open score for ssw alignment (should be positive int)
            gap_extend_penalty: gap extend score for ssw alignment (should be positive int)
        """        
        list_letters = ['A', 'C', 'G', 'T', 'N']
        dict_letter_lookup = {}
        for i,letter in enumerate(list_letters):
            dict_letter_lookup[letter] = i
            dict_letter_lookup[letter.lower()] = i
        num_letters = len(list_letters)
        lScore = [0 for i in range(num_letters**2)]
        for i in range(num_letters-1):
            for j in range(num_letters-1):
                if list_letters[i] == list_letters[j]:
                    lScore[i*num_letters+j] = match_score
                else:
                    lScore[i*num_letters+j] = -mismatch_penalty

        # translate score matrix to ctypes
        matrix = (len(lScore) * ctypes.c_int8) ()
        matrix[:] = lScore

        qNum = to_int(primer_seq, list_letters, dict_letter_lookup)
        score_size = 2 #estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
        qProfile = ssw_obj.ssw_init(qNum, ctypes.c_int32(len(primer_seq)), matrix, len(list_letters), score_size)

        self.ssw_obj = ssw_obj
        self.list_letters = list_letters
        self.dict_letter_lookup = dict_letter_lookup
        self.nFlag = nFlag
        self.matrix = matrix
        self.qProfile = qProfile
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open_penalty = gap_open_penalty
        self.gap_extend_penalty = gap_extend_penalty
        self.primer_seq = primer_seq
        self.qNum = qNum
        self.c_primer_seq_len = ctypes.c_int32(len(primer_seq))
        self.len_list_letters = len(list_letters)
        self.score_size = score_size


    def align(self,read):
        rNum = to_int(read, self.list_letters, self.dict_letter_lookup)
        qProfile = self.ssw_obj.ssw_init(self.qNum, self.c_primer_seq_len, self.matrix, self.len_list_letters, self.score_size)

        # the last three parameters to this function are not used (limits on scores returned, limits on location, mask length for second-best alignments)
        ssw_res = self.ssw_obj.ssw_align(self.qProfile, rNum, ctypes.c_int32(len(read)), self.gap_open_penalty, self.gap_extend_penalty, self.nFlag, 0, 0, 0)
        return ssw_res

    def destroy(self,ssw_res):
        self.ssw_obj.align_destroy(ssw_res)

def trimPrimersSingle(fastq_r1,fastq_r1_trimmed,min_primer_aln_score,min_primer_length,ssw_align_primer):
    """
    Trims the primer from single-end input reads, only keeps reads with the primer present in R1

    Args:
        root: root for written files
        fastq_r1: R1 reads to trim
        fastq_r1_trimmed: output file to write R1 to
        min_primer_aln_score: minimum score for alignment between primer/origin sequence and read sequence
        min_primer_length: minimum length of sequence that matches between the primer/origin sequence and the read sequence
        ssw_align_primer: SSWPrimerAlignment object for ssw alignment of primer

    Returns:
        post_trim_read_count: number of reads after trimming (primer found)
        too_short_read_count: number of reads where primer seq found was too short
        untrimmed_read_count: number of reads untrimmed (no primer found)
    """
    if fastq_r1.endswith('.gz'):
        f1_in = gzip.open(fastq_r1,'rt')
    else:
        f1_in = open(fastq_r1,'rt')

    f1_out = gzip.open(fastq_r1_trimmed, 'wt')

    tot_read_count = 0
    post_trim_read_count = 0
    too_short_read_count = 0
    untrimmed_read_count = 0

    while (1):
        f1_id_line   = f1_in.readline().strip()
        f1_seq_line  = f1_in.readline().strip()
        f1_plus_line = f1_in.readline()
        f1_qual_line = f1_in.readline().strip()

        if not f1_qual_line : break
        if not f1_plus_line.startswith("+"):
            raise Exception("Fastq %s cannot be parsed (%s%s%s%s) "%(fastq_r1,f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))
        tot_read_count += 1

        primer_seq, trimmed_seq, trimmed_primer_pos = trimLeftPrimerFromRead(f1_seq_line,ssw_align_primer,min_primer_aln_score)
        if len(primer_seq) > min_primer_length:
            post_trim_read_count += 1
            new_f1_qual_line = f1_qual_line[len(primer_seq):]
            new_f1_id_line = f1_id_line + ' primer_len=' + str(trimmed_primer_pos + 1) # trimmed primer len is (trimmed_primer_pos + 1)

            f1_out.write(new_f1_id_line + "\n" + trimmed_seq + "\n" + f1_plus_line + new_f1_qual_line + "\n")
        elif len(primer_seq) > 0:
            too_short_read_count += 1
        else:
            untrimmed_read_count += 1

    #finished iterating through fastq file
    f1_in.close()
    f1_out.close()
    return(post_trim_read_count, too_short_read_count, untrimmed_read_count)

def trimPrimersPair(fastq_r1,fastq_r2,fastq_r1_trimmed,fastq_r2_trimmed,min_primer_aln_score,min_primer_length,ssw_align_primer,ssw_align_primer_rc):
    """
    Trims the primer from paired input reads, only keeps reads with the primer present in R1

    Args:
        root: root for written files
        fastq_r1: R1 reads to trim
        fastq_r2: R2 reads to trim
        fastq_r1_trimmed: output file to write R1 to
        fastq_r2_trimmed: output file to write R2 to
        min_primer_aln_score: minimum score for alignment between primer/origin sequence and read sequence
        min_primer_length: minimum length of sequence that matches between the primer/origin sequence and the read sequence
        ssw_align_primer: SSWPrimerAlignment object for ssw alignment of primer
        ssw_align_primer_rc: SSWPrimerAlignment object for ssw alignment of reverse complement primer (to r2)

    Returns:
        post_trim_read_count: number of reads after trimming (primer found)
        too_short_read_count: number of reads where primer seq found was too short
        untrimmed_read_count: number of reads untrimmed (no primer found)
    """
    if fastq_r1.endswith('.gz'):
        f1_in = gzip.open(fastq_r1,'rt')
    else:
        f1_in = open(fastq_r1,'rt')

    if fastq_r2.endswith('.gz'):
        f2_in = gzip.open(fastq_r2,'rt')
    else:
        f2_in = open(fastq_r2,'rt')

    f1_out = gzip.open(fastq_r1_trimmed, 'wt')
    f2_out = gzip.open(fastq_r2_trimmed, 'wt')

    tot_read_count = 0
    post_trim_read_count = 0
    too_short_read_count = 0
    untrimmed_read_count = 0

    while (1):
        f1_id_line   = f1_in.readline().strip()
        f1_seq_line  = f1_in.readline().strip()
        f1_plus_line = f1_in.readline()
        f1_qual_line = f1_in.readline().strip()

        if not f1_qual_line : break
        if not f1_plus_line.startswith("+"):
            raise Exception("Fastq %s cannot be parsed (%s%s%s%s) "%(fastq_r1,f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))
        tot_read_count += 1

        f2_id_line   = f2_in.readline().strip()
        f2_seq_line  = f2_in.readline().strip()
        f2_plus_line = f2_in.readline()
        f2_qual_line = f2_in.readline().strip()

        primer_seq, trimmed_seq, trimmed_primer_pos = trimLeftPrimerFromRead(f1_seq_line,ssw_align_primer,min_primer_aln_score)
        if len(primer_seq) > min_primer_length:
            post_trim_read_count += 1
            new_f1_qual_line = f1_qual_line[len(primer_seq):]
            new_f1_id_line = f1_id_line + ' primer_len=' + str(trimmed_primer_pos + 1) # trimmed primer len is (trimmed_primer_pos + 1)

            f2_trimmed_seq, f2_primer_seq = trimRightPrimerFromRead(f2_seq_line,ssw_align_primer_rc,min_primer_aln_score)
            new_f2_qual_line = f2_qual_line[len(f2_primer_seq):]
            f1_out.write(new_f1_id_line + "\n" + trimmed_seq + "\n" + f1_plus_line + new_f1_qual_line + "\n")
            f2_out.write(f2_id_line + "\n" + f2_trimmed_seq + "\n" + f2_plus_line + new_f2_qual_line + "\n")
        elif len(primer_seq) > 0:
            too_short_read_count += 1
        else:
            untrimmed_read_count += 1

    #finished iterating through fastq file
    f1_in.close()
    f2_in.close()
    f1_out.close()
    f2_out.close()
    return(post_trim_read_count, too_short_read_count, untrimmed_read_count)

def trimLeftPrimerFromRead(read, ssw_align, min_score=40, debug=False):
    """
    Trims a primer from a given read
    E.g. for --primer--read-- returns the portion that aligns to the primer (and before) as the trimmed_primer_seq and to the right as the trimmed_read_seq

    Args:
        primer: string
        read: string
        ssw_align: ssw alignment object
        min_score: minimum alignment score between primer and sequence
        debug: print some debug statements

    Returns:
        trimmed_primer_seq: portion of read that was aligned to the primer
        trimmed_read_seq: portion of the read after the primer
        trimmed_primer_pos: end index of primer that was trimmed. Note the length of the primer is this value + 1
    """

    ssw_res = ssw_align.align(read)
    if debug:
        print('================')
        print('TRIM LEFT PRIMER')
        print(f'{min_score=}')
        print(f'{read=}')
        print(f'{ssw_align.primer_seq=}')
        print(f'{ssw_res.contents.nScore=}')
        print(f'{ssw_res.contents.nRefBeg=}')
        print(f'{ssw_res.contents.nRefEnd=}')
        print(f'{ssw_res.contents.nQryBeg=}')
        print(f'{ssw_res.contents.nQryEnd=}')
    if ssw_res.contents.nScore < min_score:
        ssw_align.destroy(ssw_res)
        return '',read,0
    else:
        trimmed_primer_seq = read[:ssw_res.contents.nRefEnd+1]
        trimmed_read_seq = read[ssw_res.contents.nRefEnd+1:]
        ssw_align.destroy(ssw_res)
        return trimmed_primer_seq,trimmed_read_seq,ssw_res.contents.nQryEnd

def trimRightPrimerFromRead(read, ssw_align, min_score=40, debug=False):
    """
    Trims a primer from a given read
    E.g. for --read--primer-- returns the portion that aligns to the primer (and after) as the trimmed_primer_seq and to the left as the trimmed_read_seq

    Args:
        read: string
        ssw_align: ssw alignment object
        mismatch_score: mismatch/gap score for alignment

    Returns:
        trimmed_read_seq: portion of the read after the primer
        trimmed_primer_seq: portion of read that was aligned to the primer

    """

    ssw_res = ssw_align.align(read)
    if debug:
        print('=================')
        print('TRIM RIGHT PRIMER')
        print(f'{min_score=}')
        print(f'{read=}')
        print(f'{ssw_align.primer_seq=}')
        print(f'{ssw_res.contents.nScore=}')
        print(f'{ssw_res.contents.nRefBeg=}')
        print(f'{ssw_res.contents.nRefEnd=}')
        print(f'{ssw_res.contents.nQryBeg=}')
        print(f'{ssw_res.contents.nQryEnd=}')

    if ssw_res.contents.nScore < min_score:
        ssw_align.destroy(ssw_res)
        return read,''
    else:
        trimmed_read_seq = read[:ssw_res.contents.nRefBeg]
        trimmed_primer_seq = read[ssw_res.contents.nRefBeg:]
        ssw_align.destroy(ssw_res)
        return trimmed_read_seq,trimmed_primer_seq

def to_int(seq, list_letters, dict_letter_lookup):
    """
    Translate a letter sequence into numbers for ssw alignment
    Args:
        seq: string, letters to translate
        list_letters: list of all letters
        dict_letter_lookup: dict of letter > number for lookup

    Returns:
        array of numbers representing sequence
    """
    num_decl = len(seq) * ctypes.c_int8
    num = num_decl()
    for i,letter in enumerate(seq):
        try:
            n = dict_letter_lookup[letter]
        except KeyError:
            n = dict_letter_lookup[list_letters[-1]]
        finally:
            num[i] = n

    return num

def get_guide_match_from_aln(guide_seq_aln, ref_seq_aln, cut_ind):
    """
    Count the number of matches and mismatches from two alignment strings, as well as the index of the cut site

    Args:
        guide_seq_aln (str): The guide aligned sequence e.g.     --ATTA-
        ref_seq_aln (str): The reference aligned sequence e.g.   AAATTGGGG
        cut_ind (int): The zero-based index in the the guide seq string of the base to the left where a cut would occur (e.g. for Cas9 and reverse-complement guide PAMGGGGGGGGGGG the cut_ind would be given as 5)

    Returns:
        match_count: number of bp in longest match
        mismatch_count: number of bp that mismatch
        gap_count: number of gaps in longest match
        aln_cut_ind: the base in the reference corresponding to the cut_ind in the guide
    """

    if len(guide_seq_aln) != len(ref_seq_aln):
        raise Exception('Sequence lengths must be the same')
    first_nongap_ind = -1
    last_nongap_ind = -1
    guide_ind = 0
    ref_ind = 0
    ref_ind_cut_site = -1
    for i in range(len(guide_seq_aln)):
        if guide_seq_aln[i] == '-':
            pass
        else:
            guide_ind += 1
            last_nongap_ind = i
            if first_nongap_ind == -1:
                first_nongap_ind = i
        if ref_seq_aln[i] != '-':
            ref_ind += 1
        if guide_ind == cut_ind:
            ref_ind_cut_site = ref_ind

    match_count = 0
    mismatch_count = 0
    gap_count = 0
    guide_ind = 0
    for i in range(first_nongap_ind,last_nongap_ind+1):
        if guide_seq_aln[i] == '-' or ref_seq_aln[i] == '-':
            gap_count += 1
        elif guide_seq_aln[i] == ref_seq_aln[i]:
            match_count += 1
        else:
            mismatch_count += 1

    potential_guide_in_ref = ref_seq_aln[first_nongap_ind:last_nongap_ind+1].replace("-","")

    return match_count, mismatch_count, gap_count, ref_ind_cut_site,potential_guide_in_ref



if __name__ == "__main__":
    settings = parse_settings(sys.argv)
    try:
        processCRISPRlungo(settings)
    except Exception as e:
        logger = logging.getLogger('CRISPRlungo')
        if logger.isEnabledFor(logging.DEBUG):
            logger.error(e, exc_info=True)
        else:
            logger.error(e)
        if '--debug' in sys.argv:
            raise e
        else:
            print(str(e))
            sys.exit(1)
