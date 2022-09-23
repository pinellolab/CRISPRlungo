import argparse
from collections import defaultdict
import ctypes
import gzip
import io
import json
import logging
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle,Patch
import multiprocessing as mp
import numpy as np
import os
import re
import subprocess
import sys
from lib import ssw_lib
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPResso2Align

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

__version__ = "v0.1.2"

def main(settings, logger):

    #data structures for plots for report
    summary_plot_objects=[]  # list of PlotObjects for plotting

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


    curr_r1_file = settings['fastq_r1'] #keep track of current input files (through trimming, etc.)
    curr_r2_file = settings['fastq_r2']

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

    dedup_input_UMI_count = 0
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
            use_cutadapt = settings['primer_filter_use_cutadapt'],
            keep_intermediate = settings['keep_intermediate'],
            can_use_previous_analysis = settings['can_use_previous_analysis'],
            )
    curr_r1_file = filtered_on_primer_fastq_r1
    curr_r2_file = filtered_on_primer_fastq_r2

    if filter_on_primer_plot_obj is not None:
        filter_on_primer_plot_obj.order = 1
        summary_plot_objects.append(filter_on_primer_plot_obj)

    #perform alignment
    (genome_mapped_bam
            ) = align_reads(
                root = settings['root']+'.genomeAlignment',
                fastq_r1 = curr_r1_file,
                fastq_r2 = curr_r2_file,
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

    (final_assignment_file,final_read_ids_for_cuts,final_cut_counts,cut_classification_lookup,
        chr_aln_plot_obj,tlen_plot_obj,deduplication_plot_obj,tx_order_plot_obj,tx_count_plot_obj,classification_plot_obj,classification_indel_plot_obj,
        origin_indel_hist_plot_obj,origin_inversion_hist_plot_obj,origin_deletion_hist_plot_obj,origin_indel_depth_plot_obj,
        r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,r1_r2_no_support_dist_plot_obj) = make_final_read_assignments(
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
                cut_merge_dist = settings['novel_cut_merge_distance'],
                collapse_to_homology_dist = settings['homology_cut_merge_distance'],
                arm_min_matched_start_bases = settings['arm_min_matched_start_bases'],
                genome_map_resolution = 1000000,
                suppress_dedup_input_based_on_aln_pos_and_UMI = settings['suppress_dedup_input_based_on_aln_pos_and_UMI'],
                discard_reads_without_r2_support = settings['discard_reads_without_r2_support'],
                samtools_command = settings['samtools_command'],
                keep_intermediate = settings['keep_intermediate'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )

    chr_aln_plot_obj.order=20
    summary_plot_objects.append(chr_aln_plot_obj)

    if tlen_plot_obj is not None:
        tlen_plot_obj.order=15
        summary_plot_objects.append(tlen_plot_obj)

    if tlen_plot_obj is not None:
        deduplication_plot_obj.order=16
        summary_plot_objects.append(deduplication_plot_obj)

    classification_plot_obj.order=36
    summary_plot_objects.append(classification_plot_obj)

    classification_indel_plot_obj.order=37
    summary_plot_objects.append(classification_indel_plot_obj)

    if tx_order_plot_obj is not None:
        tx_order_plot_obj.order= 38
        summary_plot_objects.append(tx_order_plot_obj)

    if tx_count_plot_obj is not None:
        tx_count_plot_obj.order= 40
        summary_plot_objects.append(tx_count_plot_obj)

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
    if r1_r2_no_support_dist_plot_obj is not None:
        r1_r2_no_support_dist_plot_obj.order=31
        summary_plot_objects.append(r1_r2_no_support_dist_plot_obj)

    #crispresso_infos: dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command
    crispresso_infos = prep_crispresso2(
                root = settings['root']+'.CRISPResso_r1',
                input_fastq_file = settings['fastq_r1'],
                read_ids_for_crispresso = final_read_ids_for_cuts,
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
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )
    if crispresso_classification_plot_obj is not None:
        crispresso_classification_plot_obj.order=40
        summary_plot_objects.append(crispresso_classification_plot_obj)

    make_report(report_file=settings['root']+".html",
            report_name = 'Report',
            crisprlungo_folder = '',
            crispresso_run_names = crispresso_results['run_names'],
            crispresso_sub_html_files = crispresso_results['run_sub_htmls'],
            summary_plot_objects = summary_plot_objects,
            )


    logger.info('Successfully completed!')

    # FINISHED

def parse_settings(args):
    """
    Parses settings from the command line
        First parses from the settings file, then parses from command line

    param:
        args: command line arguments

    returns:
        settings: dict of parsed settings
    """
    parser = argparse.ArgumentParser(description='CRISPRlungo: Analyzing unidirectional sequencing of genome editing', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version="%(prog)s "+__version__)
    parser.add_argument('settings_file', nargs='*', help='Tab-separated settings file')

    parser.add_argument('--debug', action='store_true', help='Tab-separated settings file')
    parser.add_argument('--root', type=str, default=None, help='Output directory file root')
    parser.add_argument('--keep_intermediate',action='store_true',help='If true, intermediate files are not deleted')


    parser.add_argument('--guide_sequences', nargs='*', help='Spacer sequences of guides (multiple guide sequences are separated by spaces). Spacer sequences must be provided without the PAM sequence, but oriented so the PAM would immediately follow the provided spacer sequence', default=[])
    parser.add_argument('--cuts','--cut_sites', nargs='*', help='Cut sites in the form chr1:234 (multiple cuts are separated by spaces)', default=[])
    parser.add_argument('--on_target_cut_sites', nargs='*', help='On-target cut sites in the form chr1:234 (multiple cuts are separated by spaces)', default=[])
    parser.add_argument('--cut_classification_annotations', nargs='*', help='User-customizable annotations for cut products in the form: chr1:234:left:Custom_label (multiple annotations are separated by spaces)', default=[])
    parser.add_argument('--cleavage_offset', type=int, help='Position where cleavage occurs, for in-silico off-target search (relative to end of spacer seq -- for Cas9 this is -3)', default=-3)

    parser.add_argument('--genome', help='Genome sequence file for alignment. This should point to a file ending in ".fa", and the accompanying index file (".fai") should exist.', default=None)
    parser.add_argument('--bowtie2_genome', help='Bowtie2-indexed genome file.',default=None)

    parser.add_argument('--fastq_r1', help='Input fastq r1 file. Reads in this file are primed from the provided primer sequence', default=None)
    parser.add_argument('--fastq_r2', help='Input fastq r2 file', default=None)
    parser.add_argument('--fastq_umi', help='Input fastq umi file', default=None)

    parser.add_argument('--novel_cut_merge_distance', type=int, help='Novel cut sites discovered within this distance (bp) from each other (and not within homology_cut_merge_distance to a known/provided cut site or a site with homology to guide_sequences) will be merged into a single cut site. Variation in the cut sites or in the fragments produced may produce clusters of cut sites in a certain region. This parameter will merge novel cut sites within this distance into a single cut site.', default=100)
    parser.add_argument('--homology_cut_merge_distance', type=int, help='Novel cut sites discovered within this distance (bp) with a known/provided/homologous site will be merged to that site. Homologous sites are defined as those that have homology to guide_sequences. Novel cut sites farther than homology_cut_merge_distance will be merged into novel cut sites based on the parameter novel_cut_merge_distance.', default=10000)

    #for finding offtargets with casoffinder
    ot_group = parser.add_argument_group('In silico off-target search parameters')
    ot_group.add_argument('--PAM', type=str, help='PAM for in-silico off-target search', default=None)
    ot_group.add_argument('--casoffinder_num_mismatches', type=int, help='If greater than zero, the number of Cas-OFFinder mismatches for in-silico off-target search. If this value is zero, Cas-OFFinder is not run', default=0)

    #specify primer filtering information
    p_group = parser.add_argument_group('Primer and filtering parameters and settings')
    p_group.add_argument('--suppress_primer_filtering', help='Whether to filter reads for the presence of the primer/origin sequence. If set, all reads will be considered (but alignments may not ', action='store_true')
    p_group.add_argument('--primer_seq', type=str, help='Sequence of primer',default=None)
    p_group.add_argument('--min_primer_aln_score', type=int, help='Minimum primer/origin alignment score for trimming.',default=40)
    p_group.add_argument('--min_primer_length', type=int, help='Minimum length of sequence required to match between the primer/origin and read sequence',default=30)
    p_group.add_argument('--primer_filter_use_cutadapt', help='Whether to use cutadapt to identify primer/origin sequence in reads. If true, cutadapt is used. If false, custom lookhead algorithm is used.', action='store_true')
    p_group.add_argument('--min_read_length', type=int, help='Minimum length of read after all filtering',default=30)
    p_group.add_argument('--transposase_adapter_seq', type=str, help='Transposase adapter sequence to be trimmed from reads',default='CTGTCTCTTATACACATCTGACGCTGCCGACGA')


    #min alignment cutoffs for alignment to each arm/side of read
    a_group = parser.add_argument_group('Alignment cutoff parameters')
    a_group.add_argument('--arm_min_matched_start_bases', type=int, help='Number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each "side" of the alignment. E.g. if a read aligns to a genomic location, the first and last arm_min_matched_start_bases of the read would have to match exactly to the aligned location.', default=10)
    a_group.add_argument('--ignore_n', help='If set, "N" bases will be ignored. By default (False) N bases will count as mismatches in the number of bases required to match at each arm/side of the read', action='store_true')

    #CRISPResso settings
    c_group = parser.add_argument_group('CRISPResso settings')
    c_group.add_argument('--crispresso_min_count', type=int, help='Min number of reads required to be seen at a site for it to be analyzed by CRISPResso', default=50)
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
    u_group.add_argument('--suppress_dedup_input_based_on_aln_pos_and_UMI', help='Suppress deduplication based on alignment position and UMI. If not set, input reads will be deduplicated based on alignment position and UMI', action='store_true')
    u_group.add_argument('--umi_regex', type=str, help='String specifying regex that UMI must match', default='NNWNNWNNN')

    #R1/R2 support settings
    r_group = parser.add_argument_group('R1/R2 support settings')
    r_group.add_argument('--r1_r2_support_max_distance', type=int, help='Max distance between r1 and r2 for the read pair to be classified as "supported" by r2', default=10000)
    r_group.add_argument('--discard_reads_without_r2_support', help='If set, reads without r2 support will be discarded from final analysis and counts', action='store_true')


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
    if cmd_args.root is None:
        if len(cmd_args.settings_file) > 0:
            settings['root'] = cmd_args.settings_file[0] + ".CRISPRlungo"
        else:
            settings['root'] = "CRISPRlungo"

    settings['debug'] = cmd_args.debug
    if 'debug' in settings_file_args:
        settings['debug'] = (settings_file_args['debug'].lower() == 'true')
        settings_file_args.pop('debug')

    logger = logging.getLogger('CRISPRlungo')
    logging_level = logging.INFO
    if settings['debug']:
        logging_level=logging.DEBUG

    logger.setLevel(logging.DEBUG)

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

    settings['cuts'] = cmd_args.cuts
    if 'cuts' in settings_file_args:
        settings['cuts'] = settings_file_args['cuts'].split(" ")
        settings_file_args.pop('cuts')
    if 'cut_sites' in settings_file_args:
        settings['cuts'].extend(settings_file_args['cut_sites'].split(" "))
        settings_file_args.pop('cut_sites')
    for cut in settings['cuts']:
        if ":" not in cut:
            parser.print_usage()
            raise Exception('Error: cut specification %s is in the incorrect format (must be given as chr1:234',cut)

    settings['on_target_cut_sites'] = cmd_args.on_target_cut_sites
    if 'on_target_cut_sites' in settings_file_args:
        settings['on_target_cut_sites'] = settings_file_args['on_target_cut_sites'].split(" ")
        settings_file_args.pop('on_target_cut_sites')

    settings['cut_classification_annotations'] = cmd_args.cut_classification_annotations
    if 'cut_classification_annotations' in settings_file_args:
        settings['cut_classification_annotations'] = settings_file_args['cut_classification_annotations'].split(" ")
        for idx,arg in enumerate(settings['cut_classification_annotations']):
            arg_els = arg.split(":")
            if len(arg_els) != 4: # chr:pos:direction:annotation
                parser.print_usage()
                raise Exception('Error: cut classification annotations must contain 4 elements (e.g. chr1:234:left:my_annotation). Please check the value "' + arg + '"')
            settings['cut_classification_annotations'][idx] = arg.replace("__"," ") #replace __ with space
        settings_file_args.pop('cut_classification_annotations')

    settings['PAM'] = cmd_args.PAM
    if 'PAM' in settings_file_args:
        settings['PAM'] = settings_file_args['PAM']
        settings_file_args.pop('PAM')

    settings['guide_sequences'] = cmd_args.guide_sequences
    if 'guide_sequences' in settings_file_args:
        settings['guide_sequences'] = settings_file_args['guide_sequences'].split(" ")
        settings_file_args.pop('guide_sequences')

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

    settings['min_primer_aln_score'] = cmd_args.min_primer_aln_score
    if 'min_primer_aln_score' in settings_file_args:
        settings['min_primer_aln_score'] = int(settings_file_args['min_primer_aln_score'])
        settings_file_args.pop('min_primer_aln_score')

    settings['min_primer_length'] = cmd_args.min_primer_length
    if 'min_primer_length' in settings_file_args:
        settings['min_primer_length'] = int(settings_file_args['min_primer_length'])
        settings_file_args.pop('min_primer_length')

    settings['primer_filter_use_cutadapt'] = cmd_args.primer_filter_use_cutadapt
    if 'primer_filter_use_cutadapt' in settings_file_args:
        settings['primer_filter_use_cutadapt'] = (settings_file_args['primer_filter_use_cutadapt'].lower() == 'true')
        settings_file_args.pop('primer_filter_use_cutadapt')

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

    settings['ignore_n'] = cmd_args.ignore_n
    if 'ignore_n' in settings_file_args:
        settings['ignore_n'] = (settings_file_args['ignore_n'].lower() == 'true')
        settings_file_args.pop('ignore_n')

    settings['run_crispresso_on_novel_sites'] = cmd_args.run_crispresso_on_novel_sites
    if 'run_crispresso_on_novel_sites' in settings_file_args:
        settings['run_crispresso_on_novel_sites'] = (settings_file_args['run_crispresso_on_novel_sites'].lower() == 'true')
        settings_file_args.pop('run_crispresso_on_novel_sites')

    settings['crispresso_min_count'] = cmd_args.crispresso_min_count
    if 'crispresso_min_count' in settings_file_args:
        settings['crispresso_min_count'] = int(settings_file_args['crispresso_min_count'])
        settings_file_args.pop('crispresso_min_count')

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
        settings['n_processes'] = int(settings['n_processes'])

    settings['novel_cut_merge_distance'] = cmd_args.novel_cut_merge_distance
    if 'novel_cut_merge_distance' in settings_file_args:
        settings['novel_cut_merge_distance'] = int(settings_file_args['novel_cut_merge_distance'])
        settings_file_args.pop('novel_cut_merge_distance')

    settings['homology_cut_merge_distance'] = cmd_args.homology_cut_merge_distance
    if 'homology_cut_merge_distance' in settings_file_args:
        settings['homology_cut_merge_distance'] = int(settings_file_args['homology_cut_merge_distance'])
        settings_file_args.pop('homology_cut_merge_distance')

    settings['dedup_input_on_UMI'] = cmd_args.dedup_input_on_UMI
    if 'dedup_input_on_UMI' in settings_file_args:
        settings['dedup_input_on_UMI'] = (settings_file_args['dedup_input_on_UMI'].lower() == 'true')
        settings_file_args.pop('dedup_input_on_UMI')

    settings['suppress_dedup_input_based_on_aln_pos_and_UMI'] = cmd_args.suppress_dedup_input_based_on_aln_pos_and_UMI
    if 'suppress_dedup_input_based_on_aln_pos_and_UMI' in settings_file_args:
        settings['suppress_dedup_input_based_on_aln_pos_and_UMI'] = (settings_file_args['suppress_dedup_input_based_on_aln_pos_and_UMI'].lower() == 'true')
        settings_file_args.pop('suppress_dedup_input_based_on_aln_pos_and_UMI')

    settings['umi_regex'] = cmd_args.umi_regex
    if 'umi_regex' in settings_file_args:
        settings['umi_regex'] = settings_file_args['umi_regex']
        settings_file_args.pop('umi_regex')

    settings['r1_r2_support_max_distance'] = cmd_args.r1_r2_support_max_distance
    if 'r1_r2_support_max_distance' in settings_file_args:
        settings['r1_r2_support_max_distance'] = int(settings_file_args['r1_r2_support_max_distance'])
        settings_file_args.pop('r1_r2_support_max_distance')

    settings['discard_reads_without_r2_support'] = cmd_args.discard_reads_without_r2_support
    if 'discard_reads_without_r2_support' in settings_file_args:
        settings['discard_reads_without_r2_support'] = (settings_file_args['discard_reads_without_r2_support'].lower() == 'true')
        settings_file_args.pop('discard_reads_without_r2_support')

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
        settings['fastq_r2'] = settings_file_args['fastq_r2']
        settings_file_args.pop('fastq_r2')
    if settings['fastq_r2'] and not os.path.isfile(settings['fastq_r2']):
        raise Exception('Error: fastq_r2 file %s does not exist',settings['fastq_r2'])

    settings['fastq_umi'] = cmd_args.fastq_umi
    if 'fastq_umi' in settings_file_args:
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

    return settings, logger


def assert_dependencies(cutadapt_command='cutadapt',samtools_command='samtools',bowtie2_command='bowtie2',crispresso_command='CRISPResso',casoffinder_command='cas-offinder'):
    """
    Asserts the presence of required software (faidx, bowtie2, casoffinder)

    params:
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

    param:
        fastq: read1 file
        number_reads_to_check: the number of reads to read in

    returns:
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

    param:
        fastq: fastq file

    returns:
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

    params:
        root: root for written files
        genome: location of genome to use
        pam: PAM sequence of guide
        guides: list of sequences of on-target guides (not including PAM)
        cleavage_offset: position where cleavage occurs (relative to end of spacer seq -- for Cas9 this is -3)
        num_mismatches: number of mismatches to find
        casoffinder_command: location of casoffinder to run
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch
    returns:
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
    linked_genome = root + '.genome.fa'
    if os.path.exists(linked_genome):
        os.remove(linked_genome)
    os.symlink(genome,linked_genome)


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

def prep_input(root, primer_seq, guide_seqs, cleavage_offset, fastq_r1, samtools_command, genome, bowtie2_command, bowtie2_genome, can_use_previous_analysis=False):
    """
    Prepares primer info by identifying genomic location if possible and calculates statistics by analyzing input fastq_r1 file
    Prepares cut info by aligning guides to the genome
    Generally, primer_seq is genomic, and is used to amplify a region near the on-target cut site
    In cases where an exogenous sequence is amplified (e.g. GUIDE-seq) primer_seq may not be genomic
    The origin sequence is determined, and trimmed from input reads

    params:
        root: root for written files
        primer_seq: sequence of primer sequence.
        guide_seqs: sequences of guides used in experiment
        cleavage_offset: offset for guide where cut occurs
        fastq_r1: input fastq
        samtools_command: location of samtools to run
        genome: location of genome to extract sequence from
        bowtie2_command: location of bowtie2 to run
        bowtie2_genome: bowtie2-indexed genome to align to
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch


    returns:
        origin_seq: common amplified region at primer to be removed from input sequences
        cut_sites: list of cut locations for input guides
        cut_annotations: dict of cut_chr:cut_site->annotation for description of cut [type, origin_status, guide_direction] 
            (type=On-target, Off-target, Known, Casoffinder, etc) (origin_status=Not-origin or Origin:left or Origin:right) (guide_direction='FW' or 'RC' for which direction the guide binds)
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
            cut_annotations[key] = ['On-target','Not-origin',guide_is_rc_str,guide_seq]
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
            cut_annotations[key] = ['On-target','Not-origin',guide_is_rc_str,guide_seq]
            if primer_is_genomic and primer_chr == guide_chr:
                if closest_cut_site_loc == -1:
                    closest_cut_site_dist = abs(primer_loc - cut_loc)
                    closest_cut_site_loc = cut_loc
                elif abs(primer_loc - cut_loc) < closest_cut_site_dist:
                    closest_cut_site_dist = abs(primer_loc - cut_loc)
                    closest_cut_site_loc = cut_loc
        else:
            raise Exception('Could not find unique genomic coordinates for guide %s'%guide_seq)

    origin_seq = primer_seq
    if primer_is_genomic:
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
        if primer_is_rc:
            origin_seq = reverse_complement(origin_seq)
        logger.debug('Got origin sequence: ' + origin_seq)

    logger.debug('Getting read length and number of total reads')
    av_read_length = get_av_read_len(fastq_r1)
    num_reads_input = get_num_reads_fastq(fastq_r1)

    cut_site_count = len(cut_sites)
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

    params:
        root: root for written files
        fastq_r1: R1 reads to dedup
        fastq_r2: R2 reads to dedup
        umi_regex: string specifying regex that UMI must match
        min_umi_seen_to_keep_read: min number of times a umi must be seen to keep the umi and read (e.g. if set to 2, a read-UMI pair that is only seen once will be discarded)
        write_UMI_counts: if True, writes a file with the UMI counts
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch
    returns:
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

    params:
        root: root for written files
        fastq_r1: R1 reads to dedup
        fastq_r2: R2 reads to dedup
        fastq_umi: UMI fastq to dedup on
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch
    returns:
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
            if len(line_els) > 3:
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

def filter_on_primer(root,fastq_r1,fastq_r2,origin_seq,min_primer_aln_score,min_primer_length=10,min_read_length=30,transposase_adapter_seq='CTGTCTCTTATACACATCTGACGCTGCCGACGA',n_processes=1,cutadapt_command='cutadapt',use_cutadapt=False,keep_intermediate=False,can_use_previous_analysis=False):
    """
    Trims the primer from the input reads, only keeps reads with the primer present in R1
    Also trims the transposase adapter sequence from reads and filters reads that are too short

    params:
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
        use_cutadapt: if true, use cutadapt to identify primer sequence. If false, custom lookhead search will be used.
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
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

                logger.info('Using %d previously-filtered fastq sequences trimmed for primers'%post_trim_read_count)

                return(filtered_on_primer_fastq_r1,filtered_on_primer_fastq_r2,post_trim_read_count,filter_on_primer_plot_obj)

            else:
                logger.info('Could not recover previously-filtered results. Reprocessing.')

    logger.info('Filtering and trimming primers from reads')

    post_trim_read_count = 'NA'
    too_short_read_count = 'NA'
    untrimmed_read_count = 'NA'
    if use_cutadapt:
        cutadapt_log = root + ".cutadapt.log"
        if fastq_r2 is None:
            filtered_on_primer_fastq_r1 = root + ".trimmed.fq.gz"
            filtered_on_primer_fastq_r2 = None
            trim_command = "%s -g primerSeq=%s -o %s --rename='{id} {comment} primer={match_sequence}' --discard-untrimmed --minimum-length %d --cores %s %s > %s"%(cutadapt_command,origin_seq,filtered_on_primer_fastq_r1,min_read_length,n_processes,fastq_r1,cutadapt_log)
        else:
            filtered_on_primer_fastq_r1 = root + ".r1.has_primer.fq.gz"
            filtered_on_primer_fastq_r2 = root + ".r2.has_primer.fq.gz"
            origin_seq_rc = reverse_complement(origin_seq)
            trim_command = "%s -g primerSeq=%s -A %s -o %s -p %s --rename='{id} {comment} primer={match_sequence}' --discard-untrimmed --minimum-length %d --pair-filter=first --cores %s %s %s > %s"%(cutadapt_command,origin_seq,origin_seq_rc,filtered_on_primer_fastq_r1,filtered_on_primer_fastq_r2,min_read_length,n_processes,fastq_r1,fastq_r2,cutadapt_log)

        logger.debug(trim_command)
        trim_result = subprocess.check_output(trim_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)

        if not os.path.exists(cutadapt_log):
            logger.error('Error found while running command:\n'+trim_command+"\nOutput: "+trim_result)

        with open(cutadapt_log,'r') as fin:
            for line in fin:
                m = re.search(r'written \(passing filters\):\s+([\d,]+) \(\d.*%\)',line)
                if m:
                    post_trim_read_count = int(m.group(1).replace(',',''))

                m = re.search(r'that were too short:\s+([\d,]+) \(\d.*%\)',line)
                if m:
                    too_short_read_count = int(m.group(1).replace(',',''))

                m = re.search(r'discarded as untrimmed:\s+([\d,]+) \(\d.*%\)',line)
                if m:
                    untrimmed_read_count = int(m.group(1).replace(',',''))

        if post_trim_read_count == 'NA':
            logger.error('Could not parse trim read count from file ' + cutadapt_log)

    else: # if not use cutadapt for trimming reads
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
            logger.error('Could not parse trim read count from file ' + cutadapt_log)
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
    plot_label = "Reads were filtered to contain primer origin sequence <pre>" + origin_seq +"</pre>"
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

def align_reads(root,fastq_r1,fastq_r2,bowtie2_reference,bowtie2_command='bowtie2',bowtie2_threads=1,samtools_command='samtools',keep_intermediate=False,can_use_previous_analysis=False):
    """
    Aligns reads to the provided reference

    params:
        root: root for written files
        fastq_r1: fastq_r1 to align
        fastq_r2: fastq_r2 to align
        bowtie2_reference: bowtie2 reference to align to (either artificial targets or reference)
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        samtools_command: location of samtools to run
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        mapped_bam_file: aligned reads in bam format
        genome_aligned_reads_count: number of aligned reads using Bowtie2
    """

    logger = logging.getLogger('CRISPRlungo')

    mapped_bam_file = root + ".bam"

    info_file = root + '.info'
    if os.path.isfile(info_file) and can_use_previous_analysis:
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 1:
                mapped_bam_file = line_els[0]
                logger.info('Using previously-performed alignment')
                return mapped_bam_file
        logger.info('Could not recover previously-performed alignment. Reanalyzing.')


    bowtie_log = root + '.bowtie2Log'
    if fastq_r2 is not None: #paired-end reads
        logger.info('Aligning paired reads using %s'%(bowtie2_command))
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --very-sensitive-local --soft-clipped-unmapped-tlen --threads {bowtie2_threads} -x {bowtie2_reference} -1 {fastq_r1} -2 {fastq_r2} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
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
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --very-sensitive-local --soft-clipped-unmapped-tlen --threads {bowtie2_threads} -x {bowtie2_reference} -U {fastq_r1} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
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
    cut_merge_dist=100,collapse_to_homology_dist=10000,guide_homology_max_gaps=2,guide_homology_max_mismatches=5,
    arm_min_matched_start_bases=10,genome_map_resolution=1000000,suppress_dedup_input_based_on_aln_pos_and_UMI=False,
    discard_reads_without_r2_support=False,samtools_command='samtools',keep_intermediate=False,can_use_previous_analysis=False):
    """
    Makes final read assignments (after deduplicating based on UMI and alignment location)

    params:
        root: root for written files
        genome_mapped_bam: bam of reads aligned to genome
        origin_seq: common amplified region at primer (to first cut site)
        cut_sites: array of known/input cut sites
        cut_annotations: dict of cut_chr:cut_site->annotation for description of cut [type, origin_status, guide_direction] 
            (type=On-target, Off-target, Known, Casoffinder, etc) (origin_status=Not-origin or Origin:left or Origin:right) (guide_direction='FW' or 'RC' for which direction the guide binds)
        cut_classification_annotations: dict of cut_chr:cut_site:direction -> annotation (as provided by users as parameters)
        guide_seqs: sequences of guides used in experiment (required here for finding homology at novel cut sites)
        cleavage_offset: position where cleavage occurs (relative to end of spacer seq -- for Cas9 this is -3)
        min_primer_length: minimum length of sequence that matches between the primer/origin sequence and the read sequence
        genome: path to genome fa file (required here for finding genome sequences for guide homology)
        r1_r2_support_max_distance: max distance between r1 and r2 for the read pair to be classified as 'supported' by r2
        cut_merge_dist: cuts within this distance (bp) will be merged
        collapse_to_homology_dist: novel cuts within this distance to a site with homology to a given guide. Novel cuts are first merged to the nearest site with homology (within collapse_to_homology_dist), and remaining cuts are merged with eachother (specified by the 'cut_merge_dist' parameter).
        guide_homology_max_gaps: when searching a sequence for homology to a guide sequence, it will be homologous if his has this number of gaps or fewer in the alignment of the guide to a sequence
        guide_homology_max_mismatches: when searching a sequence for homology to a guide sequence, it will be homologous if his has this number of mismatches or fewer in the alignment of the guide to a sequence
        arm_min_matched_start_bases: number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each 'side' a read.
        genome_map_resolution: window size (bp) for reporting number of reads aligned
        dedup_based_on_aln_pos: if true, reads with the same UMI and alignment positions will be removed from analysis
        discard_reads_without_r2_support: if true, reads without r2 support will be discarded from final analysis and counts
        samtools_command: location of samtools to run
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        final_assignment_filename: filename of final assignments for each read
        final_read_ids_for_cuts: dict of readID=>cut
        final_cut_counts: dict of cutID=>count how many times each cut was seen -- when we iterate through the reads (fastq) in the next step, we check to see that the cut was seen above a threshold before printing those reads. The count is stored in this dict.
        cut_classification_lookup: dict of cutID=> type (e.g. Linear, Translocation, etc)
        chr_aln_plot_obj: plot showing alignment locations of reads
        tlen_plot_object: plot showing distribution of template lengths
        deduplication_plot_obj: plot showing how many reads were deuplicated using UMI + location
        classification_plot_obj: plot showing assignments (linear, translocation, etc)
        classification_indel_plot_obj: plot showing assignments (linear, translocation, etc) with indel classification (short indels)
        tx_order_plot_obj: plot showing ordered translocations and counts
        tx_count_plot_obj: plot showing translocations and counts
        origin_indel_hist_plot_obj: plot of indel lenghts for origin
        origin_inversion_hist_plot_obj: plot of inversion lengths at origin (reads pointing opposite directions)
        origin_deletion_hist_plot_obj: plot of deletion lengths for ontarget (reads pointing same direction/away from the primer)
        origin_indel_depth_plot_obj: plot of read depth around origin
        r1_r2_support_plot_obj: plot of how many reads for which r1 and r2 support each other
        r1_r2_support_dist_plot_obj: plot of how far apart the r1/r2 supporting reads were aligned from each other
        r1_r2_no_support_dist_plot_obj: plot of how far apart the r1/r2 reads that did NOT support each other (but were on the same chromosome)

    """
    logger = logging.getLogger('CRISPRlungo')

    final_file = root + '.final_assignments'
    info_file = root + '.info'
    #check to see if this anaylsis has been completed previously
    if os.path.isfile(info_file) and os.path.isfile(final_file) and can_use_previous_analysis:
        read_total_read_count = 0 #how many lines we read in this time by reading the final assignments file
        previous_total_read_count = -1 #how many lines were reported to be read the preivous time (if the analysis was completed previously)
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 7:
                (total_reads_processed,dups_removed,bad_alignments_removed,nonsupporting_alignments_removed,final_total_count,found_cut_point_total,final_cut_point_total) = [int(x) for x in line_els]
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

                r1_r2_no_support_dist_plot_obj_str = fin.readline().rstrip('\n')
                r1_r2_no_support_dist_plot_obj = None
                if r1_r2_no_support_dist_plot_obj_str != "" and r1_r2_no_support_dist_plot_obj_str != "None":
                    r1_r2_no_support_dist_plot_obj = PlotObject.from_json(r1_r2_no_support_dist_plot_obj_str)


        if previous_total_read_count > -1:
            #read final assignments from file
            final_read_ids_for_cuts = defaultdict(int)
            final_cut_counts = defaultdict(int)
            with open (final_file,'r') as fin:
                head_line = fin.readline().rstrip('\n')
                head_line_els = head_line.split("\t")
                read_total_read_count = 0
                final_cut_ind = 10
                if head_line_els[final_cut_ind] == "final_cut_pos":
                    logger.info('Reading previously-processed assignments')
                    #if header matches, read file
                    for line in fin:
                        read_total_read_count += 1
                        line_els = line.rstrip("\r\n").split("\t")
                        final_read_ids_for_cuts[line_els[0]] = line_els[final_cut_ind]
                        final_cut_counts[line_els[final_cut_ind]] += 1

                if final_total_count == read_total_read_count:
                    logger.info('Using previously-processed assignments for ' + str(read_total_read_count) + ' total reads')
                    return (final_file,final_read_ids_for_cuts,final_cut_counts,cut_classification_lookup,
                        chr_aln_plot_obj,tlen_plot_obj,deduplication_plot_obj,tx_order_plot_obj,tx_count_plot_obj,classification_plot_obj,classification_indel_plot_obj,
                        origin_indel_hist_plot_obj,origin_inversion_hist_plot_obj,origin_deletion_hist_plot_obj,origin_indel_depth_plot_obj,
                        r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,r1_r2_no_support_dist_plot_obj)
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
    seen_reads_for_dedup = {} # put read info into this dict for deduplicating
    total_reads_processed = 0 # includes pairs
    total_r1_processed = 0
    dups_removed = 0
    bad_alignments_removed = 0
    nonsupporting_alignments_removed = 0
    cut_points_by_chr = defaultdict(lambda: defaultdict(int)) # dict of cut points for each chr contianing observedcuts cut_points_by_chr[chr1][1559988] = count seen
    aligned_tlens = defaultdict(int)
    aligned_chr_counts = defaultdict(int)
    aln_pos_by_chr = defaultdict(lambda: defaultdict(int))

    #keep track of distance between R1/R2
    r1_r2_support_distances = defaultdict(int) #for reads that are separated by less than r1_r2_support_max_distance
    r1_r2_not_support_distances = defaultdict(int) #for reads that are separated by more than r1_r2_support_max_distance

    #counts for 'support' status
    r1_r2_support_status_counts = defaultdict(int)
    final_file_tmp = root + '.final_assignments.tmp'
    with open(final_file_tmp,'w') as af1:
        af1.write("\t".join([str(x) for x in ['read_id','curr_position','cut_pos','cut_direction','del_primer','ins_target','insert_size','r1_r2_support_status','r1_r2_support_dist','r1_r2_orientation_str']])+"\n")

        if not os.path.exists(genome_mapped_bam):
            raise Exception('Cannot find bam file at ' + genome_mapped_bam)
        for line in read_command_output('%s view -F 256 %s'%(samtools_command,genome_mapped_bam)):
            if line.strip() == "": break

            line_els = line.split("\t")
            total_reads_processed += 1

            read_has_multiple_segments = int(line_els[1]) & 0x1
            read_is_paired_read_1 = int(line_els[1]) & 0x40 # only set with bowtie if read is from paired reads

            if read_has_multiple_segments and not read_is_paired_read_1:
                continue
            total_r1_processed += 1

            (line_info, primer_info) = line_els[0].split(" primer=")

            barcode = line_info.split(":")[-1]
            aln1 = line_els[2] + ":" + line_els[3]
            aln2 = line_els[6] + ":" + line_els[7] # if unpaired, this will just be *:0, so we can still use it in our key
            tlen = abs(int(line_els[8]))
            key = "%s %s %s %s"%(barcode,aln1,aln2,tlen)
            if key in seen_reads_for_dedup:
                seen_reads_for_dedup[key] += 1
                if not suppress_dedup_input_based_on_aln_pos_and_UMI:
                    dups_removed += 1
                    continue
            else:
                seen_reads_for_dedup[key] = 1

            line_unmapped = int(line_els[1]) & 0x4

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

            if line_unmapped or \
                    (start_clipped > 0 and end_clipped > 0) or \
                    left_matches < arm_min_matched_start_bases or \
                    right_matches < arm_min_matched_start_bases:

                bad_alignments_removed += 1
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

            del_primer = len(origin_seq) - len(primer_info) #how many bases are deleted from the full primer
            ins_target = start_clipped # how many bases were clipped from alignment (insertions at target)

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
                if discard_reads_without_r2_support and r1_r2_support_status != r1_r2_support_str_supported:
                    nonsupporting_alignments_removed += 1
                    continue
            #done setting r1/r2 support status for paired reads

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

    chr_aln_plot_root = root + ".chr_alignments"
    keys = sorted(aligned_chr_counts.keys())
    vals = [aligned_chr_counts[key] for key in keys]
    with open(chr_aln_plot_root+".txt","w") as chrs:
        chrs.write('chr\tnumReads\n')
        for key in sorted(aligned_chr_counts.keys()):
            chrs.write(key + '\t' + str(aligned_chr_counts[key]) + '\n')

    fig = plt.figure(figsize=(12,12))
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

        fig = plt.figure(figsize=(12,12))
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
    deduplication_plot_obj = None
    if total_reads_processed > 0:
        labels = ['Not duplicate','Duplicate']
        values = [total_reads_processed-dups_removed,dups_removed]
        deduplication_plot_obj_root = root + ".deduplication_by_UMI_and_aln_pos"
        fig = plt.figure(figsize=(12,12))
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
        dup_counts = defaultdict(int)
        dup_keys = seen_reads_for_dedup.keys()
        for dup_key in dup_keys:
            dup_counts[seen_reads_for_dedup[dup_key]] += 1
        keys = sorted(dup_counts.keys())
        vals = [dup_counts[key] for key in keys]
        with open(deduplication_plot_obj_root+".txt","w") as summary:
            summary.write('Reads per UMI\tNumber of UMIs\n')
            for key in keys:
                summary.write(str(key) + '\t' + str(dup_counts[key]) + '\n')
        if len(vals) > 0:
            ax2.bar(keys,vals)
            ax2.set_ymargin(0.05)
        else:
            ax2.bar(0,0)
        ax2.set_ylabel('Number of UMIs')
        ax2.set_title('Reads per UMI')



        plt.savefig(deduplication_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(deduplication_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        dedup_note = ""
        if suppress_dedup_input_based_on_aln_pos_and_UMI:
            dedup_note = ". Note that deduplication of reads for final counts has been turned off by with the parameter --suppress_dedup_input_based_on_aln_pos_and_UMI."

        deduplication_plot_obj = PlotObject(
                plot_name = deduplication_plot_obj_root,
                plot_title = 'Duplicate read counts',
                plot_label = 'Number of reads that were duplicates based on alignment and UMI' + dedup_note,
                plot_datas = [('Read assignments',final_file)]
                )

    #now that we've aggregated all cut sites, find likely cut positions and record assignments to those positions in final_cut_point_lookup
    #also annotate the origin cut
    origin_chr = None
    origin_cut_pos = None
    origin_direction = None
    #cut sites are input by user (or found by casoffinder)
    known_cut_points_by_chr = {}
    for cut_site in cut_sites:
        (cut_chr,cut_pos) = cut_site.split(":")
        if cut_site in cut_annotations:
            origin_status = cut_annotations[cut_site][1]
            if origin_status == 'Origin:right':
                origin_chr = cut_chr
                origin_cut_pos = int(cut_pos)
                origin_direction = 'right'
            elif origin_status == 'Origin:left':
                origin_chr = cut_chr
                origin_cut_pos = int(cut_pos)
                origin_direction = 'left'

        if cut_chr not in known_cut_points_by_chr:
            known_cut_points_by_chr[cut_chr] = []
        known_cut_points_by_chr[cut_chr].append(int(cut_pos))
        #add known cut points to final cut point lookup as well
        if cut_chr not in cut_points_by_chr:
            cut_points_by_chr[cut_chr] = defaultdict(int)
        cut_points_by_chr[cut_chr][int(cut_pos)] += 1

    #the counts in this file include only reads that weren't marked as duplicates
    final_cut_points_by_chr = defaultdict(lambda: defaultdict(int)) #dict of final cut points by chromosome, contains count of reads with that cut point
    final_cut_point_lookup = {} #dict from old(fuzzy/imprecise) position to new

    cut_point_homology_info = {} #dict containing cut points with sufficient homology
    with open(root+".cut_homology",'w') as fout:
        fout.write("\t".join([str(x) for x in ['guide_seq','cut_chr','cut_point','is_valid_homology_site','potential_guide','best_score','match_direction','n_matches','n_mismatches','n_gaps','aln_guide','aln_ref']])+"\n")
        aln_match_score = 5
        aln_gap_extend_score = -1
        aln_gap_open_score = -5
        aln_matrix = CRISPResso2Align.make_matrix(match_score=aln_match_score)
        padded_seq_len = 25 #how many bp to extend around cut to search for guide
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

    #create dict of cut_point > closest cut point with homology
    closest_cut_point_with_homology_lookup = {}

    #first create dict with known sites and sites with homology
    cut_points_with_homology_by_chr = defaultdict(list)
    #add known sites
    for cut_chr in known_cut_points_by_chr:
        cut_points_with_homology_by_chr[cut_chr] = known_cut_points_by_chr[cut_chr][:]
    #split cut points with homology by chr to allow sorting
    for cut_point in cut_point_homology_info:
        cut_chr,cut_pos = cut_point.split(":")
        if int(cut_pos) not in cut_points_with_homology_by_chr[cut_chr]:
            cut_points_with_homology_by_chr[cut_chr].append(int(cut_pos))

    #assign each cut point to the closest site with homology if possible
    for cut_chr in cut_points_by_chr:
        these_cut_points = sorted(cut_points_by_chr[cut_chr])

        these_cut_points_with_homology = sorted(cut_points_with_homology_by_chr[cut_chr])
        #if no cut points on this chr
        if len(these_cut_points_with_homology) == 0:
            for cut_point in these_cut_points:
                closest_cut_point_with_homology_lookup[cut_chr+":"+str(cut_point)] = None
        else:
            curr_ind = 0
            for cut_point in these_cut_points:
                dist_to_curr = abs(cut_point - these_cut_points_with_homology[curr_ind])
                #if curr_ind is at last item, choose it
                if curr_ind + 1 >= len(these_cut_points_with_homology):
                    if dist_to_curr <= collapse_to_homology_dist:
                        closest_cut_point_with_homology_lookup[cut_chr+":"+str(cut_point)] = cut_chr+":"+str(these_cut_points_with_homology[curr_ind])
                    else:
                        closest_cut_point_with_homology_lookup[cut_chr+":"+str(cut_point)] = None
                #otherwise choose closest of curr or next
                else:
                    dist_to_next = abs(cut_point - these_cut_points_with_homology[curr_ind+1])
                    if dist_to_curr <= dist_to_next:
                        if dist_to_curr <= collapse_to_homology_dist:
                            closest_cut_point_with_homology_lookup[cut_chr+":"+str(cut_point)] = cut_chr+":"+str(these_cut_points_with_homology[curr_ind])
                        else:
                            closest_cut_point_with_homology_lookup[cut_chr+":"+str(cut_point)] = None
                    else:
                        curr_ind += 1
                        if dist_to_next <= collapse_to_homology_dist:
                            closest_cut_point_with_homology_lookup[cut_chr+":"+str(cut_point)] = cut_chr+":"+str(these_cut_points_with_homology[curr_ind])
                        else:
                            closest_cut_point_with_homology_lookup[cut_chr+":"+str(cut_point)] = None

    found_cut_point_total = 0
    closest_cut_point_lookup = {}
    for cut_chr in cut_points_by_chr:
        found_cut_point_total += len(cut_points_by_chr[cut_chr])
        these_cut_points = sorted(cut_points_by_chr[cut_chr])
        last_seen_point = these_cut_points[0]
        curr_points = [last_seen_point]
        for i in range(1,len(these_cut_points)+1):
            #if next point is beyond cut_merge_dist, merge and reset!
            #do weighted sum
            this_sum = 0
            this_tot = 0
            for curr_point in curr_points:
                this_sum += curr_point*cut_points_by_chr[cut_chr][curr_point]
                this_tot += cut_points_by_chr[cut_chr][curr_point]
            this_weighted_mean = this_sum/float(this_tot)
            #this_mean = int(sum(curr_points)/float(len(curr_points)))
            if i == len(these_cut_points) or abs(these_cut_points[i] - this_weighted_mean) > cut_merge_dist:
                this_pos = int(this_weighted_mean)
                min_dist = cut_merge_dist
                #once we've selected a putative cut point, make sure it isn't within min_dist from known cut points
                #and select the closest point (keep shrinking min_dist in case another one is found closer)
                if cut_chr in known_cut_points_by_chr:
                    for known_cut_point in known_cut_points_by_chr[cut_chr]:
                        this_dist = abs(known_cut_point - this_weighted_mean)
                        if this_dist <= min_dist:
                            this_pos = known_cut_point
                            min_dist = this_dist

                #add this assignment to the lookup
                for curr_point in curr_points:
                    closest_cut_point_lookup[cut_chr+":"+str(curr_point)] = cut_chr+":"+str(this_pos)

                #reset current points
                curr_points = []

            if i < len(these_cut_points):
                curr_points.append(these_cut_points[i])

    #make final assignments - either to the next homologous site (if available), or to the next closest site
    all_cut_points_count = 0
    all_cut_points_assigned_to_homology_count = 0 #how many cut sites were assigned to their nearest homologous site
    all_cut_points_collapsed_count = 0 #how many cut sites were collapsed (not assigned to homologous sites)
    for cut_chr in cut_points_by_chr:
        for cut_point in cut_points_by_chr[cut_chr]:
            all_cut_points_count += 1
            cut_key = cut_chr+":"+str(cut_point)
            if closest_cut_point_with_homology_lookup[cut_key] is not None:
                cut_assignment = closest_cut_point_with_homology_lookup[cut_key]
                all_cut_points_assigned_to_homology_count += 1
            else:
                cut_assignment = closest_cut_point_lookup[cut_key]
                all_cut_points_collapsed_count += 1

            cut_assignment_pos = int(cut_assignment.split(":")[1])
            final_cut_point_lookup[cut_key] = cut_assignment
            final_cut_points_by_chr[cut_chr][cut_assignment_pos] += cut_points_by_chr[cut_chr][cut_point]

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
                cut_classification_lookup['%s:%s:%s'%(chrom,pos,'left')] = this_anno[0]+' off-target'
                cut_classification_lookup['%s:%s:%s'%(chrom,pos,'right')] = this_anno[0] + ' off-target'
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

    logger.info('Processed ' + str(all_cut_points_count) + ' observed cut points')
    logger.info(str(all_cut_points_assigned_to_homology_count) + ' cut points were collapsed to known/given/homologous cut sites')
    logger.info(str(all_cut_points_collapsed_count) + ' cut points were collapsed to nearby cut sites')
    logger.info(str(final_cut_point_total) + ' final cut points after collapsing')

    #add cut classifications from parameters
    for cut_classification_param in cut_classification_annotations:
        cut_classification_els = cut_classification_param.split(":")
        cut_key = ":".join(cut_classification_els[0:3])
        cut_val = cut_classification_els[3]
        cut_classification_lookup[cut_key] = cut_val

    with open(root+".cut_classification",'w') as fout:
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
    final_read_ids_for_cuts = defaultdict(int) # dict readID->cut
    final_cut_counts = defaultdict(int) # dict cut->count of reads

    if origin_direction == 'left':
        origin_direction_opposite = "right"
    else:
        origin_direction_opposite = "left"
    origin_indel_lengths = defaultdict(int) # dict indel_len -> count
    origin_inversion_lengths = defaultdict(int) # dict inversion_len -> count of reads on same chr as primer and same orientation with specified distance between origin and read (primer is to the genomic left of the cut site, this read is oriented to the genomic left side)
    origin_deletion_lengths = defaultdict(int) # dict deletion_len -> count of reads on same chr as primer and opposite direction with specified distance between origin and read (primer is to the genomic left of the cut site, this read is oriented to the genomic right side)
    origin_depth_counts_100bp = np.zeros(101) #100bp count the number of reads with deletions covering each position
    origin_depth_counts_100bp_primer = np.zeros(len(origin_seq)) #count deletions covering each base of primer (for those with deletions within 100bp)
    origin_depth_counts_100bp_total = 0 #count of total reads

    origin_depth_counts_2500bp = np.zeros(2501) #2500bp count of the number of reads with deletions covering that position
    origin_depth_counts_2500bp_primer = np.zeros(len(origin_seq)) #count reads covering each base of primer (for those with deletions within 2500bp)
    origin_depth_counts_2500bp_total = 0 #count of total reads

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
            final_read_ids_for_cuts[line_id] = final_cut_point_and_direction
            final_cut_counts[final_cut_point_and_direction] += 1


            (final_cut_point_chr,final_cut_point_start) = final_cut_point.split(":")
            final_cut_point_start = int(final_cut_point_start)
            aln_start = int(aln_loc.split(":")[1].split("-")[0])

            #read extends to the right of the cut point
            if cut_point_direction == "right":
                aln_indel = -1*(final_cut_point_start - aln_start)
            else: #read extends to the left of the cut point
                aln_indel = final_cut_point_start - aln_start - 1

            del_primer = int(line_els[del_primer_ind])

            final_cut_indel = ins_target - del_primer - aln_indel

            final_classification = cut_classification_lookup[final_cut_point_and_direction]
            final_classification_counts[final_classification] += 1

            indel_str = ''
            if final_cut_indel != 0:
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

                if aln_start >= origin_cut_pos:
                    this_del_len = aln_start - origin_cut_pos
                    if this_del_len <= 100:
                        origin_depth_counts_100bp_total += 1
                        origin_depth_counts_100bp[0:this_del_len+1] += 1
                        origin_depth_counts_100bp_primer[0:del_primer+1] += 1
                    if this_del_len <= 2500:
                        origin_depth_counts_2500bp_total += 1
                        origin_depth_counts_2500bp[0:this_del_len+1] += 1
                        origin_depth_counts_2500bp_primer[0:del_primer+1] += 1

            #finish iterating through read ids
        #close final file
    logger.info('Processed %d reads.'%(final_total_count))

    tx_keys = sorted(final_cut_counts.keys())
    tx_list_report = root + ".translocation_list.txt"
    with open (tx_list_report,'w') as fout:
        fout.write('cut\tcut_annotation\tcount\n')
        for tx_key in tx_keys:
            this_cut_annotation = 'Novel'
            if tx_key in cut_annotations:
                this_cut_annotation = cut_annotations[tx_key][0]
                fout.write("%s\t%s\t%s\n"%(tx_key,this_cut_annotation,final_cut_counts[tx_key]))
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
    top_sorted_tx_list = sorted_tx_list[:min(20,len(sorted_tx_list))]

    tx_order_plot_obj = None
    tx_count_plot_obj = None
    if len(top_sorted_tx_list) > 0:

        # make tx order plot
        if origin_chr is None:
            left_label = 'Exogenous'
            left_pos = 'None'
        else:
            left_label = origin_chr+':'+str(origin_cut_pos)+':'+origin_direction
            left_pos = origin_chr+':'+str(origin_cut_pos)
        pos_list = [':'.join(x.split(':')[0:2]) for x in top_sorted_tx_list]
        pos_list.append(left_pos)
        pos_list = sorted(list(set(pos_list)))

        cut_categories = ['On-target', 'Off-target', 'Known','Casoffinder','Novel']
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
        left_cols = []
        right_cols = []
        left_outline_cols = []
        right_outline_cols = []
        for tx_order_obj in top_sorted_tx_list:
            left_labs.append(cut_classification_lookup[tx_order_obj] + ": " + left_label)
            right_labs.append(tx_order_obj)
            counts.append(final_cut_counts[tx_order_obj])
            this_left_col = color_lookup[left_pos]
            left_cols.append(this_left_col)
            this_left_outline_col = 'None'
            left_outline_cols.append(this_left_outline_col)

            this_chr_pos = ':'.join(tx_order_obj.split(':')[0:2])
            this_right_col = color_lookup[this_chr_pos]
            this_cut_anno='Novel'
            if this_chr_pos in cut_annotations:
                this_cut_anno = cut_annotations[this_chr_pos][0]

            right_cols.append(this_right_col)
            this_right_outline_col = cut_category_lookup[this_cut_anno]
            right_outline_cols.append(this_right_outline_col)

        tx_plt = makeTxCountPlot(left_labs=left_labs, right_labs=right_labs, counts=counts, left_cols=left_cols, right_cols=right_cols, left_outline_cols=left_outline_cols, right_outline_cols=right_outline_cols,outline_cols=cut_colors.colors, outline_labs=cut_categories)
        plot_name = tx_list_report_root
        tx_plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
        tx_plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')

        tx_order_plot_obj = PlotObject(plot_name = plot_name,
                plot_title = 'Translocation Summary',
                plot_label = 'Final translocation counts',
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
        for loc_key in final_cut_counts:
            loc_key_els = loc_key.split(":")
            height_at_locs[loc_key_els[0]][int(loc_key_els[1])][loc_key_els[2]] = final_cut_counts[loc_key]
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

        tx_count_plot_obj = PlotObject(plot_name = plot_name,
                plot_title = 'Translocation Counts',
                plot_label = 'Final translocation counts',
                plot_datas = [('Fragment translocation list',tx_list_report)]
                )

        classification_labels = ['Linear'] #linear goes first in order
        for key in sorted(final_classification_counts.keys()):
            if key not in classification_labels and final_classification_counts[key] > 0:
                classification_labels.append(key)

        #assignment plot
        values = [final_classification_counts[x] for x in classification_labels]
        classification_plot_obj_root = root + ".classifications"
        with open(classification_plot_obj_root+".txt",'w') as summary:
            summary.write("\t".join(classification_labels)+"\n")
            summary.write("\t".join([str(x) for x in values])+"\n")
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(classification_labels)):
            pie_values.append(values[i])
            pie_labels.append(classification_labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Read classifications')
        plt.savefig(classification_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(classification_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(classification_labels,values)])
        classification_plot_obj = PlotObject(
                plot_name = classification_plot_obj_root,
                plot_title = 'Read Classification',
                plot_label = 'Read classification<br>'+plot_count_str,
                plot_datas = [
                    ('Alignment classifications',classification_plot_obj_root + ".txt"),
                    ('Read assignments',final_file)
                    ]
                )
        noindel_pie_labels = pie_labels #for use below
        noindel_pie_values = pie_values

        #assignment (with indels) plot

        short_indel_categories = ['',' short indels']
        classification_indel_labels = []
        classification_indel_counts = []
        for label in classification_labels:
            for short_indel_category in short_indel_categories:
                new_label = label + short_indel_category
                classification_indel_labels.append(new_label)
                classification_indel_counts.append(final_classification_indel_counts[new_label])

        classification_indel_plot_obj_root = root + ".classifications_with_indels"
        with open(classification_indel_plot_obj_root+".txt",'w') as summary:
            summary.write("\t".join(classification_indel_labels)+"\n")
            summary.write("\t".join([str(x) for x in classification_indel_counts])+"\n")
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        width = 0.1
        inner_wedge_properties = {"edgecolor":"w",'linewidth': 2}
        wedge_properties = {"width":width, "edgecolor":"w",'linewidth': 2}

        ax.pie(noindel_pie_values, labels=noindel_pie_labels,
            wedgeprops=inner_wedge_properties,autopct="%1.2f%%",radius=1-width, labeldistance=0.50,pctdistance=0.25)
        ax.pie(classification_indel_counts, labels=None,
            wedgeprops=wedge_properties,autopct="%1.2f%%",colors=['silver','crimson'],pctdistance=1)
        patch1 = Patch(color='silver', label='No indels')
        patch2 = Patch(color='crimson', label='Short indels')

        plt.legend(handles=[patch1, patch2])

        ax.set_title('Read classifications')
        plt.savefig(classification_indel_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(classification_indel_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(classification_indel_labels,classification_indel_counts)])
        classification_indel_plot_obj = PlotObject(
                plot_name = classification_indel_plot_obj_root,
                plot_title = 'Read Classification including short indels',
                plot_label = 'Read classification including annotations for presence of short indels<br>'+plot_count_str,
                plot_datas = [
                    ('Alignment classifications with indels',classification_indel_plot_obj_root + ".txt"),
                    ('Alignment classifications',classification_plot_obj_root + ".txt"),
                    ('Read assignments',final_file)
                    ]
                )

    origin_indel_hist_plot_obj = None # plot of indel lengths at origin
    if len(origin_indel_lengths) > 0:
        origin_indel_hist_plot_obj_root = root + ".origin_indel_histogram"
        keys = sorted(origin_indel_lengths.keys())
        vals = [origin_indel_lengths[key] for key in keys]
        with open(origin_indel_hist_plot_obj_root+".txt","w") as summary:
            summary.write('Indel_length\tnumReads\n')
            for key in keys:
                summary.write(str(key) + '\t' + str(origin_indel_lengths[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
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

        plot_label = 'Bar plot showing lengths of insertions (positive) or deletions (negative) for linear reads supporting the origin cut site at ' + origin_chr + ': ' + str(origin_cut_pos) 
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
    if len(origin_inversion_lengths) > 0:
        origin_inversion_hist_plot_obj_root = root + ".origin_inversion_histogram"
        keys = sorted(origin_inversion_lengths.keys())
        vals = [origin_inversion_lengths[key] for key in keys]
        with open(origin_inversion_hist_plot_obj_root+".txt","w") as summary:
            summary.write('Inversion_distance\tnumReads\n')
            for key in keys:
                summary.write(str(key) + '\t' + str(origin_inversion_lengths[key]) + '\n')

        plot_keys = [x for x in keys if abs(x) < origin_dist_cutoff]
        plot_vals = [origin_inversion_lengths[key] for key in plot_keys]
        fig = plt.figure(figsize=(12,12))
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

        plot_label = 'Bar plot showing distance from origin cut of reads supporting inversions (up to 100kb). The origin cut is at ' + origin_chr + ':' + str(origin_cut_pos) +' and the origin extends ' + origin_direction + ' from the cut to the primer. Reads supporting inversions extend ' + origin_direction + ' after genomic alignment. Distances greater than ' + str(origin_dist_cutoff) + 'bp are not shown.'
        if len(keys) == 0:
            plot_label = '(No reads supporting inversions found)'

        origin_inversion_hist_plot_obj = PlotObject(
                plot_name = origin_inversion_hist_plot_obj_root,
                plot_title = 'Inversion distance relative to origin cut',
                plot_label = plot_label,
                plot_datas = [('Inversion distance relative to origin cut (all values)',origin_inversion_hist_plot_obj_root + ".txt")]
                )

    origin_deletion_hist_plot_obj = None # plot of deletion lengths at origin
    if len(origin_deletion_lengths) > 0:
        origin_deletion_hist_plot_obj_root = root + ".origin_deletion_histogram"
        keys = sorted(origin_deletion_lengths.keys())
        vals = [origin_deletion_lengths[key] for key in keys]
        with open(origin_deletion_hist_plot_obj_root+".txt","w") as summary:
            summary.write('Deletion_distance\tnumReads\n')
            for key in keys:
                summary.write(str(key) + '\t' + str(origin_deletion_lengths[key]) + '\n')

        plot_keys = [x for x in keys if abs(x) < origin_dist_cutoff]
        plot_vals = [origin_deletion_lengths[key] for key in plot_keys]
        fig = plt.figure(figsize=(12,12))
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

        plot_label = 'Bar plot showing distance from origin cut of reads supporting deletions (up to 100kb). The origin cut is at ' + origin_chr + ':' + str(origin_cut_pos) +' and the origin extends ' + origin_direction + ' from the cut to the primer. Reads supporting deletions extend ' + origin_direction_opposite + ' after genomic alignment. Distances greater than ' + str(origin_dist_cutoff) + 'bp are not shown.'
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

    # plot of coverage of 2500bp after cut
    origin_depth_counts_2500bp_plot_obj_root = root + ".origin_depth_counts_2500bp"
    origin_depth_counts_2500bp = origin_depth_counts_2500bp_total - origin_depth_counts_2500bp
    origin_depth_counts_2500bp_primer = origin_depth_counts_2500bp_total - origin_depth_counts_2500bp_primer

    xs_primer = list(np.arange(-1*(len(origin_depth_counts_2500bp_primer)-1),0,1))
    ys_primer = list(np.flip(origin_depth_counts_2500bp_primer[1:]))
    xs_post = list(np.arange(1,len(origin_depth_counts_2500bp)))
    ys_post = list(origin_depth_counts_2500bp[1:])

    xs = xs_primer + xs_post
    ys = ys_primer + ys_post
    ax = plt.subplot(212)
    if len(origin_depth_counts_2500bp) > 0:
        ax.bar(xs_primer+xs_post,ys_primer+ys_post)
        ax.set_ymargin(0.05)
    else:
        ax.bar(0,0)
    ax.set_ylabel('Read depth (Number of reads)')
    ax.set_title('Read depth of 2500bp after origin')
    line1 = ax.axvline(xs_primer[min_primer_length]+0.5,color='k',ls='dotted')
    line2 = ax.axvline(0,color='r',ls='dotted')
    ax.legend([line1,line2],['min primer length','origin location'],loc='lower right',bbox_to_anchor=(1,0))

    with open(origin_indel_depth_plot_obj_root+"_2500bp.txt","w") as summary:
        summary.write('bpFromOriginCut\treadDepth\n')
        for (xval, yval) in zip(xs,ys):
            summary.write(str(xval)+"\t"+str(yval)+"\n")

    plt.savefig(origin_indel_depth_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(origin_indel_depth_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

    plot_label = 'Read depth surround origin. The top plot shows reads with deletions smaller than 100bp (N='+str(origin_depth_counts_100bp_total)+'). The bottom plot shows reads with deletions smaller than 2500bp (N='+str(origin_depth_counts_2500bp_total)+') . The origin cut is at ' + origin_chr + ':' + str(origin_cut_pos)+'. The vertical red dotted line shows the origin location. The vertical black dotted line shows the minimum requred length (' + str(min_primer_length)+'bp) set by the parameter --min_primer_length.'
    if sum(ys) < 0:
        plot_label = '(No reads found at origin)'

    origin_indel_depth_plot_obj = PlotObject(
            plot_name = origin_indel_depth_plot_obj_root,
            plot_title = 'Read depth around origin',
            plot_label = plot_label,
            plot_datas = [('Depth around origin (100bp)',origin_indel_depth_plot_obj_root + "_100bp.txt"),
                        ('Depth around origin (2500bp)',origin_indel_depth_plot_obj_root + "_2500bp.txt")]
            )

    #make r1/r2/support plots
    r1_r2_support_plot_obj = None # plot of how many reads for which r1 and r2 support each other
    r1_r2_support_dist_plot_obj = None # plot of how far apart the r1/r2 supporting reads were aligned from each other
    r1_r2_no_support_dist_plot_obj = None # plot of how far apart the r1/r2 reads that did NOT support each other (but were on the same chromosome)
    if len(r1_r2_support_status_counts) > 0:
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
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(r1_r2_support_labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(r1_r2_support_labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('R1/R2 support status')
        plt.savefig(r1_r2_support_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(r1_r2_support_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(r1_r2_support_labels,values)])
        r1_r2_support_plot_obj = PlotObject(
                plot_name = r1_r2_support_plot_obj_root,
                plot_title = 'R1/R2 Support Classification',
                plot_label = 'R1/R2 support classification<br>'+plot_count_str,
                plot_datas = [
                    ('R1/R2 support classifications',r1_r2_support_plot_obj_root + ".txt")
                    ]
                )

        #plot of distance between R1 and R2 for R1/R2 supportive read pairs
        r1_r2_support_dist_plot_obj_root = root + ".r1_r2_support_distances"
        keys = sorted(r1_r2_support_distances.keys())
        vals = [r1_r2_support_distances[key] for key in keys]
        with open(r1_r2_support_dist_plot_obj_root+".txt","w") as summary:
            summary.write('R1_R2_distance\tnumReads\n')
            for key in keys:
                summary.write(str(key) + '\t' + str(r1_r2_support_distances[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(keys,vals)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Distance between R1/R2')
        plt.savefig(r1_r2_support_dist_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(r1_r2_support_dist_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_label = 'Bar plot showing distance between R1 and R2 reads for reads with R2 supporting the R1 alignment (max distance for support is ' + str(r1_r2_support_max_distance) + 'bp)'
        if len(keys) == 0:
            plot_label = '(No R1/R2 supported reads)'

        r1_r2_support_dist_plot_obj = PlotObject(
                plot_name = r1_r2_support_dist_plot_obj_root,
                plot_title = 'Distance between R1/R2 Supported Reads',
                plot_label = plot_label,
                plot_datas = [('R1/R2 support distance histogram',r1_r2_support_dist_plot_obj_root + ".txt")]
                )

        #plot of distance between R1 and R2 for R1/R2 NOT supportive read pairs
        r1_r2_no_support_dist_plot_obj_root = root + ".r1_r2_not_supported_distances"
        keys = sorted(r1_r2_not_support_distances.keys())
        vals = [r1_r2_not_support_distances[key] for key in keys]
        with open(r1_r2_no_support_dist_plot_obj_root+".txt","w") as summary:
            summary.write('R1_R2_distance\tnumReads\n')
            for key in keys:
                summary.write(str(key) + '\t' + str(r1_r2_not_support_distances[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(keys,vals)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Distance between R1/R2')
        plt.savefig(r1_r2_no_support_dist_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(r1_r2_no_support_dist_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_label = 'Bar plot showing distance between R1 and R2 reads for reads with R2 NOT supporting the R1 alignment due to incorrect read orientation or exceeding maximum alignment distance (max distance for support is ' + str(r1_r2_support_max_distance) + 'bp)'
        if len(keys) == 0:
            plot_label = '(No R1/R2 not-supported reads)'

        r1_r2_no_support_dist_plot_obj = PlotObject(
                plot_name = r1_r2_no_support_dist_plot_obj_root,
                plot_title = 'Distance between R1/R2 NOT Supported Reads',
                plot_label = plot_label,
                plot_datas = [('R1/R2 not-supported distance histogram',r1_r2_no_support_dist_plot_obj_root + ".txt")]
                )

    with open(info_file,'w') as fout:
        fout.write("\t".join(['total_reads_processed','discarded_reads_duplicates','discarded_reads_not_supported_by_R2','discarded_reads_bad_alignment','final_read_count','discovered_cut_point_count','final_cut_point_count'])+"\n")
        fout.write("\t".join(str(x) for x in [total_reads_processed,dups_removed,bad_alignments_removed,nonsupporting_alignments_removed,final_total_count,found_cut_point_total,final_cut_point_total])+"\n")
        classification_json_str = json.dumps(cut_classification_lookup,separators=(',',':'))
        fout.write(classification_json_str+"\n")

        chr_aln_plot_obj_str = chr_aln_plot_obj.to_json()
        fout.write(chr_aln_plot_obj_str+"\n")

        tlen_plot_obj_str = "None"
        if tlen_plot_obj is not None:
            tlen_plot_obj_str = tlen_plot_obj.to_json()
        fout.write(tlen_plot_obj_str+"\n")

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

        r1_r2_no_support_dist_plot_obj_str = "None"
        if r1_r2_no_support_dist_plot_obj is not None:
            r1_r2_no_support_dist_plot_obj_str = r1_r2_no_support_dist_plot_obj.to_json()
        fout.write(r1_r2_no_support_dist_plot_obj_str+"\n")

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

    return (final_file,final_read_ids_for_cuts,final_cut_counts,cut_classification_lookup,
        chr_aln_plot_obj,tlen_plot_obj,deduplication_plot_obj,tx_order_plot_obj,tx_count_plot_obj,classification_plot_obj,classification_indel_plot_obj,
        origin_indel_hist_plot_obj,origin_inversion_hist_plot_obj,origin_deletion_hist_plot_obj,origin_indel_depth_plot_obj,
        r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,r1_r2_no_support_dist_plot_obj)


def read_command_output(command):
    """
    Runs a shell command and returns an iter to read the output

    param:
        command: shell command to run

    returns:
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
        run_crispresso_on_novel_sites,samtools_command,crispresso_command,n_processes=1):
    """
    Prepares reads for analysis by crispresso2
    Frequently-aligned locations with a min number of reads (crispresso_min_count) are identified
    Reads from these locations are extracted and prepared for analysis by CRISPResso2

    params:
        root: root for written files
        input_fastq_file: input fastq with read sequences
        read_ids_for_crispresso: dict of readID=>cut assignment
        cut_counts_for_crispresso: dict of cut=>count for # reads assigned to each cut
        cut_classification_lookup: dict of cutID=> type (e.g. Linear, Translocation, etc)
        cut_annotations: dict of cut_chr:cut_site->annotation for description of cut [type, origin_status, guide_direction] 
            (type=On-target, Off-target, Known, Casoffinder, etc) (origin_status=Not-origin or Origin:left or Origin:right) (guide_direction='FW' or 'RC' for which direction the guide binds)
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

    returns:
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

        if cut_counts_for_crispresso[this_cut] < crispresso_min_count:
            continue

        if this_cut not in cut_names:
            cut_name = this_cut.replace(" ","_").replace(":","_")
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
        crispresso_cmd = "%s -o %s -n %s --default_min_aln_score %d -a %s -g %s -wc %s -w %s -r1 %s --fastq_output --quantification_window_coordinates %s %s &> %s.log"%(crispresso_command,data_dir,cut_name,crispresso_min_aln_score,
            amp_seq,guide_seq,-15,crispresso_quant_window_size,reads_file,quant_window_coords,processes_str,reads_file)
        
        run_this_one = 'True'
        if not run_crispresso_on_novel_sites and 'Novel' in cut_classification_lookup[cut]:
            run_this_one = 'False: novel location'

        if printed_cuts[cut] < crispresso_min_count:
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

def run_and_aggregate_crispresso(root,crispresso_infos,final_assignment_file,n_processes,skip_failed=True,keep_intermediate=False,can_use_previous_analysis=False):
    """
    Runs CRISPResso2 commands and aggregates output

    params:
        root: root for written files
        crispresso_infos: array of metadata information for CRISPResso
            dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command, run
        final_assignment_file: filename of final assignments for each read id pair
        n_processes: number of processes to run CRISPResso commands on
        skip_failed: if true, failed CRISPResso runs are skipped. Otherwise, if one fails, the program will fail.
        keep_intermediate: whether to keep intermediate files (if False, intermediate files including produced fastqs will be deleted)
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
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
        crispresso_info_fh.write('name\tcut_annotation\tcut_type\tcut_location\treference\treads_printed\tn_total\treads_aligned\treads_unmod\treads_mod\treads_discarded\treads_insertion\treads_deletion\treads_substitution\treads_only_insertion\treads_only_deletion\treads_only_substitution\treads_insertion_and_deletion\treads_insertion_and_substitution\treads_deletion_and_substitution\treads_insertion_and_deletion_and_substitution\tamplicon_sequence\n')
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
    with open(classification_plot_obj_root+".txt",'w') as summary:
        summary.write("\t".join(labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
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
    plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(labels,values)])

    classification_plot_obj = PlotObject(
            plot_name = classification_plot_obj_root,
            plot_title = 'Read classification including CRISPResso analysis',
            plot_label = 'CRISPResso classification<br>' + plot_count_str + '<brReads labeled (not analyzed) had to few reads for CRISPResso analysis. You can adjust the minimum number of reads to run CRISPResso with the parameter "crispresso_min_count".',
            plot_datas = [
                ('CRISPResso classification summary',classification_plot_obj_root + ".txt")
                ]
            )


    with open(crispresso_stats_file,'w') as fout:
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

    params:
        mdz_str: MD:Z string from alignment

    returns:
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
        left_cols = None,
        right_cols = None,
        left_outline_cols=None,
        right_outline_cols=None,
        outline_cols=None,
        outline_labs=None
        ):
    """
    Make a count plot for the most common types of translocations

    params:
        left_labs: (list) labels for the left side translocations
        right_labs: (list) labels for the right side of translocations
        counts: (list) read counts for each translocation
        left_cols: (list) colors for left side
        right_cols: (list) colors for right side
        left_outline_cols: (list) outline colors for left side
        right_outline_cols: (list) outline colors for right side
        outline_cols: (list) colors for outline color legend
        outline_labs: (list) labels for outline color legend

        The length of all params should be the same - one for each translocation to plot
    returns:
        plt: matplotlib plot object
    """

    num_boxes = len(left_labs)

    if left_cols is None:
        left_cols = ['b']*num_boxes

    if right_cols is None:
        right_cols = ['b']*num_boxes

    if left_outline_cols is None:
        left_outline_cols = ['None']*num_boxes

    if right_outline_cols is None:
        right_outline_cols = ['None']*num_boxes

    x_width = 2 - 0.1
    y_height = 1 - 0.1
    ys = range(0,num_boxes)[::-1]
    x1_start = 0
    x2_start = 2
    count_bar_ydiff = 0.3

    fig, (ax1, ax2) = plt.subplots(1,2,sharey=True, figsize=(12,8))
    left_boxes = [Rectangle((x1_start, y), x_width, y_height)
                      for y in ys]
    left_pc = PatchCollection(left_boxes, facecolor=left_cols, edgecolor=left_outline_cols, linewidth=2)

    right_boxes = [Rectangle((x2_start, y), x_width, y_height)
                      for y in ys]

    right_pc = PatchCollection(right_boxes, facecolor=right_cols, edgecolor=right_outline_cols, linewidth=2)

    # Add collection to axes
    ax1.add_collection(left_pc)
    ax1.add_collection(right_pc)
    ax1.set_ylim(0,num_boxes)

    max_right_len = max([len(lab) for lab in right_labs])

    ax1.set_xlim(-2,10)

    for ind in range(num_boxes):
        ax1.text(-0.1,ys[ind]+y_height/2,left_labs[ind],ha='right',va='center')
        ax1.text(4,ys[ind]+y_height/2,right_labs[ind],ha='left',va='center')

    ax1.axis('off')

    if outline_labs is not None:
        legend_patches = []
        for col,lab in zip(outline_cols,outline_labs):
            legend_patches.append(Patch(facecolor='None',edgecolor=col,label=lab))
        ax1.legend(handles=legend_patches,loc="lower center", bbox_to_anchor=(0.5, -0.2))

    rects = []
    for ind in range(num_boxes):
        rects.append(Rectangle((0,ys[ind]+count_bar_ydiff),counts[ind],y_height-(count_bar_ydiff*2)))
        ax2.text(x=counts[ind],y=ys[ind]+y_height/2,s=" " + str(counts[ind]),ha='left',va='center')

    pc = PatchCollection(rects)
    ax2.add_collection(pc)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_frame_on(False)

    ax2.set_xscale("log")
    ax2.set_xlim(1,max(5,max(counts)))
    ax2.set_xlabel('Number of reads')

    return plt


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
        plot_str = "<div class='card text-center mb-2' id='plot"+str(idx)+"'>\n\t<div class='card-header'>\n"
        plot_str += "<h5>"+plot_obj.title+"</h5>\n"
        plot_str += "</div>\n"
        plot_str += "<div class='card-body'>\n"
        plot_str += "<a href='"+data_path+plot_obj.name+".pdf'><img src='"+data_path + plot_obj.name + ".png' width='80%' ></a>\n"
        plot_str += "<label>"+plot_obj.label+"</label>\n"
        for (plot_data_label,plot_data_path) in plot_obj.datas:
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
        args:
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

    params:
        root: root for written files
        fastq_r1: R1 reads to trim
        fastq_r1_trimmed: output file to write R1 to
        min_primer_aln_score: minimum score for alignment between primer/origin sequence and read sequence
        min_primer_length: minimum length of sequence that matches between the primer/origin sequence and the read sequence
        ssw_align_primer: SSWPrimerAlignment object for ssw alignment of primer

    returns:
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

        primer_seq, trimmed_seq = trimLeftPrimerFromRead(f1_seq_line,ssw_align_primer,min_primer_aln_score)
        if len(primer_seq) > min_primer_length:
            post_trim_read_count += 1
            new_f1_qual_line = f1_qual_line[len(primer_seq):]
            new_f1_id_line = f1_id_line + ' primer=' + primer_seq

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

    params:
        root: root for written files
        fastq_r1: R1 reads to trim
        fastq_r2: R2 reads to trim
        fastq_r1_trimmed: output file to write R1 to
        fastq_r2_trimmed: output file to write R2 to
        min_primer_aln_score: minimum score for alignment between primer/origin sequence and read sequence
        min_primer_length: minimum length of sequence that matches between the primer/origin sequence and the read sequence
        ssw_align_primer: SSWPrimerAlignment object for ssw alignment of primer
        ssw_align_primer_rc: SSWPrimerAlignment object for ssw alignment of reverse complement primer (to r2)

    returns:
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

        primer_seq, trimmed_seq = trimLeftPrimerFromRead(f1_seq_line,ssw_align_primer,min_primer_aln_score)
        if len(primer_seq) > min_primer_length:
            post_trim_read_count += 1
            new_f1_qual_line = f1_qual_line[len(primer_seq):]
            new_f1_id_line = f1_id_line + ' primer=' + primer_seq

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

def trimLeftPrimerFromRead(read, ssw_align, min_score=40):
    """
    Trims a primer from a given read
    E.g. for --primer--read-- returns the portion that aligns to the primer (and before) as the trimmed_primer_seq and to the right as the trimmed_read_seq

    params:
        primer: string
        read: string
        ssw_align: ssw alignment object
        min_score: minimum alignment score between primer and sequence

    returns:
        trimmed_primer_seq: portion of read that was aligned to the primer
        trimmed_read_seq: portion of the read after the primer

    """

    ssw_res = ssw_align.align(read)
    if ssw_res.contents.nScore < min_score:
        ssw_align.destroy(ssw_res)
        return '',read
    else:
        trimmed_primer_seq = read[:ssw_res.contents.nQryEnd+1]
        trimmed_read_seq = read[ssw_res.contents.nQryEnd+1:]
        ssw_align.destroy(ssw_res)
        return trimmed_primer_seq,trimmed_read_seq

def trimRightPrimerFromRead(read, ssw_align, min_score=40):
    """
    Trims a primer from a given read
    E.g. for --read--primer-- returns the portion that aligns to the primer (and after) as the trimmed_primer_seq and to the left as the trimmed_read_seq

    params:
        read: string
        ssw_align: ssw alignment object
        mismatch_score: mismatch/gap score for alignment

    returns:
        trimmed_read_seq: portion of the read after the primer
        trimmed_primer_seq: portion of read that was aligned to the primer

    """

    ssw_res = ssw_align.align(read)
#    print('================')
#    print(f'{min_score=}')
#    print(f'{read=}')
#    print(f'{ssw_align.primer_seq=}')
#    print(f'{ssw_res.contents.nScore=}')
#    print(f'{ssw_res.contents.nRefBeg=}')
#    print(f'{ssw_res.contents.nRefEnd=}')
#    print(f'{ssw_res.contents.nQryBeg=}')
#    print(f'{ssw_res.contents.nQryEnd=}')

    if ssw_res.contents.nScore < min_score:
        ssw_align.destroy(ssw_res)
        return read,''
    else:
        trimmed_read_seq = read[:ssw_res.contents.nQryBeg]
        trimmed_primer_seq = read[ssw_res.contents.nQryBeg:]
        ssw_align.destroy(ssw_res)
        return trimmed_read_seq,trimmed_primer_seq

def to_int(seq, list_letters, dict_letter_lookup):
    """
    Translate a letter sequence into numbers for ssw alignment
    params:
        seq: string, letters to translate
        list_letters: list of all letters
        dict_letter_lookup: dict of letter > number for lookup

    returns:
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

    params:
        guide_seq_aln (str): The guide aligned sequence e.g.     --ATTA-
        ref_seq_aln (str): The reference aligned sequence e.g.   AAATTGGGG
        cut_ind (int): The zero-based index in the the guide seq string of the base to the left where a cut would occur (e.g. for Cas9 and reverse-complement guide PAMGGGGGGGGGGG the cut_ind would be given as 5)

    returns:
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
    settings, logger = parse_settings(sys.argv)
    try:
        main(settings,logger)
    except Exception as e:
        if logger.isEnabledFor(logging.DEBUG):
            logger.error(e, exc_info=True)
        else:
            logger.error(e)
        raise e
