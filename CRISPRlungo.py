import argparse
from collections import defaultdict
import gzip
import io
import json
import logging
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import re
import subprocess
import sys
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoShared

__version__ = "v0.0.3"

def main():
    settings = parse_settings(sys.argv)

    #data structures for plots for report
    summary_plot_objects=[]  # list of PlotObjects for plotting

    assert_dependencies(
            samtools_command=settings['samtools_command'],
            bowtie2_command=settings['bowtie2_command'],
            crispresso_command=settings['crispresso_command'],
            casoffinder_command=settings['casoffinder_command']
            )

    av_read_length = get_av_read_len(settings['fastq_r1'])
    num_reads_input = get_num_reads_fastq(settings['fastq_r1'])
    logging.info('%d reads in input'%num_reads_input)

    cut_sites = settings['cuts']
    cut_annotations = {}
    if len(settings['on_target_cut_sites']) == 0:
        for cut_site in settings['cuts']:
            cut_annotations[cut_site] = 'Known'
    else:
        for cut_site in settings['cuts']:
            cut_annotations[cut_site] = 'Off-target'
        cut_sites.extend(settings['on_target_cut_sites'])
        for cut_site in settings['on_target_cut_sites']:
            cut_annotations[cut_site] = 'On-target'


    if settings['PAM'] is not None and settings['on-target_guides'] is not None:
        casoffinder_cut_sites = get_cut_sites_casoffinder(
                root = settings['root'],
                genome=settings['genome'],
                pam=settings['PAM'],
                guides=settings['on-target_guides'],
                cleavage_offset=settings['cleavage_offset'],
                num_mismatches=settings['casoffinder_num_mismatches'],
                casoffinder_command=settings['casoffinder_command']
                )
        cut_sites.extend(casoffinder_cut_sites)
        for cut_site in casoffinder_cut_sites:
            cut_annotations[cut_site] = 'Casoffinder'

    primer_chr = ""
    primer_loc = -1
    if settings['primer_site']:
        primer_info = settings['primer_site'].split("_")
        primer_chr = primer_info[0]
        primer_loc = int(primer_info[1])
    primer_seq = None
    if settings['primer_seq']:
        #attempt to align primer to genome
        primer_aln_cmd = '%s --seed 2248 --end-to-end -x %s -c %s' %(settings['bowtie2_command'],settings['bowtie2_genome'],settings['primer_seq'])
        primer_aln_results = subprocess.check_output(primer_aln_cmd,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding).split("\n")
        if len(primer_aln_results) > 5 and primer_aln_results[3].strip() == '1 (100.00%) aligned exactly 1 time' and len(primer_aln_results[-2].split("\t")) > 5:
            sam_line_els = primer_aln_results[-2].split("\t")
            primer_chr = sam_line_els[2]
            primer_loc = int(sam_line_els[3])
            logging.info('Setting genomic coordinates for primer aligned to %s:%s'%(primer_chr,primer_loc))
        primer_seq = settings['primer_seq']

    target_padding = settings['alignment_extension']
    target_length = av_read_length
    if target_padding < 0:
        target_length += target_padding
        target_padding = 0
        target_padding += av_read_length

    if len(cut_sites) > 0:
        (custom_index_fasta,target_names,target_info) = make_artificial_targets(
                root = settings['root']+'.customTargets',
                cuts=cut_sites,
                cut_annotations=cut_annotations,
                genome=settings['genome'],
                target_length=target_length,
                target_padding=target_padding,
                primer_chr=primer_chr,
                primer_loc=primer_loc,
                primer_seq=primer_seq,
                add_non_primer_cut_targets=settings['add_non_primer_cut_targets'],
                samtools_command=settings['samtools_command'],
                bowtie2_command=settings['bowtie2_command'],
                bowtie2_threads = settings['n_processes']
                )
    else:
        custom_index_fasta = None
        target_names = []
        target_info = {}


    curr_r1_file = settings['fastq_r1'] #if alignment to genome happens first, the input for artificial target mapping will be reads that don't align to the genome
    curr_r1_description = "input R1 reads"
    curr_r2_file = settings['fastq_r2']
    curr_r2_description = "input R2 reads"

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
        curr_r1_description = "UMI-added R1 reads"
        curr_r2_file = umi_r2
        curr_r1_description = "UMI-added R2 reads"

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
        curr_r1_description = "UMI-deduplicated R1 reads"
        curr_r2_file = dedup_r2
        curr_r2_description = "UMI-deduplicated R2 reads"

    filter_primer_match_count = 0
    if settings['primer_min_bp_match_pct'] > 0:
        #filter out reads that don't contain primer sequence
        (primer_filter_r1, primer_filter_r2, filter_tot_read_count, filter_count_with_primer_match, filter_count_with_primer_post_match, post_primer_filter_count, filter_on_primer_plot_obj
                ) = filter_on_primer_seq(
                        root = settings['root']+'.primer_filter',
                        fastq_r1 = curr_r1_file,
                        fastq_r2 = curr_r2_file,
                        primer_seq = settings['primer_seq'],
                        primer_min_bp_match_pct = settings['primer_min_bp_match_pct'],
                        primer_post_seq = settings['primer_post_seq'],
                        primer_post_min_bp_match_pct = settings['primer_post_min_bp_match_pct'],
                        can_use_previous_analysis = settings['can_use_previous_analysis']
                    )
        curr_r1_file = primer_filter_r1
        curr_r1_description = "primer-filtered R1 reads"
        curr_r2_file = primer_filter_r2
        curr_r2_description = "primer-filtered R2 reads"
        if filter_on_primer_plot_obj is not None:
            filter_on_primer_plot_obj.order = 1
            summary_plot_objects.append(filter_on_primer_plot_obj)

    genome_aligned_count = 0 # number of reads aligned to genome
    custom_aligned_count = 0 # number of reads aligned to custom targets
    frag_linear_count = 0 # number of reads fragmented and classified as linear
    frag_translocation_count = 0 # number of reads fragmented and classified as translocations
    frag_large_deletion_count = 0 # number of reads fragmented and classified as large deletions
    frag_large_inversion_count = 0 # number of reads fragmented and classified as large inversions
    frag_unidentified_count = 0 # number of reads fragmented and classified as unidentified

    #keep track of where every read is assigned
    #these lists keep track of the files containing assignments from the separate analyses
    r1_assignment_files = []
    r2_assignment_files = []

    if settings['perform_unbiased_analysis']:
        (genome_r1_assignments,genome_r2_assignments,genome_unmapped_r1, genome_unmapped_r2, genome_aligned_count, genome_mapped_bam_file,genome_chr_aln_plot_obj,genome_tlen_plot_object
            ) = align_reads(
                    root = settings['root']+'.genomeAlignment',
                    fastq_r1 = curr_r1_file,
                    fastq_r2 = curr_r2_file,
                    reads_name = curr_r1_description.replace(" R1",""),
                    bowtie2_reference = settings['bowtie2_genome'],
                    reference_name = 'Genome',
                    target_info = target_info,
                    bowtie2_command = settings['bowtie2_command'],
                    bowtie2_threads = settings['n_processes'],
                    samtools_command = settings['samtools_command'],
                    keep_intermediate = settings['keep_intermediate'],
                    can_use_previous_analysis = settings['can_use_previous_analysis']
                    )
        r1_assignment_files.append(genome_r1_assignments)
        if genome_r2_assignments is not None:
            r2_assignment_files.append(genome_r2_assignments)

        if genome_tlen_plot_object is not None:
            genome_tlen_plot_object.order = 2
            summary_plot_objects.append(genome_tlen_plot_object)
        if genome_chr_aln_plot_obj is not None:
            genome_chr_aln_plot_obj.order = 5
            summary_plot_objects.append(genome_chr_aln_plot_obj)

        curr_r1_file = genome_unmapped_r1
        curr_r1_description = "R1 reads not aligned to genome"
        curr_r2_file = genome_unmapped_r2
        curr_r2_description = "R2 reads not aligned to genome"

    if len(target_names) > 0:
        #first, align the unaligned r1s
        (custom_r1_assignments,none_file,custom_unmapped_r1, none_file2, custom_r1_aligned_count, custom_r1_mapped_bam_file,custom_r1_chr_aln_plot_obj,none_custom_tlen_plot_object
            ) = align_reads(
                    root = settings['root']+'.r1.customTargetAlignment',
                    fastq_r1 = curr_r1_file,
                    fastq_r2 = None,
                    reads_name = curr_r1_description,
                    bowtie2_reference=custom_index_fasta,
                    reference_name='Custom targets',
                    target_info = target_info,
                    arm_min_seen_bases = settings['arm_min_seen_bases'],
                    arm_min_matched_start_bases = settings['arm_min_matched_start_bases'],
                    bowtie2_command=settings['bowtie2_command'],
                    bowtie2_threads=settings['n_processes'],
                    samtools_command=settings['samtools_command'],
                    keep_intermediate = settings['keep_intermediate'],
                    can_use_previous_analysis = settings['can_use_previous_analysis']
                    )
        r1_assignment_files.append(custom_r1_assignments)
        custom_aligned_count += custom_r1_aligned_count
        curr_r1_file = custom_unmapped_r1
        curr_r1_description = "R1 reads not aligned to custom targets"

        if custom_r1_chr_aln_plot_obj is not None:
            custom_r1_chr_aln_plot_obj.order = 10
            summary_plot_objects.append(custom_r1_chr_aln_plot_obj)

        #if input is paired, align second reads to the custom targets as welll
        if curr_r2_file is not None:
            (custom_r2_assignments,none_file,custom_unmapped_r2, none_file2, custom_r2_aligned_count, custom_r2_mapped_bam_file,custom_r2_chr_aln_plot_obj,none_custom_tlen_plot_object
                ) = align_reads(
                        root = settings['root']+'.r2.customTargetAlignment',
                        fastq_r1 = curr_r2_file,
                        fastq_r2 = None,
                        reads_name = curr_r2_description,
                        bowtie2_reference=custom_index_fasta,
                        reference_name='Custom targets',
                        target_info = target_info,
                        arm_min_seen_bases = settings['arm_min_seen_bases'],
                        arm_min_matched_start_bases = settings['arm_min_matched_start_bases'],
                        bowtie2_command=settings['bowtie2_command'],
                        bowtie2_threads=settings['n_processes'],
                        samtools_command=settings['samtools_command'],
                        keep_intermediate = settings['keep_intermediate'],
                        can_use_previous_analysis = settings['can_use_previous_analysis']
                        )
            r2_assignment_files.append(custom_r2_assignments)
            custom_aligned_count += custom_r2_aligned_count

            if custom_r2_chr_aln_plot_obj is not None:
                custom_r2_chr_aln_plot_obj.order = 11
                summary_plot_objects.append(custom_r2_chr_aln_plot_obj)

            curr_r2_file = custom_unmapped_r2

    if settings['perform_unbiased_analysis']:
        #chop reads
        (frag_r1_assignments,frag_r2_assignments,frag_linear_count,frag_translocation_count,frag_large_deletion_count,frag_large_inversion_count,frag_unidentified_count,frags_plot_obj
            ) = chop_reads(
                    root = settings['root']+".frags",
                    unmapped_reads_fastq_r1 = curr_r1_file,
                    unmapped_reads_fastq_r2 = curr_r2_file,
                    bowtie2_genome = settings['bowtie2_genome'],
                    fragment_size = settings['fragment_size'],
                    fragment_step_size = settings['fragment_step_size'],
                    min_seen_frag_cutoff = settings['min_seen_frag_cutoff'],
                    fragment_min_mapq = settings['fragment_min_mapq'],
                    bowtie2_command = settings['bowtie2_command'],
                    bowtie2_threads = settings['n_processes'],
                    samtools_command = settings['samtools_command'],
                    can_use_previous_analysis = settings['can_use_previous_analysis']
                    )
        r1_assignment_files.append(frag_r1_assignments)
        if frag_r2_assignments is not None:
            r2_assignment_files.append(frag_r2_assignments)

        if frags_plot_obj is not None:
            frags_plot_obj.order = 16
            summary_plot_objects.append(frags_plot_obj)
    else:
        (unmapped_r1_assignments, unmapped_r2_assignments, unmapped_r1_count, unmapped_r2_count
                ) = prep_unmapped_assignments(
                    root = settings['root']+".unmapped",
                    unmapped_reads_fastq_r1 = curr_r1_file,
                    unmapped_reads_fastq_r2 = curr_r2_file,
                    can_use_previous_analysis = settings['can_use_previous_analysis']
                    )
        r1_assignment_files.append(unmapped_r1_assignments)
        if unmapped_r2_assignments is not None:
            r2_assignment_files.append(unmapped_r2_assignments)


    (final_assignment_file,r1_read_ids_for_crispresso,r1_cut_counts_for_crispresso,r2_read_ids_for_crispresso,r2_cut_counts_for_crispresso,deduplication_plot_obj,r1_source_plot_obj,r1_classification_plot_obj,r2_source_plot_obj,r2_classification_plot_obj,r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,r1_r2_no_support_dist_plot_obj
        ) = make_final_read_assignments(
                root = settings['root']+'.final',
                r1_assignment_files = r1_assignment_files,
                r2_assignment_files = r2_assignment_files,
                cut_sites = cut_sites,
                cut_annotations = cut_annotations,
                target_info = target_info,
                r1_r2_support_max_distance = settings['r1_r2_support_max_distance'],
                cut_merge_dist = 20,
                genome_map_resolution = 1000000,
                dedup_based_on_aln_pos = settings['dedup_input_based_on_aln_pos_and_UMI'],
                discard_reads_without_r2_support = settings['discard_reads_without_r2_support'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )

    deduplication_plot_obj.order=20
    summary_plot_objects.append(deduplication_plot_obj)
    r1_source_plot_obj.order=25
    summary_plot_objects.append(r1_source_plot_obj)
    r1_classification_plot_obj.order=26
    summary_plot_objects.append(r1_classification_plot_obj)
    if r2_source_plot_obj is not None:
        r2_source_plot_obj.order=27
        summary_plot_objects.append(r2_source_plot_obj)
    if r2_classification_plot_obj is not None:
        r2_classification_plot_obj.order=28
        summary_plot_objects.append(r2_classification_plot_obj)
    if r1_r2_support_plot_obj is not None:
        r1_r2_support_plot_obj.order=29
        summary_plot_objects.append(r1_r2_support_plot_obj)
    if r1_r2_support_dist_plot_obj is not None:
        r1_r2_support_dist_plot_obj.order=30
        summary_plot_objects.append(r1_r2_support_dist_plot_obj)
    if r1_r2_no_support_dist_plot_obj is not None:
        r1_r2_no_support_dist_plot_obj.order=31
        summary_plot_objects.append(r1_r2_no_support_dist_plot_obj)

    crispresso_infos = [] #meta info about the crispresso runs
            #dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command
    logging.info('Preparing R1 reads for analysis')
    r1_crispresso_infos = prep_crispresso2(
                root = settings['root']+'.CRISPResso_r1',
                input_fastq_file = settings['fastq_r1'],
                read_ids_for_crispresso = r1_read_ids_for_crispresso,
                cut_counts_for_crispresso = r1_cut_counts_for_crispresso,
                av_read_length = av_read_length,
                genome = settings['genome'],
                genome_len_file = settings['genome']+'.fai',
                crispresso_min_count = settings['crispresso_min_count'],
                crispresso_min_aln_score = settings['crispresso_min_aln_score'],
                samtools_command=settings['samtools_command'],
                crispresso_command=settings['crispresso_command'],
                )
    for info in r1_crispresso_infos:
        info['name'] = 'r1_'+info['name']
        crispresso_infos.append(info)

    if settings['fastq_r2'] is not None:
        logging.info('Preparing R2 reads for analysis')
        r2_crispresso_infos = prep_crispresso2(
                    root = settings['root']+'.CRISPResso_r2',
                    input_fastq_file = settings['fastq_r2'],
                    read_ids_for_crispresso = r2_read_ids_for_crispresso,
                    cut_counts_for_crispresso = r2_cut_counts_for_crispresso,
                    av_read_length = av_read_length,
                    genome = settings['genome'],
                    genome_len_file = settings['genome']+'.fai',
                    crispresso_min_count = settings['crispresso_min_count'],
                    crispresso_min_aln_score = settings['crispresso_min_aln_score'],
                    samtools_command=settings['samtools_command'],
                    crispresso_command=settings['crispresso_command'],
                    )

        for info in r2_crispresso_infos:
            info['name'] = 'r2_'+info['name']
            crispresso_infos.append(info)

    crispresso_results = run_and_aggregate_crispresso(
                root = settings['root']+".CRISPResso",
                crispresso_infos = crispresso_infos,
                n_processes = settings['n_processes'],
                min_count_to_run_crispresso =settings['crispresso_min_count'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )

    labels = ["Input Reads","Aligned Custom Targets","Aligned Genome","Fragmented Linear","Fragmented Translocations","Fragmented Large Deletions","Fragmented Large Inversions","Fragmented Unidentified"]
    values = [num_reads_input,custom_aligned_count,genome_aligned_count,frag_linear_count,frag_translocation_count,frag_large_deletion_count,frag_large_inversion_count,frag_unidentified_count]
    alignment_summary_root = settings['root']+".alignmentSummary"
    with open(alignment_summary_root+".txt",'w') as summary:
        summary.write("\t".join(labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    pie_values = []
    pie_labels = []
    for i in range(1,len(labels)):
        if values[i] > 0:
            pie_values.append(values[i])
            pie_labels.append(labels[i]+"\n("+str(values[i])+")")
    ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
    plot_name = alignment_summary_root
    plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')
    plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(labels,values)])


    summary_plot_objects.append(
            PlotObject(plot_name = plot_name,
                    plot_title = 'Alignment Summary',
                    plot_label = 'Pie chart showing final assignment of reads.<br>' + plot_count_str,
                    plot_datas = [('Alignment summary',alignment_summary_root + ".txt")]
                    ))

    make_report(report_file=settings['root']+".html",
            report_name = 'Report',
            crisprlungo_folder = '',
            crispresso_run_names = crispresso_results['run_names'],
            crispresso_sub_html_files = crispresso_results['run_sub_htmls'],
            summary_plot_objects = summary_plot_objects,
            )


    logging.info('Successfully completed!')

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



    parser.add_argument('--alignment_extension', type=int, help='Number of bp to extend beyond av read length around cut site for custom index', default=50)

    parser.add_argument('--cuts','--cut_sites', nargs='*', help='Cut sites in the form chr1:234 (multiple cuts are separated by spaces)', default=[])
    parser.add_argument('--on_target_cut_sites', nargs='*', help='On-target cut sites in the form chr1:234 (multiple cuts are separated by spaces)', default=[])

    #for finding offtargets with casoffinder
    parser.add_argument('--PAM', type=str, help='PAM for in-silico off-target search', default=None)
    parser.add_argument('--on_target_guides', nargs='*', help='Spacer sequences of guides for in-silico off-target search (multiple guide sequences are separated by spaces)', default=[])
    parser.add_argument('--casoffinder_num_mismatches', type=int, help='Number of casoffinder mismatches for in-silico off-target search', default=5)
    parser.add_argument('--cleavage_offset', type=int, help='Position where cleavage occurs, for in-silico off-target search (relative to end of spacer seq -- for Cas9 this is -3)', default=-3)

    #specify primer information
    parser.add_argument('--primer_site', type=str, help='Site of genomic primer, of the form chr1:234',default=None)
    parser.add_argument('--primer_seq', type=str, help='Sequence of primer. Either --primer_site or --primer_seq should be specified.',default=None)
    parser.add_argument('--add_non_primer_cut_targets', help='Whether to add targets for cuts without a primer. If false, only primer-cut1 pairings will be added. If true, cut1-cut2 pairings will be added.',action='store_true')
    parser.add_argument('--primer_min_bp_match_pct', type=int, help='Min percent of bases that must match primer_seq. If lower than this percent of bases match, the read will be filtered from output files.', default=0)
    parser.add_argument('--primer_post_seq', type=str, help='Sequence immedately following primer sequence to check at beginning of R1', default=None)
    parser.add_argument('--primer_post_min_bp_match_pct', type=int, help='Min percent of bases that must match primer_post_seq. If lower than this percent of bases match, the read will be filtered from output files.', default=0)

    #min alignment cutoffs for alignment to each arm/side of read
    parser.add_argument('--arm_min_seen_bases', type=int, help='Number of bases that are required to be seen on each "side" of translocated reads. E.g. if a artificial target represents a translocation between chr1 and chr2, arm_min_seen_bases would have to be seen on chr1 as well as on chr2 for the read to be counted.', default=15)
    parser.add_argument('--arm_min_matched_start_bases', type=int, help='Number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each "side" of artifical targets. E.g. if a artificial target represents a translocation between chr1 and chr2, the first arm_min_matched_start_bases of the read would have to match exactly to chr1 and the last arm_min_matched_start_bases of the read would have to match exactly to chr2', default=10)

    #CRISPResso settings
    parser.add_argument('--run_crispresso_genome_sites', help='Boolean for whether to run CRISPResso on highly-aligned genomic locations (if flag is not set (or false), CRISPResso will only be run at highly-aligned cut sites and artificial targets', action='store_true')
    parser.add_argument('--crispresso_min_count', type=int, help='Min number of reads required to be seen at a site for it to be analyzed by CRISPResso', default=50)
    parser.add_argument('--crispresso_min_aln_score', type=int, help='Min alignment score to reference sequence for quantification by CRISPResso', default=20)

    #sub-command parameters
    parser.add_argument('--samtools_command', help='Command to run samtools', default='samtools')
    parser.add_argument('--bowtie2_command', help='Command to run bowtie2', default='bowtie2')
    parser.add_argument('--crispresso_command', help='Command to run crispresso', default='CRISPResso')
    parser.add_argument('--casoffinder_command', help='Command to run casoffinder', default='cas-offinder')
    parser.add_argument('--n_processes', type=str, help='Number of processes to run on (may be set to "max")', default='1')

    #fragmentation parameters
    parser.add_argument('--perform_unbiased_analysis', help='If set, an unbiased search for translocations is performed by fragmenting and analyzing reads that do not align to the genome or to custom targets.', action='store_true')
    #these two parameters control the size and stride of fragment creation
    parser.add_argument('--fragment_size', type=int, help='Fragment size (bp) for fragment creation', default=20)
    parser.add_argument('--fragment_step_size', type=int, help='Step size (bp) between fragments for fragment creation', default=10)
    parser.add_argument('--min_seen_frag_cutoff', type=int, help='Minimum number of fragments aligned to a location for the location to be considered as a possible alignment location for a read (e.g. if a read only has 3 fragments aligning to location A and 1 fragment aligning to location B, a min_seen_frag_cutoff of 2 would not allow position B', default=2)
    parser.add_argument('--fragment_min_mapq', type=int, help='Minimum mapq (mapping quality) for a fragment to be considered. Increasing this value may reduce sensitivity for fragments aligned to repetitive regions. Any value greater than 1 will likely eliminate multi-mapped reads.', default=30)

    #umi settings
    parser.add_argument('--dedup_input_on_UMI', help='If set, input reads will be deduplicated based on UMI before alignment', action='store_true')
    parser.add_argument('--dedup_input_based_on_aln_pos_and_UMI', help='If set, input reads will be deduplicated based on alignment position and UMI', action='store_true')
    parser.add_argument('--umi_regex', type=str, help='String specifying regex that UMI must match', default='NNWNNWNNN')

    #R1/R2 support settings
    parser.add_argument('--r1_r2_support_max_distance', type=int, help='Max distance between r1 and r2 for the read pair to be classified as "supported" by r2', default=10000)
    parser.add_argument('--discard_reads_without_r2_support', help='If set, reads without r2 support will be discarded from final analysis and counts', action='store_true')

    parser.add_argument('--genome', help='Genome sequence file for alignment. This should point to a file ending in ".fa", and the accompanying index file (".fai") should exist.', default=None)
    parser.add_argument('--bowtie2_genome', help='Bowtie2-indexed genome file.',default=None)

    parser.add_argument('--fastq_r1', help='Input fastq r1 file', default=None)
    parser.add_argument('--fastq_r2', help='Input fastq r2 file', default=None)
    parser.add_argument('--fastq_umi', help='Input fastq umi file', default=None)

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

    #logger = logging.getLogger(__name__)
    logging_level = logging.INFO
    if settings['debug']:
        logging_level=logging.DEBUG

    log_formatter = logging.Formatter("%(asctime)s:%(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.basicConfig(
            level = logging_level,
            format = "%(asctime)s:%(levelname)s: %(message)s",
            filename = settings['root']+".log",
            filemode = 'w'
            )
    ch = logging.StreamHandler()
    ch.setFormatter(log_formatter)
    logging.getLogger().addHandler(ch)

    logging.info('CRISPRlungo ' + __version__)

    logging.info('Parsing settings file')

    settings['keep_intermediate'] = cmd_args.keep_intermediate
    if 'keep_intermediate' in settings_file_args:
        settings['keep_intermediate'] = (settings_file_args['keep_intermediate'].lower() == 'true')
        settings_file_args.pop('keep_intermediate')


    settings['alignment_extension'] = cmd_args.alignment_extension
    if 'alignment_extension' in settings_file_args:
        settings['alignment_extension'] = int(settings_file_args['alignment_extension'])
        settings_file_args.pop('alignment_extension')

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

    settings['PAM'] = cmd_args.PAM
    if 'PAM' in settings_file_args:
        settings['PAM'] = settings_file_args['PAM']
        settings_file_args.pop('PAM')

    settings['on_target_guides'] = cmd_args.on_target_guides
    if 'on_target_guides' in settings_file_args:
        settings['on_target_guides'] = settings_file_args['on_target_guides'].split(" ")
        settings_file_args.pop('on_target_guides')

    settings['casoffinder_num_mismatches'] = cmd_args.casoffinder_num_mismatches
    if 'casoffinder_num_mismatches' in settings_file_args:
        settings['casoffinder_num_mismatches'] = int(settings_file_args['casoffinder_num_mismatches'])
        settings_file_args.pop('casoffinder_num_mismatches')

    settings['cleavage_offset'] = cmd_args.cleavage_offset
    if 'cleavage_offset' in settings_file_args:
        settings['cleavage_offset'] = int(settings_file_args['cleavage_offset'])
        settings_file_args.pop('cleavage_offset')

    settings['primer_site'] = cmd_args.primer_site
    if 'primer_site' in settings_file_args:
        settings['primer_site'] = settings_file_args['primer_site']
        settings_file_args.pop('primer_site')

    settings['primer_seq'] = cmd_args.primer_seq
    if 'primer_seq' in settings_file_args:
        settings['primer_seq'] = settings_file_args['primer_seq']
        settings_file_args.pop('primer_seq')

    settings['add_non_primer_cut_targets'] = cmd_args.add_non_primer_cut_targets
    if 'add_non_primer_cut_targets' in settings_file_args:
        settings['add_non_primer_cut_targets'] = (settings_file_args['add_non_primer_cut_targets'].lower() == 'true')
        settings_file_args.pop('add_non_primer_cut_targets')

    # if primer is not specified, add non-primer cut targets
    if not settings['primer_seq'] and not settings['primer_site']:
        settings['add_non_primer_cut_targets'] = True

    settings['primer_min_bp_match_pct'] = cmd_args.primer_min_bp_match_pct
    if 'primer_min_bp_match_pct' in settings_file_args:
        settings['primer_min_bp_match_pct'] = int(settings_file_args['primer_min_bp_match_pct'])
        settings_file_args.pop('primer_min_bp_match_pct')

    settings['primer_post_seq'] = cmd_args.primer_post_seq
    if 'primer_post_seq' in settings_file_args:
        settings['primer_post_seq'] = settings_file_args['primer_post_seq']
        settings_file_args.pop('primer_post_seq')

    settings['primer_post_min_bp_match_pct'] = cmd_args.primer_post_min_bp_match_pct
    if 'primer_post_min_bp_match_pct' in settings_file_args:
        settings['primer_post_min_bp_match_pct'] = int(settings_file_args['primer_post_min_bp_match_pct'])
        settings_file_args.pop('primer_post_min_bp_match_pct')

    settings['arm_min_seen_bases'] = cmd_args.arm_min_seen_bases
    if 'arm_min_seen_bases' in settings_file_args:
        settings['arm_min_seen_bases'] = int(settings_file_args['arm_min_seen_bases'])
        settings_file_args.pop('arm_min_seen_bases')

    settings['arm_min_matched_start_bases'] = cmd_args.arm_min_matched_start_bases
    if 'arm_min_matched_start_bases' in settings_file_args:
        settings['arm_min_matched_start_bases'] = int(settings_file_args['arm_min_matched_start_bases'])
        settings_file_args.pop('arm_min_matched_start_bases')

    settings['run_crispresso_genome_sites'] = cmd_args.run_crispresso_genome_sites
    if 'run_crispresso_genome_sites' in settings_file_args:
        settings['run_crispresso_genome_sites'] = (settings_file_args['run_crispresso_genome_sites'].lower() == 'true')
        settings_file_args.pop('run_crispresso_genome_sites')

    settings['crispresso_min_count'] = cmd_args.crispresso_min_count
    if 'crispresso_min_count' in settings_file_args:
        settings['crispresso_min_count'] = int(settings_file_args['crispresso_min_count'])
        settings_file_args.pop('crispresso_min_count')

    settings['crispresso_min_aln_score'] = cmd_args.crispresso_min_aln_score
    if 'crispresso_min_aln_score' in settings_file_args:
        settings['crispresso_min_aln_score'] = int(settings_file_args['crispresso_min_aln_score'])
        settings_file_args.pop('crispresso_min_aln_score')


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

    settings['perform_unbiased_analysis'] = cmd_args.perform_unbiased_analysis
    if 'perform_unbiased_analysis' in settings_file_args:
        settings['perform_unbiased_analysis'] = (settings_file_args['perform_unbiased_analysis'].lower() == 'true')
        settings_file_args.pop('perform_unbiased_analysis')

    settings['fragment_size'] = cmd_args.fragment_size
    if 'fragment_size' in settings_file_args:
        settings['fragment_size'] = int(settings_file_args['fragment_size'])
        settings_file_args.pop('fragment_size')

    settings['fragment_step_size'] = cmd_args.fragment_step_size
    if 'fragment_step_size' in settings_file_args:
        settings['fragment_step_size'] = int(settings_file_args['fragment_step_size'])
        settings_file_args.pop('fragment_step_size')

    settings['min_seen_frag_cutoff'] = cmd_args.min_seen_frag_cutoff
    if 'min_seen_frag_cutoff' in settings_file_args:
        settings['min_seen_frag_cutoff'] = int(settings_file_args['min_seen_frag_cutoff'])
        settings_file_args.pop('min_seen_frag_cutoff')

    settings['fragment_min_mapq'] = cmd_args.fragment_min_mapq
    if 'fragment_min_mapq' in settings_file_args:
        settings['fragment_min_mapq'] = int(settings_file_args['fragment_min_mapq'])
        settings_file_args.pop('fragment_min_mapq')

    settings['dedup_input_on_UMI'] = cmd_args.dedup_input_on_UMI
    if 'dedup_input_on_UMI' in settings_file_args:
        settings['dedup_input_on_UMI'] = (settings_file_args['dedup_input_on_UMI'].lower() == 'true')
        settings_file_args.pop('dedup_input_on_UMI')

    settings['dedup_input_based_on_aln_pos_and_UMI'] = cmd_args.dedup_input_based_on_aln_pos_and_UMI
    if 'dedup_input_based_on_aln_pos_and_UMI' in settings_file_args:
        settings['dedup_input_based_on_aln_pos_and_UMI'] = (settings_file_args['dedup_input_based_on_aln_pos_and_UMI'].lower() == 'true')
        settings_file_args.pop('dedup_input_based_on_aln_pos_and_UMI')

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
                next
            if setting not in previous_settings:
                can_use_previous_analysis = False
                logging.info(('Not using previous analyses - got new setting %s (%s)')%(setting,settings[setting]))
                break
            elif str(settings[setting]) != str(previous_settings[setting]):
                can_use_previous_analysis = False
                logging.info(('Not using previous analyses - setting for %s has changed (%s (new) vs %s (old))')%(setting,settings[setting],previous_settings[setting]))
                break

    if can_use_previous_analysis:
        logging.info('Repeated settings detected. Using previous analyses if completed.')

    with open (settings_used_output_file,'w') as fout:
        for setting in settings:
            fout.write("%s\t%s\n"%(str(setting),str(settings[setting])))

    settings['can_use_previous_analysis'] = can_use_previous_analysis

    return settings


def assert_dependencies(samtools_command='samtools',bowtie2_command='bowtie2',crispresso_command='CRISPResso',casoffinder_command='cas-offinder'):
    """
    Asserts the presence of required software (faidx, bowtie2, casoffinder)

    params:
        samtools_command: location of samtools to run
        bowtie2_command: location of bowtie2 to run
        crispresso_command: location of crispresso to run
        casoffinder_command: location of casoffinder to run

    Raises exception if any command is not found
    """
    logging.info('Checking dependencies')

    # check faidx
    try:
        faidx_result = subprocess.check_output('%s faidx'%samtools_command, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: samtools faidx is required')
    if not 'Usage: samtools faidx' in str(faidx_result):
        raise Exception('Error: samtools faidx is required')

    #check bowtie2
    try:
        bowtie_result = subprocess.check_output('%s --version'%bowtie2_command, stderr=subprocess.STDOUT,shell=True)
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

    return(res/4)

def get_cut_sites_casoffinder(root,genome,pam,guides,cleavage_offset,num_mismatches,casoffinder_command):
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
    returns:
        casoffinder_cut_sites: list of off-target cleavage positions
    """

    logging.info('Calculating cut positions from guide sequence using Cas-OFFinder')


    #link in genome -- it throws a seg fault if it's in the bowtie directory
    linked_genome = root + '.genome.fa'
    if os.path.exists(linked_genome):
        os.remove(linked_genome)
    os.symlink(genome,linked_genome)


    casoffinder_input_file = root + '.casoffinder_input.txt'
    casoffinder_output_file = root + '.casoffinder_output.txt'
    casoffinder_log_file = root + '.casoffinder.log'
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

    logging.debug(casoffinder_cmd)
    if not os.path.exists(casoffinder_output_file):
        casoffinder_result = subprocess.check_output(casoffinder_cmd,shell=True,stderr=subprocess.STDOUT)
        logging.debug('Casoffinder output:' + casoffinder_result)
    else:
        logging.debug('Using previously-calculated offtargets')

    os.remove(linked_genome)

    casoffinder_cut_sites = []
    with open (casoffinder_output_file,'r') as cin:
        for line in cin:
            line_els = line.split("\t")
            chrom = line_els[1]
            loc = line_els[2]
            strand = line_els[4]

            if strand == "+":
                cut_loc = int(loc) + guide_len + 1 + cleavage_offset
            else: #neg strand
                cut_loc = int(loc) - cleavage_offset + 1 + pam_len
            casoffinder_cut_sites.append(chrom+":"+str(cut_loc))
    logging.info('Found %d cut sites'%(len(casoffinder_cut_sites)))

    cut_log = root + '.cut_sites.txt'
    with open (cut_log,'w') as fout:
        fout.write(' '.join(casoffinder_cut_sites))
    logging.info('Wrote cut sites to ' + cut_log)
    return casoffinder_cut_sites

def make_artificial_targets(root,cuts,cut_annotations,genome,target_length,target_padding,primer_chr ="",primer_loc=-1,primer_seq=None,add_non_primer_cut_targets=False,samtools_command='samtools',bowtie2_command='bowtie2',bowtie2_threads=1):
    """
    Generates fasta sequences surrounding cuts for alignment and bowtie2 index
    At each cut point, sequence of length target_length is generated
    Combinations of sequences at each cut site are produced
    If a primer is given, only cut sites that include the primer are used

    params:
        cuts: array of cut locations
        cut_annotations: dict of cut_site->annotation for description of cut (e.g. either On-target, Off-target, Known, Casoffinder, etc)
        genome: location of fasta genome
        target_length: how long the query fragment should be (on one side of the cut) (not including padding). If the primer_seq is given, targets with primer_seq may be shorter than read_length.
        taget_padding: sequence (bp) padding around target (no padding for primer_seq).
        primer_chr: chromosome of primer (if any)
        primer_loc: location of primer (if any)
        primer_seq: sequence of primer binding (if any) (e.g. for dsODN integration and priming). Normally either primer_seq or (primer_chr and primer_loc) are given
        add_non_primer_cut_targets: boolean for whether to add targets for cuts without a primer. If false, only primer-cut1 pairings will be added. If true, cut1-cut2 pairings will be added.
        samtools_command: location of samtools to run
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with

    returns:
        custom_index_fasta: fasta of artificial targets
        target_names: array of target names (corresponding to targets)
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targets)
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut1_anno']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
            target_info[target_name]['cut2_anno']
            target_info[target_name]['query_pos']: genomic start of query (bp)
            target_info[target_name]['target_cut_idx']: bp of cut in target
            target_info[target_name]['target_cut_str']: string designating the cut site, as well as the direction each side of the read comes off of the cut site e.g. w-chr1:50_chr2:60+c means that the left part of the read started on the left side (designated by the -) watson-strand of chr1:50, and the right part of the read started at chr2:60 and extended right on the crick-strand (complement of reference)
    """
    logging.info('Making artificial targets')
    target_names = []
    target_info = {}
    info_file = root + '.info'
    if os.path.isfile(info_file):
        target_count = -1
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 3:
                (custom_index_fasta,target_count_str,target_list_file) = line_els
                target_count = int(target_count_str)
                read_target_count = 0
                if target_count > 0 and os.path.isfile(target_list_file):
                    target_cut_str_index = 0
                    with open(target_list_file,'r') as fin:
                        head_line = fin.readline().rstrip('\n')
                        head_line_els = head_line.split("\t")
                        for line in fin:
                            read_target_count += 1
                            (target_name,target_cut_str,target_class,cut1_chr,cut1_site_str,cut1_anno,cut2_chr,cut2_site_str,cut2_anno,query_pos,target_cut_idx,sequence) = line.rstrip('\n').split("\t")

                            cut1_site = int(cut1_site_str)
                            cut2_site = int(cut2_site_str)

                            target_names.append(target_name)
                            target_info[target_name] = {
                                'sequence':sequence,
                                'class':target_class,
                                'cut1_chr':cut1_chr,
                                'cut1_site':int(cut1_site),
                                'cut1_anno':cut1_anno,
                                'cut2_chr':cut2_chr,
                                'cut2_site':int(cut2_site),
                                'cut2_anno':cut2_anno,
                                'query_pos':int(query_pos),
                                'target_cut_idx':int(target_cut_idx),
                                'target_cut_str':target_cut_str
                                }
                    if read_target_count == target_count:
                        logging.info('Using ' + str(target_count) + ' previously-created targets')
                        return custom_index_fasta,target_names,target_info
                    else:
                        logging.info('In attempting to recover previuosly-created custom targets, expecting ' + str(target_count) + ' targets, but only read ' + str(read_target_count) + ' from ' + target_list_file + '. Recreating.')
                else:
                    logging.info('Could not recover previously-created targets. Recreating.')

    fasta_cache = {}

    #first, if the primer_seq is given (for dsODN) add targets between itself and all other cuts
    if primer_seq is not None and primer_seq != "":
        left_bit_A = primer_seq

        LALA = left_bit_A + reverse(left_bit_A)
        target_name = 'CRISPRlungo_PP' #primer primer
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': LALA,
                'class': 'Primed',
                'cut1_chr':'Primer',
                'cut1_site':len(primer_seq),
                'cut1_anno':'Primer',
                'cut2_chr':'Primer',
                'cut2_site':len(primer_seq),
                'cut2_anno':'Primer',
                'query_pos':0,
                'target_cut_idx':len(primer_seq),
                'target_cut_str':"w-%s:%s~%s:%s-w"%('Primer',len(primer_seq),'Primer',len(primer_seq))
                }


        LALAc = left_bit_A + reverse_complement(left_bit_A)
        target_name = 'CRISPRlungo_PPc' #primer primer complement
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': LALAc,
                'class': 'Primed',
                'cut1_chr':'Primer',
                'cut1_site':len(primer_seq),
                'cut1_anno':'Primer',
                'cut2_chr':'Primer',
                'cut2_site':len(primer_seq),
                'cut2_anno':'Primer',
                'query_pos':0,
                'target_cut_idx':len(primer_seq),
                'target_cut_str':"w-%s:%s~%s:%s-c"%('Primer',len(primer_seq),'Primer',len(primer_seq))
                }

        for j,cut in enumerate(cuts):
            chr_els = cut.split(":")
            chr_B = chr_els[0]
            site_B = int(chr_els[1])
            cut_start_B = site_B - (target_length + target_padding)
            cut_start_B_stop = site_B - 1
            cut_end_B = site_B + (target_length + target_padding)

            left_bit_B_key = '%s %s %d %d'%(genome,chr_B,cut_start_B,cut_start_B_stop)
            if left_bit_B_key not in fasta_cache:
                fasta_cache[left_bit_B_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,cut_start_B,cut_start_B_stop),shell=True).decode(sys.stdout.encoding).strip()
            left_bit_B = fasta_cache[left_bit_B_key]

            right_bit_B_key = '%s %s %d %d'%(genome,chr_B,site_B,cut_end_B)
            if right_bit_B_key not in fasta_cache:
                fasta_cache[right_bit_B_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,site_B,cut_end_B),shell=True).decode(sys.stdout.encoding).strip()
            right_bit_B = fasta_cache[right_bit_B_key]

            LARB = left_bit_A + right_bit_B
            target_name = 'CRISPRlungo_P' + 'R' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LARB,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':len(primer_seq),
                    'cut1_anno':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cut],
                    'query_pos':0,
                    'target_cut_idx':len(primer_seq),
                    'target_cut_str':"w-%s:%s~%s:%s+w"%('Primer',len(primer_seq),chr_B,site_B)
                    }

            LARBc = left_bit_A + complement(right_bit_B)
            target_name = 'CRISPRlungo_P' + 'Rc' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LARBc,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':len(primer_seq),
                    'cut1_anno':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cut],
                    'query_pos':0,
                    'target_cut_idx':len(primer_seq),
                    'target_cut_str':"w-%s:%s~%s:%s+c"%('Primer',len(primer_seq),chr_B,site_B)
                    }

            LALB = left_bit_A + reverse(left_bit_B)
            target_name = 'CRISPRlungo_P' + 'L' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LALB,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':len(primer_seq),
                    'cut1_anno':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cut],
                    'query_pos':0,
                    'target_cut_idx':len(primer_seq),
                    'target_cut_str':"w-%s:%s~%s:%s-w"%('Primer',len(primer_seq),chr_B,site_B)
                    }

            LALBc = left_bit_A + reverse_complement(left_bit_B)
            target_name = 'CRISPRlungo_P' + 'Lc' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LALBc,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':len(primer_seq),
                    'cut1_anno':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cut],
                    'query_pos':0,
                    'target_cut_idx':len(primer_seq),
                    'target_cut_str':"w-%s:%s~%s:%s-c"%('Primer',len(primer_seq),chr_B,site_B)
                    }


    #next, add targets for each cut and itself as well as between cuts
    if add_non_primer_cut_targets:
        for i,cut in enumerate(cuts):
            chr_els = cut.split(":")
            chr_A = chr_els[0]
            site_A = int(chr_els[1])
            cut_start_A = site_A - (target_length + target_padding)
            cut_start_A_stop = site_A - 1
            cut_end_A = site_A + (target_length + target_padding)

            primer_is_in_left_bit_A = (primer_chr != "" and primer_chr == chr_A and primer_loc >= cut_start_A and primer_loc <= cut_start_A_stop)
            primer_is_in_right_bit_A = (primer_chr != "" and primer_chr == chr_A and primer_loc >= site_A and primer_loc <= cut_end_A)

            left_bit_A_key = '%s %s %d %d'%(genome,chr_A,cut_start_A,cut_start_A_stop)

            if left_bit_A_key not in fasta_cache:
                fasta_cache[left_bit_A_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_A,cut_start_A,cut_start_A_stop),shell=True).decode(sys.stdout.encoding).strip()
            left_bit_A = fasta_cache[left_bit_A_key]

            right_bit_A_key = '%s %s %d %d'%(genome,chr_A,site_A,cut_end_A)
            if right_bit_A_key not in fasta_cache:
                fasta_cache[right_bit_A_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_A,site_A,cut_end_A),shell=True).decode(sys.stdout.encoding).strip()
            right_bit_A = fasta_cache[right_bit_A_key]

            #only add both sides if primer_chr is not given -- otherwise, every read has to originate on one end from the primer site.
            wt_seq = left_bit_A + right_bit_A
            if primer_chr == "" or primer_is_in_left_bit_A or primer_is_in_right_bit_A:
                target_name = 'CRISPRlungo_WT'+str(i)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': wt_seq,
                        'class': 'Linear',
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut1_anno':cut_annotations[cut],
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'cut2_anno':cut_annotations[cut],
                        'query_pos':cut_start_A,
                        'target_cut_idx':target_padding + target_length,
                        'target_cut_str':"w-%s:%s~%s:%s+w"%(chr_A,site_A,chr_A,site_A)
                        }

            if primer_chr == "" or primer_is_in_left_bit_A:
                LALA = left_bit_A + reverse(left_bit_A)
                target_name = 'CRISPRlungo_L' + str(i) + 'L' + str(i)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': LALA,
                        'class': 'Chimera',
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut1_anno':cut_annotations[cut],
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'cut2_anno':cut_annotations[cut],
                        'query_pos':cut_start_A,
                        'target_cut_idx':target_padding + target_length,
                        'target_cut_str':"w-%s:%s~%s:%s-w"%(chr_A,site_A,chr_A,site_A)
                        }


                LALAc = left_bit_A + reverse_complement(left_bit_A)
                target_name = 'CRISPRlungo_L' + str(i) + 'Lc' + str(i)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': LALAc,
                        'class': 'Chimera',
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut1_anno':cut_annotations[cut],
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'cut2_anno':cut_annotations[cut],
                        'query_pos':cut_start_A,
                        'target_cut_idx':target_padding + target_length,
                        'target_cut_str':"w-%s:%s~%s:%s-c"%(chr_A,site_A,chr_A,site_A)
                        }

            if primer_chr == "" or primer_is_in_right_bit_A:
                RARA = right_bit_A + reverse(right_bit_A)
                target_name = 'CRISPRlungo_R' + str(i) + 'R' + str(i)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': RARA,
                        'class': 'Chimera',
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut1_anno':cut_annotations[cut],
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'cut2_anno':cut_annotations[cut],
                        'query_pos':cut_start_A,
                        'target_cut_idx':target_padding + target_length,
                        'target_cut_str':"w+%s:%s~%s:%s+w"%(chr_A,site_A,chr_A,site_A)
                        }

                RARAc = right_bit_A + reverse_complement(right_bit_A)
                target_name = 'CRISPRlungo_R' + str(i) + 'Rc' + str(i)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': RARAc,
                        'class': 'Chimera',
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut1_anno':cut_annotations[cut],
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'cut2_anno':cut_annotations[cut],
                        'query_pos':cut_start_A,
                        'target_cut_idx':target_padding + target_length,
                        'target_cut_str':"w+%s:%s~%s:%s+c"%(chr_A,site_A,chr_A,site_A)
                        }


            for j in range(i+1,len(cuts)):
                chr_els = cuts[j].split(":")
                chr_B = chr_els[0]
                site_B = int(chr_els[1])
                cut_start_B = site_B - (target_length + target_padding)
                cut_start_B_stop = site_B - 1
                cut_end_B = site_B + (target_length + target_padding)

                primer_is_in_left_bit_B = (primer_chr != "" and primer_chr == chr_B and primer_loc >= cut_start_B and primer_loc <= cut_start_B_stop)
                primer_is_in_right_bit_B = (primer_chr != "" and primer_chr == chr_B and primer_loc >= site_B and primer_loc <= cut_end_B)

                left_bit_B_key = '%s %s %d %d'%(genome,chr_B,cut_start_B,cut_start_B_stop)
                if left_bit_B_key not in fasta_cache:
                    fasta_cache[left_bit_B_key] = subprocess.check_output(
                        '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,cut_start_B,cut_start_B_stop),shell=True).decode(sys.stdout.encoding).strip()
                left_bit_B = fasta_cache[left_bit_B_key]

                right_bit_B_key = '%s %s %d %d'%(genome,chr_B,site_B,cut_end_B)
                if right_bit_B_key not in fasta_cache:
                    fasta_cache[right_bit_B_key] = subprocess.check_output(
                        '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,site_B,cut_end_B),shell=True).decode(sys.stdout.encoding).strip()
                right_bit_B = fasta_cache[right_bit_B_key]


                if primer_chr == "" or primer_is_in_left_bit_A or primer_is_in_right_bit_B:
                    LARB = left_bit_A + right_bit_B
                    target_name = 'CRISPRlungo_L' + str(i) + 'R' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LARB,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut1_anno':cut_annotations[cuts[i]],
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'cut2_anno':cut_annotations[cuts[j]],
                            'query_pos':cut_start_A,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"w-%s:%s~%s:%s+w"%(chr_A,site_A,chr_B,site_B)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                    LARBc = left_bit_A + complement(right_bit_B)
                    target_name = 'CRISPRlungo_L' + str(i) + 'Rc' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LARBc,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut1_anno':cut_annotations[cuts[i]],
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'cut2_anno':cut_annotations[cuts[j]],
                            'query_pos':cut_start_A,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"w-%s:%s~%s:%s+c"%(chr_A,site_A,chr_B,site_B)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                if primer_chr == "" or primer_is_in_left_bit_B or primer_is_in_right_bit_A:
                    LBRA = left_bit_B + right_bit_A
                    target_name = 'CRISPRlungo_L' + str(j) + 'R' + str(i)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LBRA,
                            'cut1_chr':chr_B,
                            'cut1_site':site_B,
                            'cut1_anno':cut_annotations[cuts[j]],
                            'cut2_chr':chr_A,
                            'cut2_site':site_A,
                            'cut2_anno':cut_annotations[cuts[i]],
                            'query_pos':cut_start_B,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"w-%s:%s~%s:%s+w"%(chr_B,site_B,chr_A,site_A)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                    LBRAc = left_bit_B + complement(right_bit_A)
                    target_name = 'CRISPRlungo_L' + str(j) + 'Rc' + str(i)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LBRAc,
                            'cut1_chr':chr_B,
                            'cut1_site':site_B,
                            'cut1_anno':cut_annotations[cuts[j]],
                            'cut2_chr':chr_A,
                            'cut2_site':site_A,
                            'cut2_anno':cut_annotations[cuts[i]],
                            'query_pos':cut_start_B,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"w-%s:%s~%s:%s+c"%(chr_B,site_B,chr_A,site_A)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                if primer_chr == "" or primer_is_in_left_bit_A or primer_is_in_left_bit_B:
                    LALB = left_bit_A + reverse(left_bit_B)
                    target_name = 'CRISPRlungo_L' + str(i) + 'L' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LALB,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut1_anno':cut_annotations[cuts[i]],
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'cut2_anno':cut_annotations[cuts[j]],
                            'query_pos':cut_start_A,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"w-%s:%s~%s:%s-w"%(chr_A,site_A,chr_B,site_B)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large inversion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                    LALBc = left_bit_A + reverse_complement(left_bit_B)
                    target_name = 'CRISPRlungo_L' + str(i) + 'Lc' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LALBc,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut1_anno':cut_annotations[cuts[i]],
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'cut2_anno':cut_annotations[cuts[j]],
                            'query_pos':cut_start_A,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"w-%s:%s~%s:%s-c"%(chr_A,site_A,chr_B,site_B)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large inversion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                if primer_chr == "" or primer_is_in_right_bit_A or primer_is_in_right_bit_B:
                    RARB = reverse(right_bit_A) + right_bit_B
                    target_name = 'CRISPRlungo_R' + str(i) + 'R' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': RARB,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut1_anno':cut_annotations[cuts[i]],
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'cut2_anno':cut_annotations[cuts[j]],
                            'query_pos':cut_start_A,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"w+%s:%s~%s:%s+w"%(chr_A,site_A,chr_B,site_B)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large inversion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                    RARBc = reverse_complement(right_bit_A) + right_bit_B
                    target_name = 'CRISPRlungo_Rc' + str(i) + 'R' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': RARBc,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut1_anno':cut_annotations[cuts[i]],
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'cut2_anno':cut_annotations[cuts[j]],
                            'query_pos':cut_start_A,
                            'target_cut_idx':target_padding + target_length,
                            'target_cut_str':"c+%s:%s~%s:%s+w"%(chr_A,site_A,chr_B,site_B)
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large inversion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'



    custom_index_fasta = root + '.fa'
    logging.info('Printing ' + str(len(target_names)) + ' targets to custom index (' + custom_index_fasta + ')')
    with open(custom_index_fasta,'w') as fout:
        for i in range(len(target_names)):
            fout.write('>'+target_names[i]+'\n'+target_info[target_names[i]]['sequence']+'\n')

    logging.info('Indexing custom targets using ' + bowtie2_command + '-build (' + custom_index_fasta + ')')
    this_command = bowtie2_command + '-build --offrate 3 --threads ' + str(bowtie2_threads) + ' ' + custom_index_fasta + ' ' + custom_index_fasta
    logging.debug('Bowtie build command: ' + this_command)
    index_result = subprocess.check_output(this_command, shell=True,stderr=subprocess.STDOUT)

    target_list_file = root+".txt"
    with open (target_list_file,'w') as fout:
        fout.write('\t'.join(['target_name','target_cut_str','target_class','cut1_chr','cut1_site','cut1_anno','cut2_chr','cut2_site','cut2_anno','query_pos','target_cut_idx','sequence'])+"\n")
        for target_name in target_names:
            fout.write(target_name+'\t'+'\t'.join([str(target_info[target_name][x]) for x in ['target_cut_str','class','cut1_chr','cut1_site','cut1_anno','cut2_chr','cut2_site','cut2_anno','query_pos','target_cut_idx','sequence']])+"\n")


    #write info file for restarting
    with open(info_file,'w') as fout:
        fout.write("\t".join([str(x) for x in ["custom_index_fasta","target_count","target_list_file"]])+"\n")
        fout.write("\t".join([str(x) for x in [custom_index_fasta,str(len(target_names)),target_list_file]])+"\n")
    return(custom_index_fasta,target_names,target_info)


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
        logging.info('Using previously-deduplicated fastqs with %d reads'%post_dedup_count)
        return(fastq_r1_dedup,fastq_r2_dedup, tot_read_count, count_with_regex, post_dedup_count, post_dedup_read_count)

    #otherwise perform deduplication (check if we were able to read in stats as well -- if we couldn't read them in, tot_read_count will be -1
    logging.info('Deduplicating input fastqs')
    
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
        f1_in = io.BufferedReader(gzip.open(fastq_r1,'rt'))
    else:
        f1_in = open(fastq_r1,'rt')

    if fastq_r2.endswith('.gz'):
        f2_in = io.BufferedReader(gzip.open(fastq_r2,'rt'))
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
#    print("umi_list: " + str(umi_list))

    logging.info('Read %d reads'%tot_read_count)
    logging.info("Processed " + str(umi_count) + " UMIs")
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
        logging.info('Wrote UMI counts to ' + root+".umiCounts.txt")


    logging.info('Wrote %d deduplicated reads'%post_dedup_count)
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
                logging.info('Using previously-generated fastqs with UMIs')
                (fastq_r1_dedup,fastq_r2_dedup,tot_read_count_str) = line_els
                tot_read_count = int(tot_read_count_str)

    #otherwise perform umi adding (check if we were able to read in stats as well -- if we couldn't read them in, tot_read_count will be -1
    if tot_read_count == -1:
        logging.info('Adding UMIs to input fastqs')

        tot_read_count = 0

        if fastq_r1.endswith('.gz'):
            f1_in = io.BufferedReader(gzip.open(fastq_r1,'rb'))
        else:
            f1_in = open(fastq_r1,'r')

        if fastq_r2.endswith('.gz'):
            f2_in = io.BufferedReader(gzip.open(fastq_r2,'rb'))
        else:
            f2_in = open(fastq_r2,'r')

        if fastq_umi.endswith('.gz'):
            umi_in = io.BufferedReader(gzip.open(fastq_umi,'rb'))
        else:
            umi_in = open(fastq_umi,'r')

        fastq_r1_dedup = root + '.r1.gz'
        f1_out = gzip.open(fastq_r1_dedup, 'wb')
        fastq_r2_dedup = root + '.r2.gz'
        f2_out = gzip.open(fastq_r2_dedup, 'wb')

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

        logging.info('Added UMIs fo %d reads'%tot_read_count)
        with open(dedup_stats_file,'w') as fout:
            fout.write("\t".join(["fastq_r1_dedup","fastq_r2_dedup","tot_read_count"])+"\n")
            fout.write("\t".join([str(x) for x in [fastq_r1_dedup,fastq_r2_dedup,tot_read_count]])+"\n")

    #done processing just plotting now
    return(fastq_r1_dedup,fastq_r2_dedup)

def hamming_distance(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def filter_on_primer_seq(root,fastq_r1,fastq_r2,primer_seq,primer_min_bp_match_pct=0,primer_post_seq=None,primer_post_min_bp_match_pct=0,can_use_previous_analysis=False):
    """
    Filters input reads based on whether they contain the primer (and post-primer sequence) with sufficient matching.
    The primer sequence are the bases that should match the primer. The post-primer sequence is the target sequence immediately following the primer and makes sure that the correct genomic region was amplified.

    params:
        root: root for written files
        fastq_r1: R1 reads to filter
        fastq_r2: R2 reads to filter
        primer_seq: sequence of primer to check at beginning of R1
        primer_min_bp_match_pct: min percent of bases that must match primer_seq. If lower than this percent of bases match, the read will be filtered from output files.
        primer_post_seq: sequence immedately following primer sequence to check at beginning of R1
        primer_post_min_bp_match_pct: min percent of bases that must match primer_post_seq. If lower than this percent is matched, the read will be filtered from output files.
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        fastq_r1_filtered: fastq_r1 file of filtered reads
        fastq_r2_filtered: fastq_r2 file of filtered reads
        tot_read_count: number of total reads read
        count_with_primer_match: number of reads with at least primer_min_bp_match_pct matching the primer
        count_with_primer_post_match: number of reads with at least primer_post_min_bp_match_pct matching the primer_post sequence
        count_post_filtering: number of reads post filtering
        filter_on_primer_plot_obj: plot object showing number of reads filtered
    """

    filter_stats_file = root + ".info"
    fastq_r1_filtered = None
    fastq_r2_filtered = None
    tot_read_count = -1
    count_with_primer_match = -1
    count_with_primer_post_match = -1
    count_post_filtering = -1
    #if already finished, attempt to read in stats
    if os.path.isfile(filter_stats_file) and can_use_previous_analysis:
        with open(filter_stats_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) > 3:
                (fastq_r1_filtered_str,fastq_r2_filtered_str,tot_read_count_str,count_with_primer_match_str,count_with_primer_post_match_str,count_post_filtering_str) = line_els
                fastq_r1_filtered = None if fastq_r1_filtered_str == "None" else fastq_r1_filtered_str
                fastq_r2_filtered = None if fastq_r2_filtered_str == "None" else fastq_r2_filtered_str
                tot_read_count = int(tot_read_count_str)
                count_with_primer_match = int(count_with_primer_match_str)
                count_with_primer_post_match = int(count_with_primer_post_match_str)
                count_post_filtering = int(count_post_filtering_str)

            filter_on_primer_plot_obj_str = fin.readline().rstrip('\n')
            filter_on_primer_plot_obj = None
            if filter_on_primer_plot_obj_str != "" and filter_on_primer_plot_obj_str != "None":
                filter_on_primer_plot_obj = PlotObject.from_json(filter_on_primer_plot_obj_str)

        if tot_read_count != -1:
            logging.info('Using previously-filtered fastqs with %d reads'%count_post_filtering)
            return(fastq_r1_filtered, fastq_r2_filtered, tot_read_count, count_with_primer_match, count_with_primer_post_match, count_post_filtering,filter_on_primer_plot_obj)
        else:
            logging.info('Could not recover previously-filtered files. Reanalyzing.')

    #otherwise perform filtering
    logging.info('Filtering input fastqs on primer sequence')

    tot_read_count = 0
    count_with_primer_match = 0
    count_with_primer_post_match = 0
    count_post_filtering = 0

    primer_len = len(primer_seq)
    max_primer_mismatch_count = primer_len - (primer_min_bp_match_pct/float(100) * primer_len)

    primer_post_len = 0
    max_primer_post_mismatch_count = 0
    if primer_post_seq is not None:
        primer_post_len = len(primer_post_seq)
        max_primer_post_mismatch_count = primer_post_len - (primer_post_min_bp_match_pct/float(100) * primer_post_len)

    min_primer_read_len = primer_len + primer_post_len

    if fastq_r1.endswith('.gz'):
        f1_in = io.BufferedReader(gzip.open(fastq_r1,'rt'))
    else:
        f1_in = open(fastq_r1,'rt')

    fastq_r1_filtered = root + '.r1.gz'
    f1_out = gzip.open(fastq_r1_filtered, 'wt')

    #r1 only (single end)
    if fastq_r2 is None:
        while (1):
            f1_id_line   = f1_in.readline()
            f1_seq_line  = f1_in.readline().strip()
            f1_plus_line = f1_in.readline()
            f1_qual_line = f1_in.readline()

            if not f1_qual_line : break
            if not f1_plus_line.startswith('+'):
                raise Exception("Fastq %s cannot be parsed (%s%s%s%s) "%(fastq_r1,f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))
            tot_read_count += 1

            read_is_good = True

            this_mismatches = hamming_distance(primer_seq,f1_seq_line[0:primer_len])
            if this_mismatches > max_primer_mismatch_count:
                read_is_good = False
            else:
                count_with_primer_match += 1

            if primer_post_len > 0:
                this_mismatches = hamming_distance(primer_post_seq,f1_seq_line[primer_len:primer_post_len])
                if this_mismatches > max_primer_post_mismatch_count:
                    read_is_good = False
                else:
                    count_with_primer_post_match += 1

            if read_is_good:
                f1_out.write("%s%s\n%s%s"%(f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))
                count_post_filtering += 1

    #paired end
    else:
        if fastq_r2.endswith('.gz'):
            f2_in = io.BufferedReader(gzip.open(fastq_r2,'rt'))
        else:
            f2_in = open(fastq_r2,'rt')

        fastq_r2_filtered = root + '.r2.gz'
        f2_out = gzip.open(fastq_r2_filtered, 'wt')

        while (1):
            f1_id_line   = f1_in.readline()
            f1_seq_line  = f1_in.readline().strip()
            f1_plus_line = f1_in.readline()
            f1_qual_line = f1_in.readline()

            if not f1_qual_line : break
            if not f1_plus_line.startswith('+'):
                raise Exception("Fastq %s cannot be parsed (%s%s%s%s) "%(fastq_r1,f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))

            f2_id_line   = f2_in.readline()
            f2_seq_line  = f2_in.readline().strip()
            f2_plus_line = f2_in.readline()
            f2_qual_line = f2_in.readline()

            tot_read_count += 1

            read_is_good = True

            this_mismatches = hamming_distance(primer_seq,f1_seq_line[0:primer_len])
            if this_mismatches > max_primer_mismatch_count:
                read_is_good = False
            else:
                count_with_primer_match += 1

            if primer_post_len > 0:
                this_mismatches = hamming_distance(primer_post_seq,f1_seq_line[primer_len:primer_post_len])
                if this_mismatches > max_primer_post_mismatch_count:
                    read_is_good = False
                else:
                    count_with_primer_post_match += 1

            if read_is_good:
                f1_out.write("%s%s\n%s%s"%(f1_id_line,f1_seq_line,f1_plus_line,f1_qual_line))
                f2_out.write("%s%s\n%s%s"%(f2_id_line,f2_seq_line,f2_plus_line,f2_qual_line))
                count_post_filtering += 1

        f2_in.close()
        f2_out.close()

    f1_out.close()
    f1_in.close()

    logging.info('Read %d reads'%tot_read_count)
    logging.info('Wrote %d filtered reads'%count_post_filtering)

    #plot summary
    filter_labels = ['Discarded','Retained']
    values = [tot_read_count-count_post_filtering,count_post_filtering]
    filter_on_primer_plot_obj_root = root + ".filter_on_primer_counts"
    with open(filter_on_primer_plot_obj_root+".txt",'w') as summary:
        summary.write("\t".join(filter_labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    pie_values = []
    pie_labels = []
    for i in range(len(source_labels)):
        if values[i] > 0:
            pie_values.append(values[i])
            pie_labels.append(source_labels[i]+"\n("+str(values[i])+")")
    ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
    ax.set_title('Read filtering based on primer presence')
    plt.savefig(filter_on_primer_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(filter_on_primer_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

    plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(filter_labels,values)])
    plot_label = "Reads were filtered to contain primer sequence '" + primer_seq +"' (with "+str(primer_min_bp_match_pct)+"% bases matching)"
    if primer_post_len > 0:
        plot_label += " and post-primer sequence '" + primer_post_seq + "' (with "+str(primer_post_min_bp_match_pct)+"% bases matching)"

    filter_on_primer_plot_obj = PlotObject(
            plot_name = filter_on_primer_plot_obj_root,
            plot_title = 'Read filtering based on primer presence',
            plot_label = plot_label + '<br>'+plot_count_str,
            plot_datas = [
                ('Primer filter summary',filter_on_primer_plot_obj_root + ".txt")
                ]
            )


    with open(filter_stats_file,'w') as fout:
        fout.write("\t".join(["fastq_r1_filtered","fastq_r2_filtered","tot_read_count","count_with_primer_match","count_with_primer_post_match","count_post_filtering"])+"\n")
        fout.write("\t".join([str(x) for x in [fastq_r1_filtered, fastq_r2_filtered, tot_read_count, count_with_primer_match, count_with_primer_post_match, count_post_filtering]])+"\n")

        filter_on_primer_plot_obj_str = "None"
        if filter_on_primer_plot_obj is not None:
            filter_on_primer_plot_obj_str = filter_on_primer_plot_obj.to_json()
        fout.write(filter_on_primer_plot_obj_str+"\n")

    return(fastq_r1_filtered, fastq_r2_filtered, tot_read_count, count_with_primer_match, count_with_primer_post_match, count_post_filtering, filter_on_primer_plot_obj)


def reverse_complement_cut_str(target_cut_str):
    """
    Reverse-complements a cut string of the form 'w-chr1:50~chr2:90+c' to 'w+chr2:90~chr1:50-c'
    Useful for changing cut if a read aligned to the reverse complement of a target

    params:
        target_cut_str: cut string in form 'w-chr1:50~chr2:90+c'
    returns:
        rc_target_cut_str: reverse complement target_cut_str
    """
    (left_whole,right_whole) = target_cut_str.split("~")
    left_strand = left_whole[0]
    left_dir = left_whole[1]
    left_pos = left_whole[2:]

    right_strand = right_whole[-1]
    right_dir = right_whole[-2:-1]
    right_pos = right_whole[0:-2]

    new_left_strand = 'w' if right_strand == 'c' else 'c'
    new_right_strand = 'w' if left_strand == 'c' else 'c'

    return('%s%s%s~%s%s%s'%(new_left_strand,right_dir,right_pos,left_pos,left_dir,new_right_strand))
    

def align_reads(root,fastq_r1,fastq_r2,reads_name,bowtie2_reference,reference_name,target_info,arm_min_seen_bases=10,arm_min_matched_start_bases=5,bowtie2_command='bowtie2',bowtie2_threads=1,samtools_command='samtools',keep_intermediate=False,can_use_previous_analysis=False):
    """
    Aligns reads to the provided reference (either artificial targets or genome)

    params:
        root: root for written files
        fastq_r1: fastq_r1 to align
        fastq_r2: fastq_r2 to align
        reads_name: Name displayed to user about read origin (e.g. 'R1' or 'R2')
        bowtie2_reference: bowtie2 reference to align to (either artificial targets or reference)
        reference_name: Name displayed to user for updates (e.g. 'Genome' or 'Artificial Targets')
        target_info: hash of information for each target_name
        arm_min_seen_bases: number of bases that are required to be seen on each 'side' of artifical targets. E.g. if a artificial target represents a translocation between chr1 and chr2, arm_min_seen_bases would have to be seen on chr1 as well as on chr2 for the read to be counted.
        arm_min_matched_start_bases: number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each 'side' of artifical targets. E.g. if a artificial target represents a translocation between chr1 and chr2, the first arm_min_matched_start_bases of the read would have to match exactly to chr1 and the last arm_min_matched_start_bases of the read would have to match exactly to chr2
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        samtools_command: location of samtools to run
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        r1_assignments_file: file showing locations of assignments for r1
        r2_assignments_file: file showing locations of assignments for r2
        unmapped_fastq_r1_file: fastq_r1 file of reads not aligning
        unmapped_fastq_r2_file: fastq_r2 file of reads not aligning
        aligned_count: number of reads aligned
        mapped_bam_file: aligned reads
        chr_aln_plot_obj: plot object showing chr locations of alignments
        tlen_plot_object: plot object showing insert sizes for paired alignments
    """

    logging.info('Aligning %s to %s'%(reads_name,reference_name.lower()))

    mapped_bam_file = root + ".bam"

    info_file = root + '.info'
    if os.path.isfile(info_file) and can_use_previous_analysis:
        read_count = -1
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 11:
                (read_count_str,r1_count_str,r2_count_str,unmapped_r1_count_str,unmapped_r2_count_str,aligned_count_str,r1_assignments_file_str,r2_assignments_file_str,unmapped_fastq_r1_file_str,unmapped_fastq_r2_file_str,mapped_bam_file_str) = line_els
                r1_assignments_file = None if r1_assignments_file_str == "None" else r1_assignments_file_str
                r2_assignments_file = None if r2_assignments_file_str == "None" else r2_assignments_file_str
                r1_count = int(r1_count_str)
                r2_count = int(r2_count_str)
                unmapped_r1_count = int(unmapped_r1_count_str)
                unmapped_r2_count = int(unmapped_r2_count_str)
                unmapped_fastq_r1_file = None if unmapped_fastq_r1_file_str == "None" else unmapped_fastq_r1_file_str
                unmapped_fastq_r2_file = None if unmapped_fastq_r2_file_str == "None" else unmapped_fastq_r2_file_str
                read_count  = int(read_count_str)
                aligned_count  = int(aligned_count_str)
                mapped_bam_file = None if mapped_bam_file == "None" else mapped_bam_file_str

                chr_aln_plot_obj_str = fin.readline().rstrip('\n')
                chr_aln_plot_obj = None
                if chr_aln_plot_obj_str != "" and chr_aln_plot_obj_str != "None":
                    chr_aln_plot_obj = PlotObject.from_json(chr_aln_plot_obj_str)

                tlen_plot_obj_str = fin.readline().rstrip('\n')
                tlen_plot_obj = None
                if tlen_plot_obj_str != "" and tlen_plot_obj_str != "None":
                    tlen_plot_obj = PlotObject.from_json(tlen_plot_obj_str)
                if read_count > -1:
                    r1_aln_count = r1_count-unmapped_r1_count
                    r2_aln_count = r2_count-unmapped_r2_count
                    if r2_count > -1:
                        logging.info('Using previously-processed alignment of %d/%d R1 and %d/%d R2 %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,r2_aln_count,r2_count,reads_name,read_count,reference_name.lower()))
                    else:
                        logging.info('Using previously-processed alignment of %d/%d %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,reads_name,read_count,reference_name.lower()))
                    return r1_assignments_file,r2_assignments_file,unmapped_fastq_r1_file, unmapped_fastq_r2_file, aligned_count, mapped_bam_file,chr_aln_plot_obj,tlen_plot_obj
                else:
                    logging.info('Could not recover previously-analyzed alignments. Reanalyzing.')


    bowtie_log = root + '.bowtie2Log'
    if fastq_r2 is not None: #paired-end reads
        logging.info('Aligning paired reads to %s using %s'%(reference_name.lower(),bowtie2_command))
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --end-to-end --threads {bowtie2_threads} -x {bowtie2_reference} -1 {fastq_r1} -2 {fastq_r2} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logging.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        if 'error' in aln_result.lower():
            logging.error('Error found while running command:\n'+aln_command+"\nOutput: "+aln_result)
            raise Exception('Alignment error: ' + aln_result)

        logging.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning paired reads to %s\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(reference_name.lower(),aln_command,aln_result))
    #unpaired reads
    else:
        logging.info('Aligning single-end reads to %s using %s'%(reference_name.lower(),bowtie2_command))
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --end-to-end --threads {bowtie2_threads} -x {bowtie2_reference} -U {fastq_r1} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logging.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        if 'error' in aln_result.lower():
            logging.error('Error found while running command:\n'+aln_command+"\nOutput: "+aln_result)
            raise Exception('Alignment error: ' + aln_result)

        logging.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning single-end reads to %s\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(reference_name.lower(),aln_command,aln_result))

    logging.info('Analyzing reads aligned to %s'%(reference_name.lower()))

    #analyze alignments
    # - if read aligned, store read id and loc in dict r1_assignments
    # - if read unaligned, write read to unaligned fastq
    mapped_tlens = defaultdict(int) # observed fragment lengths
    mapped_chrs = {} #stores the counts aligning to chrs or to templates
    aligned_locs = {} # aligned reads will not be chopped, but keep track of where they aligned for CRISPResso output
    aligned_chr_counts = defaultdict(int)

    unmapped_fastq_r1_file = root + '.unmapped.fq'
    unmapped_fastq_r2_file = None
    r1_assignments_file = root + '.assignments.txt'
    r2_assignments_file = None
    if fastq_r2 is not None: #paired-end reads
        r1_assignments_file = root + '.assignments_r1.txt'
        unmapped_fastq_r1_file = root + '.unmapped_r1.fq'
        r2_assignments_file = root + '.assignments_r2.txt'
        unmapped_fastq_r2_file = root + '.unmapped_r2.fq'
        uf2 = open(unmapped_fastq_r2_file,'w')
        af2 = open(r2_assignments_file,'w')
    uf1 = open(unmapped_fastq_r1_file,'w')
    af1 = open(r1_assignments_file,'w')

    read_count = 0
    r1_count = 0
    r2_count = 0
    unmapped_r1_count = 0
    unmapped_r2_count = 0

    #-F 256 - not primary aligment (filter secondary alignments)
    for line in read_command_output('%s view -F 256 %s'%(samtools_command,mapped_bam_file)):
        if line.strip() == "": break

        line_els = line.split("\t")
        read_count += 1

        read_has_multiple_segments = int(line_els[1]) & 0x1
        read_is_paired_read_1 = int(line_els[1]) & 0x40 # only set with bowtie if read is from paired reads

        read_is_r1 = True
        if read_has_multiple_segments and not read_is_paired_read_1:
            r2_count += 1
            read_is_r1 = False
        else:
            r1_count += 1

        read_id = line_els[0]


        is_bad_alignment = False
        line_unmapped = int(line_els[1]) & 0x4

        left_matches = 0
        right_matches = 0
        mez_val = ""
        for line_el in line_els:
            if line_el.startswith('MD:Z'):
                mez_val = line_el[5:]
                left_matches,right_matches = getLeftRightMismatchesMDZ(mez_val)
                break

        seq = line_els[9]
        is_rc = int(line_els[1]) & 0x10

        if line_unmapped or \
                left_matches < arm_min_matched_start_bases or \
                right_matches < arm_min_matched_start_bases:

            qual = line_els[10]

            if is_rc:
                seq = reverse_complement(seq)

            if read_is_r1:
                uf1.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                unmapped_r1_count += 1
            else:
                uf2.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                unmapped_r2_count += 1
            continue

        line_chr = line_els[2]
        line_mapq = line_els[5]
        line_start = int(line_els[3])
        line_end = line_start+(len(seq)-1)
        if is_rc:
            tmp = line_end
            line_end = line_start
            line_start = tmp

        curr_classification = 'Linear'
        curr_position = "%s:%s-%s"%(line_chr,line_start,line_end)
        curr_annotation = ''
        curr_cut = 'NA'
        if line_chr in target_info: #if this aligned to a custom chromosome
            (left_read_bases_count, left_ref_bases_count, left_all_match_count, left_start_match_count, right_read_bases_count, right_ref_bases_count, right_all_match_count, right_start_match_count) = getMatchLeftRightOfCut(target_info[line_chr]['sequence'],line,target_info[line_chr]['target_cut_idx'])

            #if read doesn't sufficiently align to both parts of the artificial target, print it to unmapped and continue
            if left_read_bases_count < arm_min_seen_bases or right_read_bases_count < arm_min_seen_bases or left_start_match_count < arm_min_matched_start_bases or right_start_match_count < arm_min_matched_start_bases:

                qual = line_els[10]

                if is_rc:
                    seq = reverse_complement(seq)

                if read_is_r1:
                    uf1.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                    unmapped_r1_count += 1
                else:
                    uf2.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                    unmapped_r2_count += 1
                continue

            curr_classification = target_info[line_chr]['class']
            curr_annotation = target_info[line_chr]['target_cut_str'] + ' ' + line_chr
            if is_rc:
                curr_annotation = reverse_complement_cut_str(target_info[line_chr]['target_cut_str']) + ' ' + line_chr

            cut1 =  target_info[line_chr]['cut1_chr'] + ":" +  str(target_info[line_chr]['cut1_site'])
            cut2 =  target_info[line_chr]['cut2_chr'] + ":" +  str(target_info[line_chr]['cut2_site'])

            left_dir = target_info[line_chr]['target_cut_str'][1:2]
            if left_dir == '-': #arm extends left in genome space after cut
                left_arm_end = target_info[line_chr]['cut1_site']
                left_arm_start = left_arm_end - left_ref_bases_count
            else: #arm extends right in genome space after cut
                left_arm_start = target_info[line_chr]['cut1_site']
                left_arm_end = left_arm_start + left_ref_bases_count

            right_dir = target_info[line_chr]['target_cut_str'][-2:-1]
            if right_dir == '+': #arm extends right in genome space after cut
                right_arm_start = target_info[line_chr]['cut2_site']
                right_arm_end = right_arm_start + right_ref_bases_count
            else:
                right_arm_end = target_info[line_chr]['cut2_site']
                right_arm_start = right_arm_end - right_ref_bases_count


            if curr_classification == 'Linear':
                #if read aligned to custom seqs in the forward direction
                if not is_rc:
                    curr_cut = cut1
                    curr_position = '%s:%s-%s'%(
                            target_info[line_chr]['cut1_chr'],left_arm_start,right_arm_end - 1
                            )
                else:
                    curr_cut = cut1
                    curr_position = '%s:%s-%s'%(
                            target_info[line_chr]['cut1_chr'],right_arm_end - 1,left_arm_start
                            )

            #if custom ref is not linear, set cuts and genome aln locations
            else:
                #if read aligned to custom seqs in the forward direction
                if not is_rc:
                    curr_cut = cut1 + "~" + cut2
                    curr_position = '%s:%s-%s~%s:%s-%s'%(
                            target_info[line_chr]['cut1_chr'],left_arm_start,left_arm_end,
                            target_info[line_chr]['cut2_chr'],right_arm_start,right_arm_end - 1
                            )
                #if read aligned to custom seqs in the reverse direction
                else:
                    curr_cut = cut2 + "~" + cut1
                    curr_position = '%s:%s-%s~%s:%s-%s'%(
                            target_info[line_chr]['cut2_chr'],right_arm_end - 1,right_arm_start,
                            target_info[line_chr]['cut1_chr'],left_arm_end,left_arm_start
                            )

        if read_is_r1:
            af1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(read_id,reference_name,curr_classification,curr_annotation,curr_position,curr_cut))
        else:
            af2.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(read_id,reference_name,curr_classification,curr_annotation,curr_position,curr_cut))

        if read_is_paired_read_1: #only set if paired
            insert_size = int(line_els[8])
            mapped_tlens[insert_size] += 1

        if not line_chr in mapped_chrs:
            mapped_chrs[line_chr] = 0
        mapped_chrs[line_chr] += 1

        if line_chr not in aligned_locs:
            aligned_locs[line_chr] = {}
        if line_start not in aligned_locs[line_chr]:
            aligned_locs[line_chr][line_start] = 0

        aligned_locs[line_chr][line_start] += 1
        aligned_chr_counts[line_chr] += 1


    uf1.close()
    af1.close()
    if unmapped_fastq_r2_file is not None:
        uf2.close()
        af2.close()

    chr_aln_plot_root = root + ".chrs"
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
    ax.set_title('Location of ' + reads_name + ' aligned to the '+reference_name.lower())
    plt.savefig(chr_aln_plot_root+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(chr_aln_plot_root+".png",pad_inches=1,bbox_inches='tight')

    plot_label = 'Bar plot showing alignment location of ' + reads_name + ' aligned to ' + reference_name.lower()
    if unmapped_fastq_r2_file is not None:
        plot_label += ' (Note that this plot shows the number of reads not readpairs)'
    if len(keys) == 0:
        plot_label = '(No ' + reads_name + ' aligned to ' + reference_name.lower() + ')'

    chr_aln_plot_obj = PlotObject(
            plot_name = chr_aln_plot_root,
            plot_title = reference_name.title() + ' alignment summary',
            plot_label = plot_label,
            plot_datas = [(reference_name.capitalize() + ' alignment summary',chr_aln_plot_root + ".txt")]
            )

    tlen_plot_obj = None
    if fastq_r2 is not None: #paired-end reads
        tlen_plot_root = root + ".insertSizes"
        keys = sorted(mapped_tlens.keys())
        vals = [mapped_tlens[key] for key in keys]
        with open(tlen_plot_root+".txt","w") as fout:
            fout.write('insertSize\tnumReads\n')
            for key in keys:
                fout.write(str(key) + '\t' + str(mapped_tlens[key]) + '\n')

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
                plot_title = reference_name.capitalize() + ' Alignment Insert Size Summary',
                plot_label = 'Bar plot showing insert size of reads aligned to ' + reference_name,
                plot_datas = [(reference_name.capitalize() + ' alignment insert size summary',tlen_plot_root + ".txt")]
                )

    aligned_count = read_count - (unmapped_r1_count + unmapped_r2_count)

    with open(info_file,'w') as fout:
        fout.write("\t".join([str(x) for x in ["read_count","r1_count","r2_count","unmapped_r1_count","unmapped_r2_count","aligned_count","r1_assignments_file","r2_assignments_file","unmapped_fastq_r1_file","unmapped_fastq_r2_file","mapped_bam_file"]])+"\n")
        fout.write("\t".join([str(x) for x in [read_count,r1_count,r2_count,unmapped_r1_count,unmapped_r2_count,aligned_count,r1_assignments_file,r2_assignments_file,unmapped_fastq_r1_file,unmapped_fastq_r2_file,mapped_bam_file]])+"\n")

        chr_aln_plot_obj_str = "None"
        if chr_aln_plot_obj is not None:
            chr_aln_plot_obj_str = chr_aln_plot_obj.to_json()
        fout.write(chr_aln_plot_obj_str+"\n")
        tlen_plot_obj_str = "None"
        if tlen_plot_obj is not None:
            tlen_plot_obj_str = tlen_plot_obj.to_json()
        fout.write(tlen_plot_obj_str+"\n")

    r1_aln_count = r1_count-unmapped_r1_count
    r2_aln_count = r2_count-unmapped_r2_count
    if r2_aln_count > 0:
        logging.info('Finished analysis of %d/%d R1 and %d/%d R2 %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,r2_aln_count,r2_count,reads_name,read_count,reference_name.lower()))
    else:
        logging.info('Finished analysis of %d/%d %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,reads_name,read_count,reference_name.lower()))
    return(r1_assignments_file,r2_assignments_file,unmapped_fastq_r1_file,unmapped_fastq_r2_file,aligned_count,mapped_bam_file,chr_aln_plot_obj,tlen_plot_obj)

def make_final_read_assignments(root,r1_assignment_files,r2_assignment_files,cut_sites,cut_annotations,target_info,r1_r2_support_max_distance=100000,cut_merge_dist=20,genome_map_resolution=1000000,dedup_based_on_aln_pos=True,discard_reads_without_r2_support=True,can_use_previous_analysis=False):
    """
    Makes final read assignments after:
        Pairing R1 and R2 if they were processed in separate steps
        Deduplicating based on UMI and alignment location

    params:
        root: root for written files
        r1_assignment_files: list of r1 assignment files
        r2_assignment_files: list of r2 assignment files
        cut_sites: array of known cut sites
        cut_annotations: dict of cut_site->annotation for description of cut (e.g. either On-target, Off-target, Known, Casoffinder, etc)
        target_info: hash of information for each target_name
        r1_r2_support_max_distance: max distance between r1 and r2 for the read pair to be classified as 'supported' by r2
        cut_merge_dist: cuts within this distance (bp) will be merged
        genome_map_resolution: window size (bp) for reporting number of reads aligned
        dedup_based_on_aln_pos: if true, reads with the same UMI and alignment positions will be removed from analysis
        discard_reads_without_r2_support: if true, reads without r2 support will be discarded from final analysis and counts
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        final_assignment_filename: filename of final assignments for each read id pair
        r1_read_ids_for_crispresso: dict of readID=>list of cutIDs for that read
        r1_cut_counts_for_crispresso: dict of cutID=>count how many times each cut was seen -- when we iterate through the reads (fastq) in the next step, we check to see that the cut was seen above a threshold before printing those reads. The count is stored in this dict.
        r2_read_ids_for_crispresso: dict of readID=>cut assignment
        r2_cut_points_for_crispresso
        deduplication_plot_obj: plot showing how many reads were deuplicated
        r1_source_plot_obj: plot for R1 sources (genome/custom/frags)
        r1_classification_plot_obj: plot for R1 assignments (linear, translocation, etc)
        r2_source_plot_obj: plot for R2 sources (genome/custom/frags)
        r2_classification_plot_obj: plot for R2 assignments (linear, translocation, etc)
    """

    final_file = root + '.final_assignments'
    info_file = root + '.info'
    #check to see if this anaylsis has been completed previously
    if os.path.isfile(info_file) and os.path.isfile(final_file) and can_use_previous_analysis:
        read_total_read_count = 0 #how many lines we read in this time by reading the final assignments file
        previous_total_read_count = -1 #how many lines were reported to be read the preivous time (if the analysis was completed previously)
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 5:
                (previous_total_read_count_str,final_duplicate_count_str,final_nonsupporting_count_removed_str, final_processed_count_str, final_observed_cuts_count_str) = line_els
                previous_total_read_count = int(previous_total_read_count_str)
                final_duplicate_count = int(final_duplicate_count_str)
                final_nonsupporting_count_removed = int(final_nonsupporting_count_removed_str)
                final_processed_count = int(final_processed_count_str)
                final_observed_cuts_count = int(final_observed_cuts_count_str)

            #load plots
            deduplication_plot_obj_str = fin.readline().rstrip('\n')
            deduplication_plot_obj = None
            if deduplication_plot_obj_str != "" and deduplication_plot_obj_str != "None":
                deduplication_plot_obj = PlotObject.from_json(deduplication_plot_obj_str)
            r1_source_plot_obj_str = fin.readline().rstrip('\n')
            r1_source_plot_obj = None
            if r1_source_plot_obj_str != "" and r1_source_plot_obj_str != "None":
                r1_source_plot_obj = PlotObject.from_json(r1_source_plot_obj_str)
            r1_classification_plot_obj_str = fin.readline().rstrip('\n')
            r1_classification_plot_obj = None
            if r1_classification_plot_obj_str != "" and r1_classification_plot_obj_str != "None":
                r1_classification_plot_obj = PlotObject.from_json(r1_classification_plot_obj_str)

            r2_source_plot_obj_str = fin.readline().rstrip('\n')
            r2_source_plot_obj = None
            if r2_source_plot_obj_str != "" and r2_source_plot_obj_str != "None":
                r2_source_plot_obj = PlotObject.from_json(r2_source_plot_obj_str)
            r2_classification_plot_obj_str = fin.readline().rstrip('\n')
            r2_classification_plot_obj = None
            if r2_classification_plot_obj_str != "" and r2_classification_plot_obj_str != "None":
                r2_classification_plot_obj = PlotObject.from_json(r2_classification_plot_obj_str)

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
            r1_read_ids_for_crispresso = defaultdict(list)
            r1_cut_counts_for_crispresso = defaultdict(int)
            r2_read_ids_for_crispresso = defaultdict(list)
            r2_cut_counts_for_crispresso = defaultdict(int)
            with open (final_file,'r') as fin:
                head_line = fin.readline().rstrip('\n')
                head_line_els = head_line.split("\t")
                read_total_read_count = 0
                r1_final_cut_ind = 18 #for paired-end
                r2_final_cut_ind = 20 #for paired-end
                single_r1_final_cut_ind = 9 #for single-end
                #make sure header matches
                #paired-end mode
                if len(head_line_els) == 21 and head_line_els[r1_final_cut_ind] == "r1_final_cut_keys" and head_line_els[r2_final_cut_ind] == "r2_final_cut_keys":
                    #if header matches, read file
                    for line in fin:
                        read_total_read_count += 1
                        line_els = line.rstrip("\r\n").split("\t")
                        if line_els[r1_final_cut_ind] != "NA":
                            r1_keys = line_els[r1_final_cut_ind].split(",")
                            r1_read_ids_for_crispresso[line_els[0]] = r1_keys
                            for key in r1_keys:
                                r1_cut_counts_for_crispresso[key] += 1
                        if line_els[r2_final_cut_ind] != "NA":
                            r2_keys = line_els[r2_final_cut_ind].split(",")
                            r2_read_ids_for_crispresso[line_els[0]] = r2_keys
                            for key in r2_keys:
                                r2_cut_counts_for_crispresso[key] += 1

                elif len(head_line_els) == 10 and head_line_els[single_r1_final_cut_ind] == "final_cut_keys":
                    r2_read_ids_for_crispresso = None
                    r2_cut_counts_for_crispresso = None
                    #if header matches, read file
                    for line in fin:
                        read_total_read_count += 1
                        line_els = line.rstrip("\r\n").split("\t")
                        if line_els[single_r1_final_cut_ind] != "NA":
                            r1_keys = line_els[single_r1_final_cut_ind].split(",")
                            r1_read_ids_for_crispresso[line_els[0]] = r1_keys
                            for key in r1_keys:
                                r1_cut_counts_for_crispresso[key] += 1

                if previous_total_read_count == read_total_read_count:
                    logging.info('Using previously-processed assignments for ' + str(final_processed_count) + '/' + str(read_total_read_count) + ' total reads')
                    return final_file,r1_read_ids_for_crispresso,r1_cut_counts_for_crispresso,r2_read_ids_for_crispresso,r2_cut_counts_for_crispresso,deduplication_plot_obj,r1_source_plot_obj,r1_classification_plot_obj,r2_source_plot_obj,r2_classification_plot_obj,r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,r1_r2_no_support_dist_plot_obj
                else:
                    logging.info('In attempting to recover previously-processed assignments, expecting ' + str(previous_total_read_count) + ' reads, but only read ' + str(read_total_read_count) + ' from ' + final_file + '. Reprocessing.')

    logging.info('Making final read assignments')

    sorted_r1_file = root + '.r1.sorted'
    sorted_r2_file = root + '.r2.sorted'
    final_file_tmp = root + '.final_assignments.tmp'
    cat_and_sort_cmd = 'cat ' + " ".join(r1_assignment_files) + ' | LC_COLLATE=C sort > ' + sorted_r1_file
#    cat_and_sort_cmd = 'cat ' + " ".join(r1_assignment_files) + ' | sort > ' + sorted_r1_file
    logging.debug('r1 cat command: ' + str(cat_and_sort_cmd))
    cat_and_sort_result = subprocess.check_output(cat_and_sort_cmd, shell=True,stderr=subprocess.STDOUT)

    #if single end
    if len(r2_assignment_files) == 0:
        sorted_r2_file = None
    else:
        #cat_and_sort_cmd = 'cat ' + " ".join(r2_assignment_files) + ' | sort > ' + sorted_r2_file
        cat_and_sort_cmd = 'cat ' + " ".join(r2_assignment_files) + ' | LC_COLLATE=C sort > ' + sorted_r2_file
        logging.debug('r2 cat command: ' + str(cat_and_sort_cmd))
        cat_and_sort_result = subprocess.check_output(cat_and_sort_cmd, shell=True,stderr=subprocess.STDOUT)

    #first pass through file:
    # -aggregate cut sites that will be later aggregated into final cut sites
    # -identify duplicates based on r1/r2 alignment
    # -classify reads based on r1/r2 support

    id_ind = 0
    source_ind = 1
    classification_ind = 2
    annotation_ind = 3
    alignment_ind = 4
    cut_point_ind = 5

    seen_reads = {} #put into this dict if it's been seen for deduplicating
    total_reads_processed = 0
    dups_seen = 0
    cut_points_by_chr = {} # dict of cut points for each chr contianing observedcuts cut_points_by_chr[chr1][1559988] = count seen

    #these classificaiton strings are shown to the user and can be modified in one place here
    r1_r2_support_str_supported = "Supported by R2"
    r1_r2_support_str_not_supported_orientation = "Not supported by R2 - different R1/R2 orientations"
    r1_r2_support_str_not_supported_distance = "Not supported by R2 - alignment separation more than maximum"
    r1_r2_support_str_not_supported_diff_chrs = "Not supported by R2 - different chrs"
    r1_r2_support_str_not_supported_not_aln = "Not supported by R2 - unaligned"

    if sorted_r2_file is None:
        #single end pre-processing of cut sites
        with open(sorted_r1_file,'r') as sf1, open(final_file_tmp,'w') as ff:
            for r1_line in sf1:
                r1_line = r1_line.rstrip('\n')
                total_reads_processed += 1

                r1_line_els = r1_line.split("\t")
                r1_id = r1_line_els[0]
                r1_umi = r1_line_els[0].split(":")[-1]

                this_key = r1_umi + " " + r1_line_els[alignment_ind]
                #if this is a duplicate, skip it
                if this_key in seen_reads:
                    dups_seen += 1
                    ff.write(r1_line + "\tduplicate\t"+seen_reads[this_key]+"\n")
                    continue

                #if not a duplicate:
                seen_reads[this_key] = r1_id
                ff.write(r1_line + "\tnot_duplicate\tNA\n")

                #parse out cut positions and record them
                r1_cut_points = r1_line_els[cut_point_ind]
                if r1_line_els[classification_ind] != 'Linear':
                    for cut_point in r1_cut_points.split("~"):
                        (cut_chr,cut_pos) = cut_point.split(":")
                        if cut_chr != "*":
                            if cut_chr not in cut_points_by_chr:
                                cut_points_by_chr[cut_chr] = defaultdict(int)
                            cut_points_by_chr[cut_chr][int(cut_pos)] += 1
    else:
        #paired end pre-processing of cut sites
        with open(sorted_r1_file,'r') as sf1, open(sorted_r2_file,'r') as sf2, open(final_file_tmp,'w') as ff:
            for r1_line in sf1:
                r1_line = r1_line.rstrip('\n')
                r2_line = sf2.readline().rstrip('\n')
                total_reads_processed += 1

                r1_line_els = r1_line.split("\t")
                r1_id = r1_line_els[0]
                r1_umi = r1_line_els[0].split(":")[-1]

                r2_line_els = r2_line.split("\t")
                r2_id = r2_line_els[0]
                r2_umi = r2_line_els[0].split(":")[-1]


                if r1_id.split(" ")[0] != r2_id.split(" ")[0]:
                    raise Exception("Some reads were lost -- read ids don't match up:\n1:"+r1_line+"\n2:"+r2_line)

                # assign support for this read
                # these variables will be printed out
                r1_r2_support_status = None
                r1_r2_support_dist = None # distance between R1 and R2
                r2_r2_orientation_str = None

                #these variables are for assignment of this read
                r1_last_aln_chr_for_support = None #values to see if r1 and r2 agree
                r1_last_aln_pos_for_support = None
                r1_last_orientation_for_support = None
                r2_first_aln_chr_for_support = None
                r2_first_aln_pos_for_support = None
                r2_first_orientation_for_support = None

                r1_alignment = r1_line_els[alignment_ind]
                for r1_aln_pos in r1_alignment.split("~"):
                    if '*' not in r1_aln_pos:
                        (aln_chr,aln_range) = r1_aln_pos.split(":")
                        (aln_pos,aln_end) = aln_range.split("-")

                        r1_last_aln_chr_for_support = aln_chr
                        r1_last_aln_pos_for_support = aln_end
                        if (aln_pos < aln_end):
                            r1_last_orientation_for_support = "+"
                        else:
                            r1_last_orientation_for_support = "-"
                    else:
                        r1_last_aln_chr_for_support = "*"
                        r1_last_aln_pos_for_support = "*"
                        r1_last_orientation_for_support = None

                # make assignment for R1/R2 support
                r2_alignment = r2_line_els[alignment_ind]
                r2_aln_pos = r2_alignment.split("~")[0]
                if '*' not in r2_aln_pos:
                    (aln_chr,aln_range) = r2_aln_pos.split(":")
                    (aln_pos,aln_end) = aln_range.split("-")
                    r2_first_aln_chr_for_support = aln_chr
                    r2_first_aln_pos_for_support = aln_pos
                    if (aln_pos < aln_end):
                        r2_first_orientation_for_support = "+"
                    else:
                        r2_first_orientation_for_support = "-"
                else:
                    r2_first_aln_chr_for_support = "*"
                    r2_first_aln_pos_for_support = "*"

                #determine R1/R2 support
                if r1_last_aln_chr_for_support == "*" or r2_first_aln_chr_for_support == "*":
                    r1_r2_support_status = r1_r2_support_str_not_supported_not_aln
                elif r1_last_aln_chr_for_support == r2_first_aln_chr_for_support:
                    r1_r2_support_dist = abs(int(r1_last_aln_pos_for_support) - int(r2_first_aln_pos_for_support))
                    #distance must be within cutoff
                    if r1_r2_support_dist <= r1_r2_support_max_distance:
                        #must have different orientations (R1 is pos, R2 is rc/egative)
                        #for standard illumina sequencing, R1 is forward and R2 is reverse-complement
                        if r1_last_orientation_for_support != r2_first_orientation_for_support:
                            r1_r2_support_status = r1_r2_support_str_supported
                        else:
                            r1_r2_support_status = r1_r2_support_str_not_supported_orientation
                    else:
                        r1_r2_support_status = r1_r2_support_str_not_supported_distance
                else:
                    r1_r2_support_status = r1_r2_support_str_not_supported_diff_chrs

                r1_r2_orientation_str = "%s/%s"%(r1_last_orientation_for_support,r2_first_orientation_for_support)


                this_key = r1_umi + " " + r1_line_els[alignment_ind] + " " + r2_line_els[alignment_ind]
                #if this is a duplicate, skip it
                if this_key in seen_reads:
                    dups_seen += 1
                    ff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(r1_line,r2_line,'duplicate',seen_reads[this_key],r1_r2_support_status,r1_r2_support_dist,r1_r2_orientation_str))
                    continue

                #if not a duplicate:
                seen_reads[this_key] = r1_id
                ff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(r1_line,r2_line,'not_duplicate','NA',r1_r2_support_status,r1_r2_support_dist,r1_r2_orientation_str))


                #parse out cut positions and record them
                r1_cut_points = r1_line_els[cut_point_ind]
                r2_cut_points = r2_line_els[cut_point_ind]
                if r1_line_els[classification_ind] != 'Linear':
                    for cut_point in r1_cut_points.split("~"):
                        (cut_chr,cut_pos) = cut_point.split(":")
                        if cut_chr != "*":
                            if cut_chr not in cut_points_by_chr:
                                cut_points_by_chr[cut_chr] = defaultdict(int)
                            cut_points_by_chr[cut_chr][int(cut_pos)] += 1
                if r2_line_els[classification_ind] != 'Linear':
                    for cut_point in r2_cut_points.split("~"):
                        (cut_chr,cut_pos) = cut_point.split(":")
                        if cut_chr != "*":
                            if cut_chr not in cut_points_by_chr:
                                cut_points_by_chr[cut_chr] = defaultdict(int)
                            cut_points_by_chr[cut_chr][int(cut_pos)] += 1
        #end paired end pre-processing of cut sites

    #now that we've aggregated all cut sites, find likely cut positions and record assignments to those positions in final_cut_point_lookup
    #cut sites are input by user
    known_cut_points_by_chr = {}
    for cut_site in cut_sites:
        (cut_chr,cut_pos) = cut_site.split(":")
        if cut_chr not in known_cut_points_by_chr:
            known_cut_points_by_chr[cut_chr] = []
        known_cut_points_by_chr[cut_chr].append(int(cut_pos))

    #the counts in this file include only reads that weren't marked as duplicates, and reads that weren't Linear
    #each non-duplicate read could contribute at least 2 cut points (perhaps more?), for a total of at least 2*(non-duplicated,non-linear)
    with open(root+'.seen_cut_points','w') as fout:
        for cut_chr in cut_points_by_chr:
            these_cut_points = sorted(cut_points_by_chr[cut_chr].keys())
            for cut_point in these_cut_points:
                fout.write(cut_chr + ':' + str(cut_point)+"\t"+str(cut_points_by_chr[cut_chr][cut_point])+"\n")

    final_cut_points_by_chr = {} #dict of final cut points by chromosome, contains count of reads with that cut point
    final_cut_point_lookup = {} #dict from old(fuzzy/imprecise) position to new
    r1_translocation_cut_counts = {} #dict of counts for translocations between cut1 and cut2 e.g. r1_translcoation_cut_counts[cut1][cut2] = 5
    r1_translocation_cut_counts['*'] = defaultdict(int)
    r2_translocation_cut_counts = {}
    r2_translocation_cut_counts['*'] = defaultdict(int)

    r1_fragment_translocation_cut_counts = {} #dict of counts for translocations as in r1_translocation_cut_counts except the keys have a the watson/crick and direction annotation (e.g. w-chr2:60455)
    r1_fragment_translocation_classifications = {} #classification for each cut
    r1_fragment_translocation_cut_counts['u+*'] = defaultdict(int)
    r2_fragment_translocation_cut_counts = {}
    r2_fragment_translocation_cut_counts['u+*'] = defaultdict(int)


    for cut_chr in cut_points_by_chr:
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
                if cut_chr in known_cut_points_by_chr:
                    for known_cut_point in known_cut_points_by_chr[cut_chr]:
                        this_dist = abs(known_cut_point - this_weighted_mean)
                        if this_dist <= min_dist:
                            this_pos = known_cut_point
                            min_dist = this_dist

                #initialize defaultdict for this final_cut_point
                r1_translocation_cut_counts[cut_chr+':'+str(this_pos)] = defaultdict(int)
                r2_translocation_cut_counts[cut_chr+':'+str(this_pos)] = defaultdict(int)
                r1_fragment_translocation_cut_counts['w-'+cut_chr+':'+str(this_pos)] = defaultdict(int)
                r1_fragment_translocation_cut_counts['c-'+cut_chr+':'+str(this_pos)] = defaultdict(int)
                r1_fragment_translocation_cut_counts['w+'+cut_chr+':'+str(this_pos)] = defaultdict(int)
                r1_fragment_translocation_cut_counts['c+'+cut_chr+':'+str(this_pos)] = defaultdict(int)
                r2_fragment_translocation_cut_counts['w-'+cut_chr+':'+str(this_pos)] = defaultdict(int)
                r2_fragment_translocation_cut_counts['c-'+cut_chr+':'+str(this_pos)] = defaultdict(int)
                r2_fragment_translocation_cut_counts['w+'+cut_chr+':'+str(this_pos)] = defaultdict(int)
                r2_fragment_translocation_cut_counts['c+'+cut_chr+':'+str(this_pos)] = defaultdict(int)

                read_counts_at_curr_point = 0
                #add this assignment to the lookup
                for curr_point in curr_points:
                    final_cut_point_lookup[cut_chr+":"+str(curr_point)] = cut_chr+":"+str(this_pos)
                    read_counts_at_curr_point += cut_points_by_chr[cut_chr][curr_point]

                if cut_chr not in final_cut_points_by_chr:
                    final_cut_points_by_chr[cut_chr] = {}
                final_cut_points_by_chr[cut_chr][this_pos] = read_counts_at_curr_point

                #reset current points
                curr_points = []

            if i < len(these_cut_points):
                curr_points.append(these_cut_points[i])

    with open(root + ".final_cut_points.txt",'w') as fout:
        fout.write('chr\tpos\tcount\tannotation\n')
        for chrom in sorted(final_cut_points_by_chr.keys()):
            for pos in sorted(final_cut_points_by_chr[chrom]):
                this_anno='Novel'
                this_str = chrom+":"+str(pos)
                if this_str in cut_annotations:
                    this_anno = cut_annotations[this_str]
                fout.write("%s\t%d\t%d\t%s\n"%(chrom,pos,final_cut_points_by_chr[chrom][pos],this_anno))

    with open(root+".final_cut_points_lookup.txt",'w') as fout:
        fout.write('discovered cut point\tfinal mapped cut point\tcount\tannotation\n');
        for key in sorted(final_cut_point_lookup.keys(),key=lambda k: (k.split(":")[0],int(k.split(":")[1]))):
            (chrom,pos)=key.split(":")
            this_anno='Novel'
            this_final = final_cut_point_lookup[key]
            if this_final in cut_annotations:
                this_anno = cut_annotations[this_final]
            fout.write("%s\t%s\t%d\t%s\n"%(key,final_cut_point_lookup[key],cut_points_by_chr[chrom][int(pos)],this_anno))


    final_total_count = 0
    final_duplicate_count = 0
    final_nonsupporting_count_removed = 0
    final_processed_count = 0
    r1_final_source_counts= defaultdict(int) # source is Genome, Custom targets or Fragmented
    r1_final_classification_counts = defaultdict(int) # classification is Primed Linear Chimera Translocation or Unidentified
    r1_final_source_classification_counts = defaultdict(int)
    r2_final_source_counts = defaultdict(int)
    r2_final_classification_counts = defaultdict(int)
    r2_final_source_classification_counts = defaultdict(int)


    #keep track of alignments to genome
    r1_aln_pos_counts_by_chr = {}
    r1_aln_pos_genome_counts_by_chr = {}
    r1_aln_pos_custom_aln_counts_by_chr = {}
    r1_aln_pos_frag_counts_by_chr = {}

    r2_aln_pos_counts_by_chr = {}
    r2_aln_pos_genome_counts_by_chr = {}
    r2_aln_pos_custom_aln_counts_by_chr = {}
    r2_aln_pos_frag_counts_by_chr = {}

    #keep track of alignments to custom chrs
    r1_custom_aln_counts = {}

    r1_uncut_counts = defaultdict(int)
    r1_uncut_genome_counts = defaultdict(int)
    r1_uncut_custom_aln_counts = defaultdict(int)
    r1_uncut_frag_counts = defaultdict(int)
    r1_cut_counts = defaultdict(int)
    r1_cut_custom_aln_counts = defaultdict(int)
    r1_cut_frag_counts = defaultdict(int)

    r2_uncut_counts = defaultdict(int)
    r2_uncut_genome_counts = defaultdict(int)
    r2_uncut_custom_aln_counts = defaultdict(int)
    r2_uncut_frag_counts = defaultdict(int)
    r2_cut_counts = defaultdict(int)
    r2_cut_custom_aln_counts = defaultdict(int)
    r2_cut_frag_counts = defaultdict(int)

    #keep track of read ids aligning to cuts
    #reads with these ids will be pulled out of the fastq for analysis by CRISPResso
    r1_read_ids_for_crispresso = defaultdict(list)
    r1_cut_counts_for_crispresso = defaultdict(int)
    r2_read_ids_for_crispresso = defaultdict(list)
    r2_cut_counts_for_crispresso = defaultdict(int)

    #keep track of distance between R1/R2
    r1_r2_support_distances = defaultdict(int) #for reads that are separated by less than r1_r2_support_max_distance
    r1_r2_not_support_distances = defaultdict(int) #for reads that are separated by more than r1_r2_support_max_distance

    #counts for 'support' status
    r1_r2_support_status_counts = defaultdict(int)

    #column indices in the final file
    ind_offset = cut_point_ind + 1
    r1_id_ind = id_ind
    r1_source_ind = source_ind
    r1_annotation_ind = annotation_ind
    r1_classification_ind = classification_ind
    r1_alignment_ind = alignment_ind
    r1_cut_point_ind = cut_point_ind
    r2_id_ind = id_ind + ind_offset
    r2_source_ind = source_ind + ind_offset
    r2_annotation_ind = annotation_ind + ind_offset
    r2_classification_ind = classification_ind + ind_offset
    r2_alignment_ind = alignment_ind + ind_offset
    r2_cut_point_ind = cut_point_ind + ind_offset
    duplicate_ind = 12
    r1_r2_support_ind = 14
    r1_r2_support_dist_ind = 15
    r1_r2_support_orientation_ind = 16
    duplicate_na_str = "NA\tNA\tNA\tNA" #string describing final cut points for duplicates if paired (because they are duplicates, they are not reported)
    with open(final_file_tmp,'r') as fin, open(final_file,'w') as fout:
        if sorted_r2_file is None:
            #single end
            final_head_els = ['id','source','classification','annotation','alignment','cut_point','read_duplicate_status','read_duplicate_of','final_cut_points','final_cut_keys']
            fout.write("\t".join(final_head_els)+"\n")
            duplicate_ind = 6
            duplicate_na_str = "NA\tNA"
        else:
            #paired end
            final_head_els = ['r1_id','r1_source','r1_classification','r1_annotation','r1_alignment','r1_cut_point','r2_id','r2_source','r2_classification','r2_annotation','r2_alignment','r2_cut_point','read_duplicate_status','read_duplicate_of','r1_r2_support_status','r1_r2_support_distance','r1_r2_support_orientations','r1_final_cut_points','r1_final_cut_keys','r2_final_cut_points','r2_final_cut_keys']
            fout.write("\t".join(final_head_els)+"\n")
        for line in fin:
            line = line.rstrip('\n')
            final_total_count += 1
            line_els = line.split("\t")
            if dedup_based_on_aln_pos and line_els[duplicate_ind] == "duplicate":
                final_duplicate_count += 1
                fout.write("%s\t%s\n"%(line,duplicate_na_str))
                continue

            #now that duplicates are removed, tally support
            if sorted_r2_file is not None:
                r1_r2_support_status = line_els[r1_r2_support_ind]
                r1_r2_support_status_counts[r1_r2_support_status] += 1
                r1_r2_support_dist = line_els[r1_r2_support_dist_ind]
                if r1_r2_support_status == r1_r2_support_str_supported:
                    r1_r2_support_distances[int(r1_r2_support_dist)] += 1
                elif r1_r2_support_status == r1_r2_support_str_not_supported_distance:
                    r1_r2_not_support_distances[int(r1_r2_support_dist)] += 1
                if discard_reads_without_r2_support and r1_r2_support_status != r1_r2_support_str_supported:
                    final_nonsupporting_count_removed += 1
                    fout.write("%s\t%s\n"%(line,duplicate_na_str))
                    continue

            final_processed_count += 1
            r1_classification = line_els[r1_classification_ind]
            r1_cut_points = line_els[r1_cut_point_ind]
            r1_alignment = line_els[r1_alignment_ind]
            r1_source = line_els[r1_source_ind]

            r1_final_source_counts[r1_source] += 1
            r1_final_classification_counts[r1_classification] += 1
            r1_final_source_classification_counts[(r1_source,r1_classification)] += 1

            r1_associated_cut_points = []
            r1_cut_keys = []
            r1_print_str = ""
            r2_print_str = ""

            if r1_classification == 'Linear':
                (r1_chr,r1_range) = r1_alignment.split(":")
                (r1_start,r1_stop) = map(int,r1_range.split("-"))
                #check to see if this read overlaps with any known cutpoints
                if r1_chr in final_cut_points_by_chr:
                    for cut_point in final_cut_points_by_chr[r1_chr]:
                        #if cut point within forward-facing linear read
                        if cut_point > r1_start and cut_point < r1_stop:
                            key = r1_chr+":"+str(cut_point)
                            r1_associated_cut_points.append(key)
                            r1_uncut_counts[key] += 1
                            if r1_source == 'Genome':
                                r1_uncut_genome_counts[key] += 1
                            elif r1_source == 'Fragmented':
                                r1_uncut_frag_counts[key] += 1
                            elif r1_source == 'Custom targets':
                                r1_uncut_custom_aln_counts[key] += 1
                            r1_translocation_cut_counts[key][key] += 1
                            key1 = 'w-' + key
                            key2 = key + '+w'
                            r1_fragment_translocation_cut_counts[key1][key2] += 1
                            r1_fragment_translocation_classifications[key1 + " " + key2] = r1_classification

                            #could have multiple cut point overlaps here
                            r1_cut_key = r1_classification + " " + key
                            r1_read_ids_for_crispresso[line_els[r1_id_ind]].append(r1_cut_key)
                            r1_cut_counts_for_crispresso[r1_cut_key] += 1
                            r1_cut_keys.append(r1_cut_key)
                        #if cut point within reverse-facing linear read
                        elif cut_point > r1_stop and cut_point < r1_start:
                            key = r1_chr+":"+str(cut_point)
                            r1_associated_cut_points.append(key)
                            r1_uncut_counts[key] += 1
                            if r1_source == 'Genome':
                                r1_uncut_genome_counts[key] += 1
                            elif r1_source == 'Fragmented':
                                r1_uncut_frag_counts[key] += 1
                            elif r1_source == 'Custom targets':
                                r1_uncut_custom_aln_counts[key] += 1
                            r1_translocation_cut_counts[key][key] += 1
                            key1 = 'w+' + key
                            key2 = key + '-w' 
                            r1_fragment_translocation_cut_counts[key1][key2] += 1
                            r1_fragment_translocation_classifications[key1 + " " + key2] = r1_classification

                            #could have multiple cut point overlaps here
                            r1_cut_key = r1_classification + " " + key
                            r1_read_ids_for_crispresso[line_els[r1_id_ind]].append(r1_cut_key)
                            r1_cut_counts_for_crispresso[r1_cut_key] += 1
                            r1_cut_keys.append(r1_cut_key)

            else:
                #if it's not Linear, we've already annotated the cut points
                #if we only have one cut point, it's a chimeric cut
                cut_1 = r1_cut_points.split("~")[0]
                cut_2 = cut_1
                if '~' in r1_cut_points:
                    cut_1,cut_2 = r1_cut_points.split("~")
                if "*" not in cut_1:
                    cut_1 = final_cut_point_lookup[cut_1]
                    r1_associated_cut_points.append(cut_1)
                    r1_cut_counts[cut_1] += 1
                    if r1_source == 'Fragmented':
                        r1_cut_frag_counts[cut_1] += 1
                    elif r1_source == 'Custom targets':
                        r1_cut_custom_aln_counts[cut_1] += 1
                else:
                    cut_1 = "*"

                if "*" not in cut_2:
                    cut_2 = final_cut_point_lookup[cut_2]
                    r1_associated_cut_points.append(cut_2)
                    r1_cut_counts[cut_2] += 1
                    if r1_source == 'Fragmented':
                        r1_cut_frag_counts[cut_2] += 1
                    elif r1_source == 'Custom targets':
                        r1_cut_custom_aln_counts[cut_2] += 1
                else:
                    cut_2 = "*"


                r1_translocation_cut_counts[cut_1][cut_2] += 1

                #also keep track of read ids covering each cut for analysis by CRISPResso
                #this_anno has cut stradn information as well as artificial target name
                this_anno = line_els[r1_annotation_ind].split(" ")[0]

                #add strand and direction information from anno column to newly-identified final cuts
                frag_key_1 = this_anno[0:2]+cut_1
                frag_key_2 = cut_2+this_anno[-2:]

                r1_cut_key = "%s %s~%s"%(r1_classification,frag_key_1,frag_key_2)
                r1_read_ids_for_crispresso[line_els[r1_id_ind]].append(r1_cut_key)
                r1_cut_counts_for_crispresso[r1_cut_key] += 1
                r1_cut_keys.append(r1_cut_key)


                r1_fragment_translocation_cut_counts[frag_key_1][frag_key_2] += 1
                r1_fragment_translocation_classifications[frag_key_1 + " " + frag_key_2] = r1_classification

            #make window alignment assignments
            for r1_aln_pos in r1_alignment.split("~"):
                if '*' not in r1_aln_pos:
                    (aln_chr,aln_range) = r1_aln_pos.split(":")
                    (aln_pos,aln_end) = aln_range.split("-")
                    aln_pos_window = int(aln_pos)/genome_map_resolution

                    if aln_chr not in r1_aln_pos_counts_by_chr:
                        r1_aln_pos_counts_by_chr[aln_chr] = defaultdict(int)
                        r1_aln_pos_genome_counts_by_chr[aln_chr] = defaultdict(int)
                        r1_aln_pos_custom_aln_counts_by_chr[aln_chr] = defaultdict(int)
                        r1_aln_pos_frag_counts_by_chr[aln_chr] = defaultdict(int)
                    #assign each read to neighboring windows
                    r1_aln_pos_counts_by_chr[aln_chr][aln_pos_window] += 1
                    if r1_source == 'Genome':
                        r1_aln_pos_genome_counts_by_chr[aln_chr][aln_pos_window] += 1
                    elif r1_source == 'Fragmented':
                        r1_aln_pos_frag_counts_by_chr[aln_chr][aln_pos_window] += 1
                    elif r1_source == 'Custom targets':
                        r1_aln_pos_custom_aln_counts_by_chr[aln_chr][aln_pos_window] += 1


            r1_cut_points_str = ",".join(r1_associated_cut_points) if r1_associated_cut_points else "NA"
            r1_cut_keys_str = ",".join(r1_cut_keys) if r1_cut_keys else "NA"
            r1_print_str = r1_cut_points_str+"\t"+r1_cut_keys_str

            #now do R2
            if sorted_r2_file is not None:
                r2_classification = line_els[r2_classification_ind]
                r2_cut_points = line_els[r2_cut_point_ind]
                r2_alignment = line_els[r2_alignment_ind]
                r2_source = line_els[r2_source_ind]
                #alignment for window is only chr:pos (Linear reads have a chr:pos-end format, so we fix that with this variable)

                r2_final_source_counts[r2_source] += 1
                r2_final_classification_counts[r2_classification] += 1
                r2_final_source_classification_counts[(r2_source,r2_classification)] += 1

                r2_associated_cut_points = []
                r2_cut_keys = []

                if r2_classification == 'Linear':
                    (r2_chr,r2_range) = r2_alignment.split(":")
                    (r2_start,r2_stop) = r2_range.split("-")
                    #check to see if this read overlaps with any known cutpoints
                    if r2_chr in final_cut_points_by_chr:
                        for cut_point in final_cut_points_by_chr[r2_chr]:
                            if cut_point > int(r2_start) and cut_point < int(r2_stop):
                                key = r2_chr+":"+str(cut_point)
                                r2_associated_cut_points.append(key)
                                r2_uncut_counts[key] += 1
                                if r2_source == 'Genome':
                                    r2_uncut_genome_counts[key] += 1
                                elif r2_source == 'Fragmented':
                                    r2_uncut_frag_counts[key] += 1
                                elif r2_source == 'Custom targets':
                                    r2_uncut_custom_aln_counts[key] += 1
                                r2_translocation_cut_counts[key][key] += 1
                                key1 = 'w-' + key
                                key2 = key + '+w'
                                r2_fragment_translocation_cut_counts[key1][key2] += 1

                                r2_cut_key = r2_classification + " " + key
                                r2_read_ids_for_crispresso[line_els[r2_id_ind]].append(r2_cut_key)
                                r2_cut_counts_for_crispresso[r2_cut_key] += 1
                                r2_cut_keys.append(r2_cut_key)
                else:
                    #if it's not Linear, we've already annotated the cut points
                    #if we only have one cut point, it's a chimeric cut
                    cut_1 = r2_cut_points.split("~")[0]
                    cut_2 = cut_1
                    if '~' in r2_cut_points:
                        cut_1,cut_2 = r2_cut_points.split("~")

                    if "*" not in cut_1:
                        cut_1 = final_cut_point_lookup[cut_1]
                        r2_associated_cut_points.append(cut_1)
                        r2_cut_counts[cut_1] += 1
                        if r2_source == 'Fragmented':
                            r2_cut_frag_counts[cut_1] += 1
                        elif r2_source == 'Custom targets':
                            r2_cut_custom_aln_counts[cut_1] += 1
                    else:
                        cut_1 = "*"

                    if "*" not in cut_2:
                        cut_2 = final_cut_point_lookup[cut_2]
                        r2_associated_cut_points.append(cut_2)
                        r2_cut_counts[cut_2] += 1
                        if r2_source == 'Fragmented':
                            r2_cut_frag_counts[cut_2] += 1
                        elif r2_source == 'Custom targets':
                            r2_cut_custom_aln_counts[cut_2] += 1
                    else:
                        cut_2 = "*"


                    r2_translocation_cut_counts[cut_1][cut_2] += 1

                    #also keep track of read ids covering each cut for analysis by CRISPResso
                    this_anno = line_els[r2_annotation_ind].split(" ")[0]
                    #add strand and direction information from anno column to newly-identified funal cuts
                    frag_key_1 = this_anno[0:2]+cut_1
                    frag_key_2 = cut_2+this_anno[-2:]

                    r2_cut_key = "%s %s~%s"%(r2_classification,frag_key_1,frag_key_2)
                    r2_read_ids_for_crispresso[line_els[r2_id_ind]].append(r2_cut_key)
                    r2_cut_counts_for_crispresso[r2_cut_key] += 1
                    r2_cut_keys.append(r2_cut_key)

                    r2_fragment_translocation_cut_counts[frag_key_1][frag_key_2] += 1


                #make window alignment assignments
                for r2_aln_pos in r2_alignment.split("~"):
                    if '*' not in r2_aln_pos:
                        (aln_chr,aln_range) = r2_aln_pos.split(":")
                        aln_pos = aln_range.split("-")[0]
                        aln_pos_window = int(aln_pos)/genome_map_resolution
                        if aln_chr not in r2_aln_pos_counts_by_chr:
                            r2_aln_pos_counts_by_chr[aln_chr] = defaultdict(int)
                            r2_aln_pos_genome_counts_by_chr[aln_chr] = defaultdict(int)
                            r2_aln_pos_custom_aln_counts_by_chr[aln_chr] = defaultdict(int)
                            r2_aln_pos_frag_counts_by_chr[aln_chr] = defaultdict(int)
                        #assign each read to neighboring windows
                        r2_aln_pos_counts_by_chr[aln_chr][aln_pos_window] += 1
                        if r2_source == 'Genome':
                            r2_aln_pos_genome_counts_by_chr[aln_chr][aln_pos_window] += 1
                        elif r2_source == 'Fragmented':
                            r2_aln_pos_frag_counts_by_chr[aln_chr][aln_pos_window] += 1
                        elif r2_source == 'Custom targets':
                            r2_aln_pos_custom_aln_counts_by_chr[aln_chr][aln_pos_window] += 1
                r2_cut_points_str = ",".join(r2_associated_cut_points) if r2_associated_cut_points else "NA"
                r2_cut_keys_str = ",".join(r2_cut_keys) if r2_cut_keys else "NA"
                r2_print_str = "\t"+r2_cut_points_str+"\t"+r2_cut_keys_str


            #finished r2



            fout.write("%s\t%s%s\n"%(line,r1_print_str,r2_print_str))

    logging.info('Processed %d reads. Of these,\n    %d duplicate reads were removed\n    %d reads were discarded because they were not supported by R2\n  %d reads were used for final calculations'%(final_total_count,final_duplicate_count, final_nonsupporting_count_removed, final_processed_count))


    final_observed_cuts_count = 0
    r1_cut_report = root + ".r1_cut_report.txt"
    r2_cut_report = root + ".r2_cut_report.txt"
    if sorted_r2_file is None:
        #single-end cut report
        r2_cut_report = None
        with open (r1_cut_report,'w') as r1_cut_out:
            r1_cut_out.write('chr\tstart\tcut_tot\tuncut_tot\tuncut_genome\tcut_custom\tuncut_custom\tcut_frags\tuncut_frags\tcut_provided_as_input\n')
            for chrom in sorted(final_cut_points_by_chr.keys()):
                if chrom == "*":
                    continue
                for pos in sorted(final_cut_points_by_chr[chrom]):
                    final_observed_cuts_count += 1
                    cut_key = chrom + ":" + str(pos)
                    cut_annotation = 'Novel'
                    if cut_key in cut_annotations:
                        cut_annotation = cut_annotations[cut_key]
                    r1_cut_out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"%(chrom,pos,
                        r1_cut_counts[cut_key],
                        r1_uncut_counts[cut_key],
                        r1_uncut_genome_counts[cut_key],
                        r1_cut_custom_aln_counts[cut_key],
                        r1_uncut_custom_aln_counts[cut_key],
                        r1_cut_frag_counts[cut_key],
                        r1_uncut_frag_counts[cut_key],
                        cut_annotation))

        logging.info('Wrote cut report ' + r1_cut_report + ' for ' + str(final_observed_cuts_count) + ' cuts')
    else:
        #paired cut report
        with open (r1_cut_report,'w') as r1_cut_out, open(r2_cut_report,'w') as r2_cut_out:
            r1_cut_out.write('chr\tstart\tcut_tot\tuncut_tot\tuncut_genome\tcut_custom\tuncut_custom\tcut_frags\tuncut_frags\tcut_provided_as_input\n')
            r2_cut_out.write('chr\tstart\tcut_tot\tuncut_tot\tuncut_genome\tcut_custom\tuncut_custom\tcut_frags\tuncut_frags\tcut_provided_as_input\n')
            for chrom in sorted(final_cut_points_by_chr.keys()):
                if chrom == "*":
                    continue
                for pos in sorted(final_cut_points_by_chr[chrom]):
                    final_observed_cuts_count += 1
                    cut_key = chrom + ":" + str(pos)
                    cut_annotation = 'Novel'
                    if cut_key in cut_annotations:
                        cut_annotation = cut_annotations[cut_key]
                    r1_cut_out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"%(chrom,pos,
                        r1_cut_counts[cut_key],
                        r1_uncut_counts[cut_key],
                        r1_uncut_genome_counts[cut_key],
                        r1_cut_custom_aln_counts[cut_key],
                        r1_uncut_custom_aln_counts[cut_key],
                        r1_cut_frag_counts[cut_key],
                        r1_uncut_frag_counts[cut_key],
                        cut_annotation))
                    r2_cut_out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"%(chrom,pos,
                        r2_cut_counts[cut_key],
                        r2_uncut_counts[cut_key],
                        r2_uncut_genome_counts[cut_key],
                        r2_cut_custom_aln_counts[cut_key],
                        r2_uncut_custom_aln_counts[cut_key],
                        r2_cut_frag_counts[cut_key],
                        r2_uncut_frag_counts[cut_key],
                        cut_annotation))

        logging.info('Wrote cut reports ' + r1_cut_report + ' ' + r2_cut_report + ' for ' + str(final_observed_cuts_count) + ' cuts')

    def write_translocation_reports(read_num,translocation_cut_counts,fragment_translocation_cut_counts):
        tx_keys = sorted(translocation_cut_counts.keys())
        tx_list_report = root + "." + read_num + "_translocation_list.txt"
        with open (tx_list_report,'w') as fout:
            fout.write('cut_1\tcut_1_anno\tcut_2\tcut_2_anno\tcount\n')
            for tx_key_1 in tx_keys:
                cut_1_annotation = 'Novel'
                if tx_key_1 in cut_annotations:
                    cut_1_annotation = cut_annotations[tx_key_1]
                for tx_key_2 in sorted(translocation_cut_counts[tx_key_1].keys()):
                    cut_2_annotation = 'Novel'
                    if tx_key_2 in cut_annotations:
                        cut_2_annotation = cut_annotations[tx_key_2]
                    #only write if counts are > 0
                    if translocation_cut_counts[tx_key_1][tx_key_2] > 0:
                        fout.write("%s\t%s\t%s\t%s\t%s\n"%(tx_key_1,cut_1_annotation,tx_key_2,cut_2_annotation,translocation_cut_counts[tx_key_1][tx_key_2]))
        logging.info('Wrote ' + read_num + ' translocation list ' + tx_list_report)

        tx_report = root + "." + read_num + "_translocation_table.txt"
        with open (tx_report,'w') as fout:
            fout.write('Translocations\t'+"\t".join(tx_keys)+"\n")
            for tx_key in tx_keys:
                fout.write(tx_key+"\t"+"\t".join([str(translocation_cut_counts[tx_key][x]) for x in tx_keys])+"\n")
        logging.info('Wrote ' + read_num + ' translocation table ' + tx_report)

        tx_ontarget_report = root + '.' + read_num + "_ontarget_translocation_table.txt"
        on_targets = sorted([x for x in cut_annotations if cut_annotations[x] == 'On-target'])
        if len(on_targets) > 0:
            on_to_off_target_set = set()
            for on_target in on_targets:
                if on_target in translocation_cut_counts:
                    for tx_key_2 in translocation_cut_counts[on_target].keys():
                        on_to_off_target_set.add(tx_key_2)
            #put on-targets first
            for on_target in on_targets:
                on_to_off_target_set.discard(on_target)
            on_to_off_target_list = on_targets + sorted(list(on_to_off_target_set))

            with open (tx_ontarget_report,'w') as fout:
                fout.write('Translocations\t'+"\t".join(on_to_off_target_list)+"\n")
                for tx_key in on_targets:
                    fout.write(tx_key+"\t"+"\t".join([str(translocation_cut_counts[tx_key][x]) if tx_key in translocation_cut_counts else str(0) for x in on_to_off_target_list])+"\n")
            logging.info('Wrote ' + read_num + ' on-target translocation table ' + tx_ontarget_report)

        # fragment translocations
        tx1_keys = sorted(fragment_translocation_cut_counts.keys())
        tx2_keys = defaultdict(int)
        tx_list_report = root + '.' + read_num + "_fragment_translocation_list.txt"
        with open (tx_list_report,'w') as fout:
            fout.write('fragment_1\tfragment_1_anno\tfragment_2\tfragment_2_anno\tcount\tannotation\n')
            for tx_key_1 in tx1_keys:
                cut_1 = tx_key_1[2:]
                cut_1_annotation = 'Novel'
                if cut_1 in cut_annotations:
                    cut_1_annotation = cut_annotations[cut_1]
                for tx_key_2 in sorted(fragment_translocation_cut_counts[tx_key_1].keys()):
                    cut_2 = tx_key_2[:-2]
                    cut_2_annotation = 'Novel'
                    if cut_2 in cut_annotations:
                        cut_2_annotation = cut_annotations[cut_2]
                    tx_classification = 'Unknown'
                    tx_class_key = tx_key_1 + " " + tx_key_2
                    if tx_class_key in r1_fragment_translocation_classifications:
                        tx_classification = r1_fragment_translocation_classifications[tx_class_key]
                    # only write if counts exist
                    if fragment_translocation_cut_counts[tx_key_1][tx_key_2] > 0:
                        fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(tx_key_1,cut_1_annotation,tx_key_2,cut_2_annotation,fragment_translocation_cut_counts[tx_key_1][tx_key_2],tx_classification))
                    tx2_keys[tx_key_2] += 1
        logging.info('Wrote ' + read_num + ' fragment translocation list ' + tx_list_report)

        tx_report = root + "." + read_num + "_fragment_translocation_table.txt"
        tx2_keys_to_report = sorted(tx2_keys)
        with open (tx_report,'w') as fout:
            fout.write('Translocations\t'+"\t".join(tx2_keys_to_report)+"\n")
            for tx_key in tx1_keys:
                fout.write(tx_key+"\t"+"\t".join([str(fragment_translocation_cut_counts[tx_key][x]) if tx_key in fragment_translocation_cut_counts else str(0) for x in tx2_keys_to_report])+"\n")
        logging.info('Wrote ' + read_num + ' fragment translocation table ' + tx_report)

        tx_ontarget_report = root + "." + read_num + "_ontarget_fragment_translocation_table.txt"
        on_targets = sorted([x for x in cut_annotations if cut_annotations[x] == 'On-target'])
        on_targets_1_to_report = defaultdict(int)
        on_targets_2_to_report = defaultdict(int)
        if len(on_targets) > 0:
            on_to_off_target_set = set()
            for on_target in on_targets:
                for prefix in ['w-','w+','c-','c+']:
                    tx_key_1 = prefix+on_target
                    if tx_key_1 in fragment_translocation_cut_counts:
                        on_targets_1_to_report[tx_key_1] += 1
                        for tx_key_2 in fragment_translocation_cut_counts[tx_key_1].keys():
                            on_targets_2_to_report[tx_key_2] += 1

            on_targets_1_list = sorted(on_targets_1_to_report)
            on_targets_2_list = sorted(on_targets_2_to_report)
            with open (tx_ontarget_report,'w') as fout:
                fout.write('Translocations\t'+"\t".join(on_targets_2_list)+"\n")
                for tx_key in on_targets_1_to_report:
                    fout.write(tx_key+"\t"+"\t".join([str(fragment_translocation_cut_counts[tx_key][x]) if tx_key in fragment_translocation_cut_counts else str(0) for x in on_targets_2_list])+"\n")
            logging.info('Wrote ' + read_num + ' on-target fragment translocation table ' + tx_ontarget_report)

    #r1 translocations
    write_translocation_reports('r1',r1_translocation_cut_counts,r1_fragment_translocation_cut_counts)
    #r2 translocations
    if sorted_r2_file is not None:
        write_translocation_reports('r2',r2_translocation_cut_counts,r2_fragment_translocation_cut_counts)

    #deduplication plot
    deduplication_plot_obj = None
    if total_reads_processed > 0:
        labels = ['Not duplicate','Duplicate']
        values = [total_reads_processed-dups_seen,dups_seen]
        deduplication_plot_obj_root = root + ".deduplication"
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Duplicate read counts')
        plt.savefig(deduplication_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(deduplication_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        dedup_note = ""
        if not dedup_based_on_aln_pos:
            dedup_note = ". Note that deduplication of read for final counts is off by default but can be turned on by running with the parameter --dedup_input_based_on_aln_pos_and_UMI."

        deduplication_plot_obj = PlotObject(
                plot_name = deduplication_plot_obj_root,
                plot_title = 'Duplicate read counts',
                plot_label = 'Number of reads that were counted as deduplication based on alignment and UMI' + dedup_note,
                plot_datas = [('Read assignments',final_file)]
                )

    def write_aln_reports_and_make_plots(read_num,aln_pos_counts_by_chr,aln_pos_genome_counts_by_chr,aln_pos_custom_aln_counts_by_chr,aln_pos_frag_counts_by_chr,final_file,final_source_counts,final_classification_counts,final_source_classification_counts):
        """
        Create reports and plots for alignment, source, and classification of reads

        params:
            read_num: either 'r1' or 'r2' depending on reads to process
            aln_pos_counts_by_chr: dict of chr->count of reads
            aln_pos_genome_counts_by_chr: dict of chr->count of reads aligned to genome
            aln_pos_custom_aln_counts_by_chr: dict of chr->count of reads aligned via alignment to custom targets
            aln_pos_frag_counts_by_chr: dict of chr->count of reads aligned via fragmentation
            final_file: final read assignment filename
            final_source_counts: dict of source->count where source is Genome, Custom targets, or Fragmented
            final_classification_counts: dict of classification->count where classification is Primed, Linear, Chimera, Translocation, Unidentified
            final_source_classification_counts: dict of [classification][source]->count

        returns:
            classification_plot_obj: plot_obj showing classification of reads
            source_plot_obj: plot_obj showing source of reads
        """
        read_num_capitalized = read_num.capitalize()

        #alignment report
        align_report = root + "." + read_num + "_align_report.txt"
        with open (align_report,'w') as fout:
            fout.write('chr\tstart\tread_tot\tread_genome\tread_custom\tread_frags\n')
            for chrom in sorted(aln_pos_counts_by_chr.keys()):
                for pos in sorted(aln_pos_counts_by_chr[chrom].keys()):
                    fout.write("%s\t%d\t%d\t%d\t%d\t%d\n"%(chrom,pos*genome_map_resolution,
                        aln_pos_counts_by_chr[chrom][pos],
                        aln_pos_genome_counts_by_chr[chrom][pos],
                        aln_pos_custom_aln_counts_by_chr[chrom][pos],
                        aln_pos_frag_counts_by_chr[chrom][pos]))

        source_labels = ['Genome','Custom targets','Fragmented','Unmapped']
        classification_labels = ['Primed','Linear','Chimera','Large deletion','Large inversion','Translocation','Unidentified']

        #check to make sure labels match
        for k in final_source_counts.keys():
            if k not in source_labels:
                raise Exception('Source label %s not found in %s'%(k,str(source_labels)))

        source_classification_file = root + "." + read_num + "_source_classification.txt"
        with open(source_classification_file,'w') as summary:
            summary.write("\t"+"\t".join(classification_labels)+"\n")
            for source_label in source_labels:
                summary.write(source_label + "\t" + "\t".join([str(final_source_classification_counts[(source_label,classification_label)]) for classification_label in classification_labels]) + "\n")


        #sources plot
        values = [final_source_counts[x] for x in source_labels]
        source_plot_obj_root = root + "." + read_num + "_sources"
        with open(source_plot_obj_root+".txt",'w') as summary:
            summary.write("\t".join(source_labels)+"\n")
            summary.write("\t".join([str(x) for x in values])+"\n")
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(source_labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(source_labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Read alignment')
        plt.savefig(source_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(source_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(source_labels,values)])
        source_plot_obj = PlotObject(
                plot_name = source_plot_obj_root,
                plot_title = read_num_capitalized + ' Read Alignments',
                plot_label = 'Method used to determine the status of each ' + read_num_capitalized + ' read<br>'+plot_count_str,
                plot_datas = [
                    (read_num_capitalized + ' alignment sources',source_plot_obj_root + ".txt"),
                    (read_num_capitalized + ' alignment sources by classification',source_classification_file),
                    ('Read assignments',final_file)
                    ]
                )

        #assignment plot
        values = [final_classification_counts[x] for x in classification_labels]
        classification_plot_obj_root = root + "." + read_num + "_classifications"
        with open(classification_plot_obj_root+".txt",'w') as summary:
            summary.write("\t".join(classification_labels)+"\n")
            summary.write("\t".join([str(x) for x in values])+"\n")
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        pie_values = []
        pie_labels = []
        for i in range(len(classification_labels)):
            if values[i] > 0:
                pie_values.append(values[i])
                pie_labels.append(classification_labels[i]+"\n("+str(values[i])+")")
        ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
        ax.set_title('Read alignment')
        plt.savefig(classification_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(classification_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')

        plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(classification_labels,values)])
        classification_plot_obj = PlotObject(
                plot_name = classification_plot_obj_root,
                plot_title = read_num_capitalized + ' Read Classification',
                plot_label = read_num_capitalized + ' Read classification<br>'+plot_count_str,
                plot_datas = [
                    (read_num_capitalized + ' alignment classifications',classification_plot_obj_root + ".txt"),
                    (read_num_capitalized + ' alignment sources by classification',source_classification_file),
                    ('Read assignments',final_file)
                    ]
                )
        return(source_plot_obj,classification_plot_obj)

    r1_source_plot_obj,r1_classification_plot_obj = write_aln_reports_and_make_plots('r1',r1_aln_pos_counts_by_chr,r1_aln_pos_genome_counts_by_chr,r1_aln_pos_custom_aln_counts_by_chr,r1_aln_pos_frag_counts_by_chr,final_file,r1_final_source_counts,r1_final_classification_counts,r1_final_source_classification_counts)

    r2_source_plot_obj = None #source of alignment for each read (genome/custom/frags)
    r2_classification_plot_obj = None #classification of each read
    #r2 alignment report
    if sorted_r2_file is not None:
        r2_source_plot_obj,r2_classification_plot_obj = write_aln_reports_and_make_plots('r2',r2_aln_pos_counts_by_chr,r2_aln_pos_genome_counts_by_chr,r2_aln_pos_custom_aln_counts_by_chr,r2_aln_pos_frag_counts_by_chr,final_file,r2_final_source_counts,r2_final_classification_counts,r2_final_source_classification_counts)

    r1_r2_support_plot_obj = None # plot of how many reads for which r1 and r2 support each other
    r1_r2_support_dist_plot_obj = None # plot of how far apart the r1/r2 supporting reads were aligned from each other
    r1_r2_no_support_dist_plot_obj = None # plot of how far apart the r1/r2 reads that did NOT support each other (but were on the same chromosome)
    if sorted_r2_file is not None:
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
            for key in sorted(r1_r2_support_distances.keys()):
                summary.write(str(key) + '\t' + str(r1_r2_support_distances[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(range(len(keys)),vals,tick_label=keys)
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
            for key in sorted(r1_r2_not_support_distances.keys()):
                summary.write(str(key) + '\t' + str(r1_r2_not_support_distances[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(range(len(keys)),vals,tick_label=keys)
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
        fout.write("\t".join(['total_reads_processed','duplicate_reads_discarded','discarded_reads_not_supported_by_R2','final_read_count','unique_cuts_observed_count'])+"\n")
        fout.write("\t".join(str(x) for x in [final_total_count,final_duplicate_count,final_nonsupporting_count_removed, final_processed_count, final_observed_cuts_count])+"\n")

        deduplication_plot_obj_str = deduplication_plot_obj.to_json()
        fout.write(deduplication_plot_obj_str+"\n")

        r1_source_plot_obj_str = r1_source_plot_obj.to_json()
        fout.write(r1_source_plot_obj_str+"\n")

        r1_classification_plot_obj_str = r1_classification_plot_obj.to_json()
        fout.write(r1_classification_plot_obj_str+"\n")

        r2_source_plot_obj_str = "None"
        if r2_source_plot_obj is not None:
            r2_source_plot_obj_str = r2_source_plot_obj.to_json()
        fout.write(r2_source_plot_obj_str+"\n")

        r2_classification_plot_obj_str = "None"
        if r2_classification_plot_obj is not None:
            r2_classification_plot_obj_str = r2_classification_plot_obj.to_json()
        fout.write(r2_classification_plot_obj_str+"\n")

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


    return final_file,r1_read_ids_for_crispresso,r1_cut_counts_for_crispresso,r2_read_ids_for_crispresso,r2_cut_counts_for_crispresso,deduplication_plot_obj,r1_source_plot_obj,r1_classification_plot_obj,r2_source_plot_obj,r2_classification_plot_obj,r1_r2_support_plot_obj,r1_r2_support_dist_plot_obj,r1_r2_no_support_dist_plot_obj


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

def chop_reads(root,unmapped_reads_fastq_r1,unmapped_reads_fastq_r2,bowtie2_genome,fragment_size,fragment_step_size,min_seen_frag_cutoff=2,fragment_min_mapq=30,bowtie2_command='bowtie2',bowtie2_threads=1,samtools_command='samtools',can_use_previous_analysis=False):
    """
    Creates chopped reads
    Keeps track of locations where reads map
    For reads that do not map, creates chopped fragments which are to be mapped

    params:
        root: root for written files
        unmapped_reads_fastq_r1: fastq file of r1 unmappable reads to chop
        unmapped_reads_fastq_r2: fastq file of r2 unmappable reads to chop
        bowtie2_genome: location of bowtie2 indexed reference genome
        fragment_size: size of resulting chopped fragments.
        fragment_step_size: stride between chopped fragments. 5 results in a new fragment every 5bp
        min_seen_frag_cutoff: minimum number of fragments aligned to a location for the location to be considered as a possible alignment location for a read (e.g. if a read only has 3 fragments aligning to location A and 1 fragment aligning to location B, a min_seen_frag_cutoff of 2 would not allow position B
        fragment_min_mapq: min mapq (mapping quality) for a fragment to be considered
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        samtools_command: location of samtools to run
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch


    returns:
        r1_assignments_file: file of location assignments for each r1 read
        r2_assignments_file: file of location assignments for each r2 read (None if single-end)
        linear_count: number of reads that were identified as linear (for some reason they didn't map using global mapping, but they map to a linear piece of DNA - eg there are bits on the end of the read that don't map)
        translocation_count: number of reads that were identified as translocations
        large_deletion_count: number of reads that were identified as large deletion (same chr, different orientation(e.g. one extending left, one extending right))
        large_inversion_count: number of reads that were identified as large inversions (same chr, same orientation (e.g. both reading left))
        unidentified_count: number of reads that couldn't be identified as translocations
        frags_plot_obj: plot object summarizing fragments
    """
    info_file = root + '.info'
    if os.path.isfile(info_file) and can_use_previous_analysis:
        frags_mapped_count = -1
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 8:
                (r1_assignments_file_str,r2_assignments_file_str,linear_count_str,translocation_count_str,large_deletion_count_str,large_inversion_count_str,unidentified_count_str,frags_mapped_count_str) = line_els
                r1_assignments_file = None if r1_assignments_file_str == "None" else r1_assignments_file_str
                r2_assignments_file = None if r2_assignments_file_str == "None" else r2_assignments_file_str
                linear_count = int(linear_count_str)
                translocation_count  = int(translocation_count_str)
                large_deletion_count  = int(large_deletion_count_str)
                large_inversion_count  = int(large_inversion_count_str)
                unidentified_count  = int(unidentified_count_str)
                frags_mapped_count  = int(frags_mapped_count_str)

                frags_plot_obj_str = fin.readline().rstrip('\n')
                frags_plot_obj = None
                if frags_plot_obj_str != "" and frags_plot_obj_str != "None":
                    frags_plot_obj = PlotObject.from_json(frags_plot_obj_str)

                if frags_mapped_count >= 0:
                    logging.info('Using previously-processed fragment analysis for %d mapped fragments (%d linear reads, %d translocations, %d large deletions, %d large inversions, and %d unidentified reads)'%(frags_mapped_count,linear_count,translocation_count,large_deletion_count,large_inversion_count,unidentified_count))
                    return (r1_assignments_file,r2_assignments_file,linear_count,translocation_count,large_deletion_count,large_inversion_count,unidentified_count,frags_plot_obj)
                else:
                    logging.info('Could not recover previously-analyzed fragments. Reanalyzing.')

    logging.info('Creating read fragments')

    unmapped_ids = {}
    unmapped_id = 0
    max_frags = 0 #max frags produced for a read
    frags_per_read = {} #number of fragments created per read

    read_files_to_frag = []
    if unmapped_reads_fastq_r1 is not None:
        read_files_to_frag.append(unmapped_reads_fastq_r1)
    if unmapped_reads_fastq_r2 is not None:
        read_files_to_frag.append(unmapped_reads_fastq_r2)
    unmapped_frag_file = root + ".to_map.fq"
    with open(unmapped_frag_file,"w") as unmapped_fastq:
        for (read_idx,unmapped_reads_fastq) in enumerate(read_files_to_frag):
            with open(unmapped_reads_fastq,'r') as unmapped_reads:
                line_id = unmapped_reads.readline()
                while(line_id):
                    line_id = line_id.rstrip('\n')
                    line_seq = unmapped_reads.readline().rstrip('\n')
                    line_plus = unmapped_reads.readline().rstrip('\n')
                    line_qual = unmapped_reads.readline().rstrip('\n')

                    this_id = line_id

                    frag_num = 0
                    offset = 0
                    while offset + fragment_size < len(line_seq):
                        new_id = "%s.CLR%d.CLID%s.CLO%s.CLF%s"%(this_id,read_idx+1,unmapped_id,offset,frag_num)
                        new_seq =  line_seq[offset:offset + fragment_size]
                        new_qual = line_qual[offset:offset + fragment_size]
                        unmapped_fastq.write("%s\n%s\n+\n%s\n"%(new_id,new_seq,new_qual))
                        frag_num += 1
                        offset += fragment_step_size

                    #write last frag
                    if len(line_seq) > fragment_size:
                        offset = len(line_seq)-fragment_size
                        new_id = "%s.CLR%d.CLID%s.CLO%s.CLF%s"%(this_id,read_idx+1,unmapped_id,offset,frag_num)
                        new_seq =  line_seq[offset:offset + fragment_size]
                        new_qual = line_qual[offset:offset + fragment_size]
                        unmapped_fastq.write("%s\n%s\n+\n%s\n"%(new_id,new_seq,new_qual))
                        frag_num += 1

                    #write frags that are too short
                    if len(line_seq) <= fragment_size:
                        new_id = "%s.CLR%d.CLID%s.CLO%s.CLF%s"%(this_id,read_idx+1,unmapped_id,0,0)
                        unmapped_fastq.write("%s\n%s\n+\n%s\n"%(new_id,line_seq,line_qual))
                        frag_num += 1

                    if frag_num > max_frags:
                        max_frags = frag_num
                    if frag_num not in frags_per_read:
                        frags_per_read[frag_num] = 0
                    frags_per_read[frag_num] += 1

                    unmapped_id += 1

                    line_id = unmapped_reads.readline()

    number_unmapped_reads_chopped = unmapped_id
    logging.info('Created chopped reads for ' + str(number_unmapped_reads_chopped) + ' reads')

    frags_per_unaligned_read_root = root + ".fragsPerUnalignedRead"
    keys = sorted(frags_per_read.keys())
    if len(keys) > 0:
        vals = [frags_per_read[key] for key in keys]
        with open(frags_per_unaligned_read_root + ".txt","w") as frags:
            frags.write('numFragments\tnumReads\n')
            for key in keys:
                frags.write(str(key) + '\t' + str(frags_per_read[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        ax.bar(range(len(keys)),vals,tick_label=keys)
        ax.set_ymargin(0.05)
        ax.set_xlabel('Number of fragments')
        ax.set_ylabel('Number of Reads')
        ax.set_title('Number of fragments per unaligned read')

        plt.savefig(frags_per_unaligned_read_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(frags_per_unaligned_read_root+".png",pad_inches=1,bbox_inches='tight')

        frags_plot_obj = PlotObject(
                plot_name = frags_per_unaligned_read_root,
                plot_title = 'Fragments per unaligned read',
                plot_label = 'Bar plot showing number of fragments produced per unaligned read',
                plot_datas = [('Fragments per unaligned read',frags_per_unaligned_read_root + ".txt")]
                )
    else:
        frags_plot_obj = None


    mapped_chopped_sam_file = root + ".mapped.sam"
    chopped_bowtie_log = root + ".mapped.bowtie2Log"
    chopped_aln_command = '%s --seed 2248 --sam-no-qname-trunc --reorder --no-hd --threads %d --end-to-end -x %s -U %s -S %s' %(bowtie2_command,bowtie2_threads,bowtie2_genome,unmapped_frag_file,mapped_chopped_sam_file)

    logging.info('Aligning chopped reads')

    logging.debug(chopped_aln_command)
    aln_result = subprocess.check_output(chopped_aln_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
    logging.debug('Alignment of chopped fragments to genome: ' + aln_result)
    with open (chopped_bowtie_log,'w') as lout:
        lout.write('Alignment of chopped fragments to genome\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(chopped_aln_command,aln_result))

    logging.info('Analyzing chopped reads')

    def analyze_curr_id_vals(curr_id_vals,curr_id_aln_pos,fragment_size,fragment_step_size,max_frags,min_seen_frag_cutoff=2):
        """
        Helper function to analyze a list of curr_id_vals
        A midpoint is found that splits the aligned fragments such that the number of locations in each split is minimized.
        A string describing this read and the aligned fragments is also created
        I pulled it out into a separate function because it's used in two spots and may need to be changed

        params:
            curr_id_vals: dict of fragNum:alnPos for each fragment in a specific read - this is the start of the read that would have produced this fragment
            curr_id_aln_pos: dict of fragNum:alnPos for each fragment in a specific read - this is the actual alignment position of the fragment
            fragment_size: size of resulting chopped fragments.
            fragment_step_size: stride between chopped fragments. 5 results in a new fragment every 5bp
            max_frags: the max number of frags to print (for constructing the outline to be printed)
            min_seen_frag_cutoff: minimum number of fragments aligned to a location for the location to be considered as a possible alignment location for a read (e.g. if a read only has 3 fragments aligning to location A and 1 fragment aligning to location B, a min_seen_frag_cutoff of 2 would not allow position B
        returns:
            outline: String to be printed to output with locationa and max_frags number of columns describing where each fragment aligned
            curr_chr_count: how many chrs fragments from this read aligned to
            left_chr: Chromosome of location where the left part of this read was determined to come from
            left_cut_pos: Basepair location where the left part could have resulted from a cut at
            left_read_start: Basepair location where the left fragment could have started at
            left_pos: Basepair location where the left part of the read aligned
            left_orientation: orientation (+ or -) for whether the read extends to the left or the right away from the cut
            right_cut_pos
            right_pos
        """
        # first create arrays from the dict of positions
        curr_id_chrs = defaultdict(int) #dict of where fragments align chr->count
        pos_list = [None]*max_frags #chr and pos where each frag mapped
        val_list = [] #for printing the outline
        mapped_frag_len = 0 # length of array for which all frags are mapped frag (some reads are shorter, and don't have any fragments on the left side). I used the length so I can perform range and indexing operations (as opposed to keeping the index of the last mapped frag)

        val_counts = defaultdict(int)
        for val in curr_id_vals.values():
            val_counts[val] += 1

        for i in range(max_frags):
            val = ""
            if i in curr_id_vals:
                mapped_frag_len = i+1
                val = curr_id_vals[i]
                val_chr,val_pos,val_strand = curr_id_vals[i].split(" ")
                if val_counts[val] >= min_seen_frag_cutoff and val_chr != "*":
                    curr_id_chrs[val_chr] += 1
                    pos_list[i] = val

            val_list.append(val)

        all_list = [x for x in pos_list if x is not None]
        num_pos_all = len(set(all_list)) #could also be zero if completely unmapped

        #if there is only one mapped location
        if num_pos_all == 1:
            one_chr,one_pos,one_strand = all_list[0].split(" ")
            left_aln = -1
            right_aln = -1
            for idx in range(0,mapped_frag_len):
                val = pos_list[idx]
                if val is None:
                    continue
                this_chr,this_pos,this_strand = val.split(" ")
                if left_aln == -1:
                    left_aln = curr_id_aln_pos[idx]
                #just set the right one every time as we progress through the array
                right_aln = curr_id_aln_pos[idx]

            # if aligned in the forward orientation, the left will be less than right
            if left_aln < right_aln:
                read_start = left_aln
                read_end = right_aln + fragment_size
                orientation = '+'
            # if aligned in the reverse orientation, right is less than left, and the start of the alignment will be the left + fragment size
            else:
                read_start = left_aln + fragment_size
                read_end = right_aln
                orientation = '-'

            curr_annotation = "Linear"
            left_chr = one_chr
            left_cut_pos = read_start
            left_read_start = read_start
            left_pos = read_start
            left_orientation = orientation
            right_chr = one_chr
            right_cut_pos = read_end
            right_read_end = read_end
            right_pos = read_end
            right_orientation = orientation
            min_i = None

            outline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(curr_annotation,left_chr,left_cut_pos,left_read_start,left_pos,left_orientation,right_chr,right_cut_pos,right_read_end,right_pos,right_orientation,min_i,"\t".join(val_list))

            curr_chr_count = 1

            curr_loc = "%s:%s-%s"%(one_chr,read_start,read_end)
            curr_cuts = "NA"

            return outline,curr_chr_count,curr_loc,curr_cuts,curr_annotation,left_chr,left_cut_pos,left_read_start,left_orientation,right_chr,right_cut_pos,right_read_end,right_orientation

        #if there are multiple mapped locations
        else:
            min_pos_unique_sum = max_frags
            max_majority_sum = 0
            min_i = None #  keep track of the right-most i that maximizes the partition
            #create partitions of the fragments, choose the one that:
            # - minimizes the number of unique locations seen in each partition
            # - maximizes the number of fragments that are the majority in their respective partitions
    #        print('pos list: ' + str(pos_list))
    #        print('last frag ind: ' + str(mapped_frag_len))
            for i in range(1,mapped_frag_len):
                left_list = [x for x in pos_list[0:i] if x is not None]
                num_pos_left = len(set(left_list))

                right_list = [x for x in pos_list[i:mapped_frag_len] if x is not None]
                len_right_list = len(right_list)
                num_pos_right = len(set(right_list))

                #if one of the partitions has no pos (all are None), continue
                if num_pos_right == 0 or num_pos_left == 0:
                    continue

#                if val_list[0] == "chr1 223763905" and val_list[1] == "chr2 60498399":
#                    print("i:"+str(i)+"/"+str(mapped_frag_len))
#                    print('left bit: ' + str(pos_list[0:i]))
#                    print('right bit: ' + str(pos_list[i:mapped_frag_len]))
#                    print('num_pos_left: ' + str(num_pos_left))
#                    print('num_pos_right: ' + str(num_pos_right))
#                    print('ratio left: ' + str(ratio_left))
#                    print('ratio right: ' + str(ratio_right))
#                    if i > 5:
#                        toggle = 1
#                    if toggle == 1 and i < 5:
#                        asdf()

                this_unique_sum = num_pos_left + num_pos_right
                #if the total num of unique positions in each partition is minimized here, set this to be the best partition so far
                #or if it is a tie, calculate the majority ratio and choose the partitioning that maximizes that
                if this_unique_sum <= min_pos_unique_sum:
                    #calculate the ratio of the most-frequent position observed out of all positions in that partition for left and right partitions
                    left_counts = {}
                    left_val = None;
                    left_majority_count = 0
                    for val in left_list:
                        left_counts[val] = left_counts.get(val,0) + 1
                        if left_counts[val] > left_majority_count and val is not None:
                            left_majority_count = left_counts[val]
                            left_val = val

                    right_counts = {}
                    right_val = None;
                    right_majority_count = 0
                    for val in right_list:
                        right_counts[val] = right_counts.get(val,0) + 1
                        if right_counts[val] > right_majority_count and val is not None and val != left_val:
                            right_majority_count = right_counts[val]
                            right_val = val

                    this_majority_sum = left_majority_count + right_majority_count

                    #if the unique sum is less, then just set this i to min_i
                    if this_unique_sum < min_pos_unique_sum:
                        min_pos_unique_sum = this_unique_sum
                        max_majority_count = this_majority_sum
                        min_i = i
                    elif this_majority_sum > max_majority_sum: #this_unique_sum == min_pos_unique_sum from if statement above
                        min_pos_unique_sum = this_unique_sum
                        max_majority_sum = this_majority_sum
                        min_i = i

            if min_i is not None:
                left_counts = {}
                left_val = None;
                left_count = 0
                for val in pos_list[0:min_i]:
                    left_counts[val] = left_counts.get(val,0) + 1
                    if left_counts[val] > left_count and not val is None:
                        left_count = left_counts[val]
                        left_val = val
                if left_count > 0:
                    if left_counts[left_val] >= min_seen_frag_cutoff:
                        left_chr,left_pos,left_strand = left_val.split(" ")
                    else: #an alignment was seen, but below the cutoff
                        left_chr = "*"
                        left_pos = "-1"
                        left_strand='u'
                else: #no alignments were seen
                    left_chr = "*"
                    left_pos = "-2"
                    left_strand='u'

                right_counts = {}
                right_val = None;
                right_count = 0
                for val in pos_list[min_i:mapped_frag_len]:
                    right_counts[val] = right_counts.get(val,0) + 1
                    if right_counts[val] > right_count and not val is None:
                        right_count = right_counts[val]
                        right_val = val
                if right_count > 0:
                    if right_counts[right_val] >= min_seen_frag_cutoff:
                        right_chr,right_pos,right_strand = right_val.split(" ")
                    else: #an alignment was seen, but below the cutoff
                        right_chr = "*"
                        right_pos = "-1"
                        right_strand='u'
                else: #no alignments were seen
                    right_chr = "*"
                    right_pos = "-2"
                    right_strand='u'
            else: #no alignments were seen, min_i is None
                left_chr = "*"
                left_pos = "-3"
                left_strand='u'
                right_chr = "*"
                right_pos = "-3"
                right_strand='u'
                min_i = -1

            # estimate position where cuts may have occurred to create this read
            left_cut_pos = "*"
            left_read_start = left_pos
            if left_chr != "*":
                left_left_aln = -1; #for the left part of the read, identify the left and right parts of it -- the right part will be the cut position
                left_right_aln = -1;
                for idx in range(0,min_i):
                    val = pos_list[idx]
                    if val is None:
                        continue
                    this_chr,this_pos,this_strand = val.split(" ")
                    #make sure this fragment is aligned to the correct(identified) spot
                    if this_chr == left_chr and this_pos == left_pos:
                        if left_left_aln == -1:
                            left_left_aln = curr_id_aln_pos[idx]
                        #just set the left_right one every time as we progress through the array
                        left_right_aln = curr_id_aln_pos[idx]

                # if aligned in the forward orientation, the cut will be at the right of the last aligned position
                if left_left_aln < left_right_aln:
                    left_cut_pos = left_right_aln + fragment_size
                    left_read_start = left_left_aln
                # if aligned in the reverse orientation, the cut will be exactly at the left_right aligned position, but the start of the alignment will be the left_left + fragment size
                else:
                    left_cut_pos = left_right_aln
                    left_read_start = left_left_aln + fragment_size


            right_cut_pos = "*"
            right_read_end = right_pos
            if right_chr != "*":
                right_left_aln = -1; #for the right part of the read, identify the left and right parts of it -- the left part will be the cut position
                right_right_aln = -1;
                for idx in range(min_i,mapped_frag_len):
                    val = pos_list[idx]
                    if val is None:
                        continue
                    this_chr,this_pos,this_strand = val.split(" ")
                    #make sure this fragment is aligned to the correct(identified) spot
                    if this_chr == right_chr and this_pos == right_pos:
                        if right_left_aln == -1:
                            right_left_aln = curr_id_aln_pos[idx]
                        #just set the right_right one every time as we progress through the array
                        right_right_aln = curr_id_aln_pos[idx]

                # if aligned in the forward orientation, the cut will be at the right_left aligned position, and the alignment start will be right_right+fragment size
                if right_left_aln < right_right_aln:
                    right_read_end = right_right_aln + fragment_size
                    right_cut_pos = right_left_aln
                # if aligned in the reverse orientation, the cut will be at the right_left aligned position+fragment_size
                else:
                    right_cut_pos = right_left_aln + fragment_size
                    right_read_end = right_right_aln

            #create annotation for this read
            #if left read is to the left of the cut (forward alignment of read to the genome)
            left_anno = "%s-%s:%s"%(left_strand,left_chr,left_cut_pos)
            left_orientation = '-'
            #if left read is to the right of the cut
            if left_read_start > left_cut_pos:
                left_anno = "%s+%s:%s"%(left_strand,left_chr,left_cut_pos)
                left_orientation = '+'

            #if right read is to the right of the cut (forward alignment of the read to the genome)
            right_anno = "%s:%s+%s"%(right_chr,right_cut_pos,right_strand)
            right_orientation = '+'
            #if right read is to the left of the cut
            if right_read_end < right_cut_pos:
                right_anno = "%s:%s-%s"%(right_chr,right_cut_pos,right_strand)
                right_orientation = '-'

            curr_annotation = left_anno+"_"+right_anno

#            if curr_annotation == 'c+chr18:80075325_chr18:80075365-c':
#                print('left pos: ' + str(left_pos))
#                print('left start: ' + str(left_read_start))
#                print('left cut: ' + str(left_cut_pos))
#                print('right cut: ' + str(right_cut_pos))
#                print('right read_end: ' + str(right_read_end))
#                print('right pos: ' + str(right_pos))
#                print('left: ' + str(pos_list[0:min_i]))
#                print('left_left_aln: ' + str(left_left_aln))
#                print('left_right_aln: ' + str(left_right_aln))
#                print('right: ' + str(pos_list[min_i:mapped_frag_len]))
#                print('right_left_aln: ' + str(right_left_aln))
#                print('right_right_aln: ' + str(right_right_aln))


            outline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(curr_annotation,left_chr,left_cut_pos,left_read_start,left_pos,left_orientation,right_chr,right_cut_pos,right_read_end,right_pos,right_orientation,min_i,"\t".join(val_list))

            curr_chr_count = len(curr_id_chrs.keys())

            curr_loc = "%s:%s-%s~%s:%s-%s"%(left_chr,left_pos,left_cut_pos,right_chr,right_cut_pos,right_pos)
            curr_cuts = left_chr + ":" + str(left_cut_pos) + "~" + right_chr + ":" + str(right_cut_pos)

        return outline,curr_chr_count,curr_loc,curr_cuts,curr_annotation,left_chr,left_cut_pos,left_read_start,left_orientation,right_chr,right_cut_pos,right_read_end,right_orientation
        ### end of helper function

    frag_meta_file = root + ".meta.txt"

    r1_assignments_file = root + '.assignments_r1.txt'
    r2_assignments_file = root + '.assignments_r2.txt'
    frags_mapped_chrs = defaultdict(int)
    with open(frag_meta_file,"w") as frag_file, open(r1_assignments_file,'w') as af1, open(r2_assignments_file,'w') as af2:
        head = "ID\tcut_annotation\tleft_chr\tleft_cut\tleft_aln_start\tleft_pos\tleft_orientation\tright_chr\tright_cut\tright_aln_end\tright_pos\tright_orientation\tbreak_index"
        for i in range(max_frags):
            head += "\t"+str(i)
        frag_file.write(head+"\n")

        translocations = defaultdict(int) # hash of all translocation locations
        linear_count = 0
        translocation_count = 0
        large_deletion_count = 0
        large_inversion_count = 0
        unidentified_count = 0
        chroms_per_frag_read_count = defaultdict(int) # keep track for all read how many chrs the fragments mapped to
        frags_mapped_count = 0

        curr_id = "" # whole id of the read (alignments will be aggregated based on this id)
        curr_idx = "" # read index (either 1 for r1 or 2 for r2
        curr_id_vals = {} # keep track of where each frag maps to (this is a hash but it's like an array, so we can check for uninitialized values)
        curr_id_aln_pos = {} # keep track of the location where each frag aligns - curr_id_vals keeps track of where a read having a fragment in that position would align.
        with open(mapped_chopped_sam_file,'r') as aln_frags:
            for line in aln_frags:
                line_els = line.split("\t")
                line_chr = line_els[2]
                frags_mapped_chrs[line_chr] += 1
                frags_mapped_count += 1

                line_mapq = line_els[4]
                line_unmapped = int(line_els[1]) & 0x4
                line_reverse = int(line_els[1]) & 0x10
                line_start = int(line_els[3])-1
                line_strand = 'w'
                if line_reverse:
                    line_strand = 'c'

                match = re.match('(.*)\.CLR(\d)\.CLID(\d+)\.CLO(-?\d+)\.CLF(-?\d+)$',line_els[0])
                if not match:
                    raise Exception('Cannot parse id %s from line %s in %s. Perhaps line was trimmed?\n'%(line_els[0],line,mapped_chopped_sam_file))
                (next_id,read_idx,lungo_id,lungo_offset,lungo_frag) = match.groups()

                #write results from last id if this line is for a new id
                if curr_id != next_id and curr_id != "":
                    (outline,curr_chr_count,curr_loc,curr_cuts,curr_annotation,left_chr,left_cut,left_pos,left_orientation,right_chr,right_cut,right_pos,right_orientation) = analyze_curr_id_vals(curr_id_vals,curr_id_aln_pos,fragment_size,fragment_step_size,max_frags,min_seen_frag_cutoff)
                    frag_file.write(curr_id+"\t"+outline)

                    curr_classification = "NA"
                    if left_chr == "*" or right_chr == "*":
                        unidentified_count += 1
                        curr_classification = "Unidentified"
                    elif curr_annotation == 'Linear':
                        linear_count += 1
                        curr_classification = 'Linear'

                    else:
                        key = left_chr + " " + str(left_pos) + " " + right_chr + " " + str(right_pos)
                        translocations[key] += 1
                        if left_chr == right_chr:
                            if left_orientation == right_orientation:
                                if left_pos == right_cut or left_cut == right_pos:
                                    #vals = [outline,curr_chr_count,curr_loc,curr_cuts,curr_annotation,left_chr,left_cut,left_pos,left_orientation,right_chr,right_cut,right_pos,right_orientation]
                                    #names = ["outline","curr_chr_count","curr_loc","curr_cuts","curr_annotation","left_chr","left_cut","left_pos","left_orientation","right_chr","right_cut","right_pos","right_orientation"]
                                    #print("\n".join([str(x)+":"+str(y) for x,y in zip(names,vals)]))
                                    #print('currid: ' + curr_id)
                                    raise Exception('Left and right have same position!')
                                else:
                                    large_inversion_count += 1
                                    curr_classification = "Large inversion"
                            else:
                                large_deletion_count += 1
                                curr_classification = "Large deletion"
                        else:
                            translocation_count += 1
                            curr_classification = "Translocation"

                    chroms_per_frag_read_count[curr_chr_count] += 1
                    if int(curr_idx) == 1:
                        af1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'Fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))
                    else:
                        af2.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'Fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))

                    curr_id_vals = {}
                    curr_id_aln_pos = {}

                curr_id = next_id
                curr_idx = read_idx

                inferred_start = line_start - int(lungo_offset)
                if line_reverse:
                    inferred_start = line_start + int(lungo_offset)

                if int(line_mapq) > fragment_min_mapq:
                    curr_id_vals[int(lungo_frag)] = line_chr + " " + str(inferred_start) + " " + line_strand
                    curr_id_aln_pos[int(lungo_frag)] = line_start
                else:
                    curr_id_vals[int(lungo_frag)] = "*" + " " + str(inferred_start) + " " + line_strand
                    curr_id_aln_pos[int(lungo_frag)] = line_start

            #finish last item
            if frags_mapped_count > 0:
                (outline,curr_chr_count,curr_loc,curr_cuts,curr_annotation,left_chr,left_cut,left_pos,left_orientation,right_chr,right_cut,right_pos,right_orientation) = analyze_curr_id_vals(curr_id_vals,curr_id_aln_pos,fragment_size,fragment_step_size,max_frags,min_seen_frag_cutoff)
                frag_file.write(curr_id+"\t"+outline)

                curr_classification = "NA"
                if left_chr == "*" or right_chr == "*":
                    unidentified_count += 1
                    curr_classification = "Unidentified"
                elif curr_annotation == 'Linear':
                    linear_count += 1
                    curr_classification = 'Linear'
                else:
                    key = left_chr + " " + str(left_cut) + " " + right_chr + " " + str(right_cut)
                    translocations[key] += 1
                    if left_chr == right_chr:
                        if left_orientation == right_orientation:
                            large_inversion_count += 1
                            curr_classification = "Large inversion"
                        else:
                            large_deletion_count += 1
                            curr_classification = "Large deletion"
                    else:
                        translocation_count += 1
                        curr_classification = "Translocation"

                chroms_per_frag_read_count[curr_chr_count] += 1

                if int(curr_idx) == 1:
                    af1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'Fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))
                else:
                    af2.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'Fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))
            #done with last one

    logging.info("Found %d linear reads, %d translocations, %d large deletions, %d large inversions, and %d unidentified reads"%(linear_count,translocation_count,large_deletion_count,large_inversion_count,unidentified_count))

    frags_aligned_chrs_root = root + ".chrs"
    keys = sorted(frags_mapped_chrs.keys())
    if len(keys) > 0:
        vals = [frags_mapped_chrs[key] for key in keys]
        with open(frags_aligned_chrs_root+".txt","w") as fout:
            fout.write("chr\tnumReads\n")
            for key in keys:
                fout.write(str(key) + '\t' + str(frags_mapped_chrs[key]) + '\n')
        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        ax.bar(range(len(keys)),vals,tick_label=keys)
        ax.set_ymargin(0.05)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Location of fragments from unaligned reads')

        plt.savefig(frags_aligned_chrs_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(frags_aligned_chrs_root+".png",pad_inches=1,bbox_inches='tight')

    frag_chroms_per_read_root = root + ".alignedChromsPerRead"
    keys = sorted(chroms_per_frag_read_count.keys())
    if len(keys) > 0:
        vals = [chroms_per_frag_read_count[key] for key in keys]
        with open(frag_chroms_per_read_root+".txt","w") as fout:
            fout.write("numChroms\tnumReads\n")
            for key in keys:
                fout.write(str(key) + '\t' + str(chroms_per_frag_read_count[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        ax.bar(range(len(keys)),vals,tick_label=keys)
        ax.set_ymargin(0.05)
        ax.set_xlabel('Number of Chromosomes')
        ax.set_ylabel('Number of Reads')
        ax.set_title('Number of different chromosomes aligned by a fragmented single read')

        plt.savefig(frag_chroms_per_read_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(frag_chroms_per_read_root+".png",pad_inches=1,bbox_inches='tight')

    # make translocation table
    # dict of count of reads aligning from chrA to chrB
    translocation_table = {}
    translocation_report_file = root + ".translocationReport.txt"
    with open(translocation_report_file,"w") as fout:
        fout.write("from\tto\tcount\n")
        for key in sorted(translocations.keys()):
            (from_chr,from_pos,to_chr,to_pos) = key.split(" ")
            fout.write(from_chr+'\t'+to_chr+ '\t' + str(translocations[key]) + '\n')
            if from_chr not in translocation_table:
                translocation_table[from_chr] = {}
            if to_chr not in translocation_table[from_chr]:
                translocation_table[from_chr][to_chr] = 0
            translocation_table[from_chr][to_chr] += translocations[key]

    translocation_report_table_file = root + ".translocationReport.table.txt"
    with open(translocation_report_table_file,"w") as fout:
        chrs = sorted(translocations.keys())
        translocation_head = "data\t"+"\t".join(chrs)
        fout.write(translocation_head+"\n")
        for key in chrs:
            line = key
            for key2 in chrs:
                val = 0
                if key in translocation_table and key2 in translocation_table[key]:
                    val = translocation_table[key][key2]
                line += "\t"+str(val)
            fout.write(line+"\n")

    if unmapped_reads_fastq_r2 is None:
        os.remove(r2_assignments_file)
        r2_assignments_file = None

    with open(info_file,'w') as fout:
        fout.write("\t".join(['r1_assignments_file','r2_assignments_file','linear_count','translocation_count','large_deletion_count','large_inversion_count','unidentified_count','read_total'])+"\n")
        fout.write("\t".join([str(x) for x in [r1_assignments_file,r2_assignments_file,linear_count,translocation_count,large_deletion_count,large_inversion_count,unidentified_count,frags_mapped_count]])+"\n")
        frags_plot_obj_str = "None"
        if frags_plot_obj is not None:
            frags_plot_obj_str = frags_plot_obj.to_json()
        fout.write(frags_plot_obj_str+"\n")

    return (r1_assignments_file,r2_assignments_file,linear_count,translocation_count,large_deletion_count,large_inversion_count,unidentified_count,frags_plot_obj)

def prep_unmapped_assignments(root,unmapped_reads_fastq_r1,unmapped_reads_fastq_r2,can_use_previous_analysis=False):
    """
    Creates assigmnets file for unmapped reads

    params:
        root: root for written files
        unmapped_reads_fastq_r1: fastq file of r1 unmappable reads
        unmapped_reads_fastq_r2: fastq file of r2 unmappable reads
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch


    returns:
        r1_assignments_file: file of location assignments for each r1 read
        r2_assignments_file: file of location assignments for each r2 read (None if single-end)
        unmapped_r1_count: number of unmapped R1 reads
        unmapped_r2_count: number of unmapped R2 reads
    """
    info_file = root + '.info'
    if os.path.isfile(info_file) and can_use_previous_analysis:
        unmapped_r1_count = 0
        unmapped_r2_count = 0
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 4:
                (r1_assignments_file_str,r2_assignments_file_str,unmapped_r1_count_str,unmapped_r2_count_str) = line_els
                r1_assignments_file = None if r1_assignments_file_str == "None" else r1_assignments_file_str
                r2_assignments_file = None if r2_assignments_file_str == "None" else r2_assignments_file_str
                unmapped_r1_count  = int(unmapped_r1_count_str)
                unmapped_r2_count  = int(unmapped_r2_count_str)

                if unmapped_r1_count >= 0 or unmapped_r2_count >= 0:
                    logging.info('Using previously-created assignments for %s unmapped R1 reads and %s unmapped R2 reads'%(unmapped_r1_count, unmapped_r2_count))
                    return (r1_assignments_file,r2_assignments_file,unmapped_r1_count, unmapped_r2_count)
                else:
                    logging.info('Could not recover previously-created assignments. Recreating.')

    logging.info('Creating assignments for unmapped reads')

    unmapped_str = "%s\t%s\t%s\t%s\t%s"%('Unmapped','Unidentified','u+*:*','*:*','*:*')

    r1_assignments_file = root + '.assignments_r1.txt'
    r2_assignments_file = None
    unmapped_r1_count = 0
    unmapped_r2_count = 0

    if unmapped_reads_fastq_r1.endswith('.gz'):
        f1_in = io.BufferedReader(gzip.open(unmapped_reads_fastq_r1,'rt'))
    else:
        f1_in = open(unmapped_reads_fastq_r1,'rt')


    with open(r1_assignments_file,'w') as af1:
        while (1):
            f1_id_line   = f1_in.readline().strip()
            f1_seq_line  = f1_in.readline().strip()
            f1_plus_line = f1_in.readline()
            f1_qual_line = f1_in.readline().strip()
            if not f1_qual_line : break

            f1_id = f1_id_line[1:]

            af1.write("%s\t%s\n"%(f1_id,unmapped_str))
            unmapped_r1_count += 1
    f1_in.close()

    if unmapped_reads_fastq_r2 is not None:
        r2_assignments_file = root + '.assignments_r2.txt'

        if unmapped_reads_fastq_r2.endswith('.gz'):
            f2_in = io.BufferedReader(gzip.open(unmapped_reads_fastq_r2,'rt'))
        else:
            f2_in = open(unmapped_reads_fastq_r2,'rt')


        with open(r2_assignments_file,'w') as af2:
            while (1):
                f2_id_line   = f2_in.readline().strip()
                f2_seq_line  = f2_in.readline()
                f2_plus_line = f2_in.readline()
                f2_qual_line = f2_in.readline()
                if not f2_qual_line : break

                f2_id = f2_id_line[1:]

                af2.write("%s\t%s\n"%(f2_id,unmapped_str))
                unmapped_r2_count += 1
        f2_in.close()


    with open(info_file,'w') as fout:
        fout.write("\t".join(['r1_assignments_file','r2_assignments_file','unmapped_r1_count','unmapped_r2_count'])+"\n")
        fout.write("\t".join([str(x) for x in [r1_assignments_file,r2_assignments_file,unmapped_r1_count,unmapped_r2_count]])+"\n")


    logging.info('Finished creation of assignments for %s unaligned R1 reads and %s unaligned R2 reads'%(unmapped_r1_count, unmapped_r2_count))

    return (r1_assignments_file,r2_assignments_file,unmapped_r1_count,unmapped_r2_count)


def prep_crispresso2(root,input_fastq_file,read_ids_for_crispresso,cut_counts_for_crispresso,av_read_length,genome,genome_len_file,crispresso_min_count,crispresso_min_aln_score,samtools_command,crispresso_command):
    """
    Prepares reads for analysis by crispresso2
    Frequently-aligned locations with a min number of reads (crispresso_min_count) are identified
    Reads from these locations are extracted and prepared for analysis by CRISPResso2

    params:
        root: root for written files
        input_fastq_file: input fastq with read sequences
        read_ids_for_crispresso: dict of readID=>cut assignment (cut assignment = cut type, cut1loc, cut2loc)
        cut_counts_for_crispresso: dict of cut=>count for # reads assigned to each cut
        av_read_length: average read length (used to calculate reference region size for CRISPResso)
        genome: path to genome fa file
        genome_len_file: path to tab-sep file with lengths of chrs in genome (ends with .fai)
        crispresso_min_count: min number of reads at site for crispresso2 processing
        crispresso_min_aln_score: minimum score for reads to align to amplicons
        samtools_command: location of samtools to run
        crispresso_command: location of crispresso to run

    returns:
        crispresso_infos: dict containing metadata about each crispresso run
            #dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command
    """
    logging.info('Preparing reads for analysis by CRISPResso2')

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
        f_in = io.BufferedReader(gzip.open(input_fastq_file,'rb'))
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

        for this_cut in read_ids_for_crispresso[read_id]:

            if '*' in this_cut:
                continue

            if cut_counts_for_crispresso[this_cut] < crispresso_min_count:
                continue

            if this_cut not in cut_names:
                cut_name = this_cut.replace(" ","_").replace("~","_").replace(':','_').replace("+","p").replace("-","m")
                cut_names[this_cut] = cut_name
            cut_name = cut_names[this_cut]

            printed_read_count += 1
            reads_file = os.path.join(data_dir,cut_name+".fq")
            if reads_file not in filehandles:
                fh = open(reads_file,'w')
                filehandles[reads_file] = fh
            filehandles[reads_file].write("%s\n%s\n%s%s\n"%(id_line,seq_line,plus_line,qual_line))

            printed_cuts[this_cut] += 1


    amp_half_length = av_read_length/1.5
    for i,cut in enumerate(printed_cuts):
        cut_els = cut.split(" ")
        cut_type = " ".join(cut_els[0:-1]) #Large deltion is two words..
        cut_loc = cut_els[-1]
        cut_name = cut_names[cut]
        amp_seq = None
        if cut_type == "Linear":
            cut_chr,cut_pos = cut_loc.split(":")
            amp_start = int(cut_pos)-amp_half_length
            amp_end = int(cut_pos)+amp_half_length
            amp_seq = subprocess.check_output(
                '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,cut_chr,amp_start,amp_end),shell=True).decode(sys.stdout.encoding).strip()
        else:
            left_strand = cut_loc[0]
            left_direction = cut_loc[1]
            left_cut,right_cut = cut_loc[2:-2].split("~")
            left_chr,left_pos = left_cut.split(":")
            left_pos = int(left_pos)
            right_chr,right_pos = right_cut.split(":")
            right_pos = int(right_pos)
            right_strand = cut_loc[-1]
            right_direction = cut_loc[-2]

            cut_loc = cut_loc[2:-2]

            if left_direction == "-":
                left_amp_end = left_pos
                left_amp_start = left_pos - amp_half_length
                left_amp_seq = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,left_chr,left_amp_start,left_amp_end - 1),shell=True).decode(sys.stdout.encoding).strip()
                if left_direction == "c":
                    left_amp_seq = complement(left_amp_seq)

            elif left_direction == "+":
                left_amp_end = left_pos + amp_half_length
                left_amp_start = left_pos

                left_amp_seq = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,left_chr,left_amp_start,left_amp_end),shell=True).decode(sys.stdout.encoding).strip()
                if left_direction == "w":
                    left_amp_seq = reverse(left_amp_seq)
                else:
                    left_amp_seq = reverse_complement(left_amp_seq)
            else:
                raise Exception('Got unexpected cut direction: ' + left_direction + ' from cut ' + cut)

            if right_direction == "+":
                right_amp_start = right_pos
                right_amp_end = right_pos + amp_half_length
                right_amp_seq = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,right_chr,right_amp_start,right_amp_end),shell=True).decode(sys.stdout.encoding).strip()
                if right_direction == "c":
                    right_amp_seq = complement(right_amp_seq)

            elif right_direction == "-":
                right_amp_end = right_pos
                right_amp_start = right_pos - amp_half_length

                right_amp_seq = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,right_chr,right_amp_start,right_amp_end - 1),shell=True).decode(sys.stdout.encoding).strip()
                if right_direction == "w":
                    right_amp_seq = reverse(right_amp_seq)
                else:
                    right_amp_seq = reverse_complement(right_amp_seq)

            else:
                raise Exception('Got unexpected cut direction: ' + right_direction + ' from cut ' + cut)

#            print('left amp: ' +left_amp_seq + ' right: ' + right_amp_seq)
#            asdf()
            amp_seq = left_amp_seq + right_amp_seq

        reads_file = os.path.join(data_dir,cut_name+".fq")
        output_folder = os.path.join(data_dir,'CRISPResso_on_'+cut_names[cut])
        crispresso_cmd = "%s -o %s -n %s --default_min_aln_score %d -a %s -r1 %s &> %s.log"%(crispresso_command,data_dir,cut_name,crispresso_min_aln_score,amp_seq,reads_file,reads_file)

        crispresso_infos.append({
                "cut":cut,
                "name":cut_names[cut],
                "type": cut_type,
                "cut_loc": cut_loc,
                "amp_seq": amp_seq,
                "output_folder":output_folder,
                "reads_file": reads_file,
                "printed_read_count": printed_cuts[cut],
                "command": crispresso_cmd
                })

    return (crispresso_infos)

def run_and_aggregate_crispresso(root,crispresso_infos,n_processes,min_count_to_run_crispresso,skip_failed=True,can_use_previous_analysis=False):
    """
    Runs CRISPResso2 commands and aggregates output

    params:
        root: root for written files
        crispresso_infos: array of metadata information for CRISPResso
            dict of: cut, name, type, cut_loc, amp_seq, output_folder, reads_file, printed_read_count, command
        n_processes: number of processes to run CRISPResso commands on
        min_count_to_run_crispresso: minimum number of reads assigned to a site to run CRISPResso
        skip_failed: if true, failed CRISPResso runs are skipped. Otherwise, if one fails, the program will fail.
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        crispresso_results: dict of run_names and run_sub_htmls for display in report

    """

    crispresso_stats_file = root + ".info"

    crispresso_commands = []
    for crispresso_info in crispresso_infos:
        if crispresso_info['printed_read_count'] >= min_count_to_run_crispresso:
            crispresso_commands.append(crispresso_info['command'])
    expected_completed_command_count = len(crispresso_commands)

    #check to see if runs have been completed previously
    crispresso_run_names = []
    crispresso_sub_htmls = {}
    read_completed_command_count = 0
    if os.path.isfile(crispresso_stats_file) and can_use_previous_analysis:
        with open(crispresso_stats_file,'r') as fin:
            head_line = fin.readline()
            for line in fin:
                line_els = line.strip().split()
                crispresso_run_names.append(line_els[0])
                crispresso_sub_htmls[line_els[0]] = line_els[1]
                read_completed_command_count += 1
        if read_completed_command_count == expected_completed_command_count:
            logging.info('Using %d previously-completed CRISPResso runs'%read_completed_command_count)
            crispresso_results = {}
            crispresso_results['run_names'] = crispresso_run_names
            crispresso_results['run_sub_htmls'] = crispresso_sub_htmls
            return(crispresso_results)
        else:
            logging.info('Could not recover previously-completed CRISPResso runs. Rerunning.')


    #otherwise run crispresso runs
    logging.info('Running and analyzing ' + str(len(crispresso_infos)) + ' alignments using CRISPResso2')


    CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_commands,n_processes,'run',skip_failed)

    crispresso_results = {}
    crispresso_results['run_names'] = []
    crispresso_results['run_sub_htmls'] = {}
    crispresso_command_file = root + '.commands.txt'
    crispresso_info_file = root + '.summary.txt'
    with open(crispresso_info_file,'w') as crispresso_info_file, open(crispresso_command_file,'w') as crispresso_command_file:
        crispresso_info_file.write('name\tcut_annotation\tcut_type\tcut_location\treference\treads_printed\tn_total\treads_aligned\treads_unmod\treads_mod\treads_discarded\treads_insertion\treads_deletion\treads_substitution\treads_only_insertion\treads_only_deletion\treads_only_substitution\treads_insertion_and_deletion\treads_insertion_and_substitution\treads_deletion_and_substitution\treads_insertion_and_deletion_and_substitution\tamplicon_sequence\n')
        for crispresso_info in crispresso_infos:
            crispresso_command_file.write(crispresso_info['command'])
            name = crispresso_info['name']
            cut = crispresso_info['cut']
            cut_type = crispresso_info['type']
            cut_loc = crispresso_info['cut_loc']
            cut_amp_seq = crispresso_info['amp_seq']
            n_printed = crispresso_info['printed_read_count']
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

            try:
                run_data = CRISPRessoShared.load_crispresso_info(crispresso_info['output_folder'])

            except:
                logging.debug('Could not load CRISPResso run information from ' + crispresso_info['name'] + ' at ' + crispresso_info['output_folder'])
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

                else:
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

                unmod_pct = "NA"
                mod_pct = "NA"
                if n_aligned > 0:
                    unmod_pct = 100*n_unmod/float(n_aligned)
                    mod_pct = 100*n_mod/float(n_aligned)
            new_vals = [name,cut,cut_type,cut_loc,ref_name,n_printed,n_total,n_aligned,n_unmod,n_mod,n_discarded,n_insertion,n_deletion,n_substitution,n_only_insertion,n_only_deletion,n_only_substitution,n_insertion_and_deletion,n_insertion_and_substitution,n_deletion_and_substitution,n_insertion_and_deletion_and_substitution,cut_amp_seq]
            crispresso_info_file.write("\t".join([str(x) for x in new_vals])+"\n")
        with open(crispresso_stats_file,'w') as fout:
            fout.write('Run Name\tRun sub-html\n')
            for run_name in crispresso_results['run_names']:
                fout.write(run_name + "\t" + crispresso_results['run_sub_htmls'][run_name]+"\n")
        return crispresso_results


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


def sam2aln(ref_seq,sam_line,include_ref_surrounding_read = True):
    """
    Creates an alignment of nucleotides from a sam line and a reference sequence by parsing the cigar string

    params:
        ref_seq: string, reference sequence
        sam_line: tab-sep string, sam line. Expected entries are alignment position (element 4) cigar (element 6), and sequence (element 10)
        include_ref_surrounding_read: boolean, whether to include ref sequence to the left or right of the sequenced read. If false, only the reference directly overlapping the read is included

    returns:
        ref_str: reference string with gaps added as appropriate
        aln_str: reference string with clipped bases removed and gap added as appropriate
        clipped_left_bp: number of bp clipped from left side of read - includes both hard and soft-clipped bases
        clipped_right_bp: number of bp clipped from right side of read

    tests from samtools spec:
    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*',False)
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r002	0	ref	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r002	0	ref	9	30	3S6M1P1I4M5S	*	0	0	AAAAGATAAGGATAGGGGG	*')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    """
    sam_line_els = sam_line.split("\t")
    cigar = sam_line_els[5]

    cigar_pattern = re.compile("(\d+)(\w)")
    cigar_els = []
    for (c_count,c_char) in re.findall(cigar_pattern,cigar):
        cigar_els.append((int(c_count),c_char))

    remaining_ref = ref_seq
    remaining_aln = sam_line_els[9]

    aln_str = ""
    ref_str = ""
    curr_start = int(sam_line_els[3])-1
    if include_ref_surrounding_read:
        aln_str = "-"*curr_start
        ref_str = ref_seq[0:curr_start]

    remaining_ref = ref_seq[curr_start:]

#    print('aln_str: ' + aln_str)
#    print('ref_str: ' + ref_str)
#    print('remaining_aln: ' + remaining_aln)
#    print('remaining_ref: ' + remaining_ref)

    clipped_left_bp = 0 #keep track of hard clipping
    clipped_right_bp = 0

    for idx,(c_count,c_char) in enumerate(cigar_els):
#        print('curr operation: ' + str(c_count) + ' ' + c_char)
        if c_char == 'M':
            aln_str += remaining_aln[0:c_count]
            remaining_aln = remaining_aln[c_count:]
            ref_str += remaining_ref[0:c_count]
            remaining_ref = remaining_ref[c_count:]
        elif c_char == 'I':
            aln_str += remaining_aln[0:c_count]
            remaining_aln = remaining_aln[c_count:]
            ref_str += '-'*c_count
        elif c_char == 'D' or c_char == 'N':
            aln_str += '-'*c_count
            ref_str += remaining_ref[0:c_count]
            remaining_ref = remaining_ref[c_count:]
        elif c_char == 'S':
            remaining_aln = remaining_aln[c_count:]
            if idx == 0:
                clipped_left_bp += c_count
            if idx == len(cigar_els)-1:
                clipped_right_bp += c_count
        elif c_char == 'H' or c_char == 'P':
            if idx == 0:
                clipped_left_bp += c_count
            if idx == len(cigar_els)-1:
                clipped_right_bp += c_count


            pass
        else:
            raise Exception('Unable to parse cigar character: ' + c_char)
#        print('aln_str: ' + aln_str)
#        print('ref_str: ' + ref_str)
#        print('remaining_aln: ' + remaining_aln)
#        print('remaining_ref: ' + remaining_ref)


    if include_ref_surrounding_read:
        aln_str += '-'*len(remaining_ref)
        ref_str += remaining_ref

#    print('Final')
#    print('aln_str: ' + aln_str)
#    print('ref_str: ' + ref_str)

    return(ref_str,aln_str,clipped_left_bp,clipped_right_bp)

def getMatchLeftRightOfCut(ref_seq,sam_line,cut_pos,debug=False):
    """
    Gets the number of mismatches at the beginning and the end of the read specified in the sam line. Left and right are defined by the cut position

    params:
        ref_seq: string, reference sequence
        sam_line: tab-sep string, sam line. Expected entries are alignment position (element 4) cigar (element 6), and sequence (element 10)
        cut_pos: the position of the cut that defines left and right. Technically, this cut happens after this many bases.
        debug: boolean whether to print debug

    returns:
        left_aln_bases_count: number of bases from this read that were aligned to the left part
        left_ref_bases_count: number of bases from the ref that were aligned to the left part
        left_all_match_count: number of total bases on the left part that matched
        left_start_match_count: number of bases at the beginning of the read that matched exactly
        right_aln_bases_count
        right_ref_bases_count
        right_all_match_count
        right_start_match_count
    """

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln(ref_seq,sam_line)
    ref_str = ref_str.upper()
    aln_str = aln_str.upper()

    if debug:
        print('cut pos: ' + str(cut_pos))
        print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    left_all_match_count = 0
    left_start_match_count = 0
    left_read_bases_count = 0
    left_ref_bases_count = 0

    seen_read = False
    seen_mismatch = False
    aln_pos = 0 # index in alignment
    ref_pos = 0 # index of ref sequence (to check against cut_pos)
    while ref_pos < cut_pos:
        if ref_str[aln_pos] != '-':
            ref_pos += 1
        if aln_str[aln_pos] != '-':
            seen_read = True
            left_read_bases_count += 1
        if seen_read:
            if ref_str[aln_pos] == aln_str[aln_pos]:
                left_all_match_count += 1
                if not seen_mismatch:
                    left_start_match_count += 1
            else:
                seen_mismatch = True
            if ref_str[aln_pos] != '-':
                left_ref_bases_count += 1
        aln_pos += 1

    if debug:
        print('left all match count: ' + str(left_all_match_count))
        print('left start match count: ' + str(left_start_match_count))
        print('left read bases count: ' + str(left_read_bases_count))
        print('left ref bases count: ' + str(left_ref_bases_count))

    right_all_match_count = 0
    right_start_match_count = 0
    right_read_bases_count = 0
    right_ref_bases_count = 0

    seen_read = False
    seen_mismatch = False
    aln_pos = len(aln_str)-1 # index in alignment
    ref_pos = len(ref_str.replace("-","")) # index of ref sequence (to check against cut_pos)
    while ref_pos > cut_pos:
        if ref_str[aln_pos] != '-':
            ref_pos -= 1
        if aln_str[aln_pos] != '-':
            seen_read = True
            right_read_bases_count += 1
        if seen_read:
            if ref_str[aln_pos] == aln_str[aln_pos]:
                right_all_match_count += 1
                if not seen_mismatch:
                    right_start_match_count += 1
            else:
                seen_mismatch = True
            if ref_str[aln_pos] != '-':
                right_ref_bases_count += 1
        aln_pos -= 1

    if debug:
        print('right all match count: ' + str(right_all_match_count))
        print('right start match count: ' + str(right_start_match_count))
        print('right read bases count: ' + str(right_read_bases_count))
        print('right ref bases count: ' + str(right_ref_bases_count))

    if clipped_left_bp > 0: #if the left side of the read was soft/hard clipped
        if left_read_bases_count > 0: #if any part of this read was on the left side
            left_read_bases_count += clipped_left_bp
            left_start_match_count = 0
        else: #if this read was actually all on the right side
            right_read_bases_count += clipped_left_bp
    if clipped_right_bp > 0: #if the right side of the read was soft/hard clipped
        if right_read_bases_count > 0: #if any part of this read was on the right side
            right_read_bases_count += clipped_right_bp
            right_start_match_count = 0
        else: #if this read was actually all on the left side
            left_read_bases_count += clipped_right_bp

    return (left_read_bases_count, left_ref_bases_count, left_all_match_count, left_start_match_count,
        right_read_bases_count, right_ref_bases_count, right_all_match_count, right_start_match_count)



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
  padding-top: 40px;
  padding-bottom: 40px;
  background-color: #f5f5f5;
}

</style>
<div class='container'>
<div class='row justify-content-md-center'>
<div class='col-8'>
    <div class='text-center pb-4'>
    <h1 class='display-3'>CRISPRlungo</h1><hr><h2>"""+report_name+"""</h2>
    </div>
"""
        data_path = ""
        if len(crispresso_run_names) > 0:
            run_string = """<div class='card text-center mb-2'>
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

        for plot_obj in ordered_plot_objects:
            plot_str = "<div class='card text-center mb-2'>\n\t<div class='card-header'>\n"
            plot_str += "<h5>"+plot_obj.title+"</h5>\n"
            plot_str += "</div>\n"
            plot_str += "<div class='card-body'>\n"
            plot_str += "<a href='"+data_path+plot_obj.name+".pdf'><img src='"+data_path + plot_obj.name + ".png' width='80%' ></a>\n"
            plot_str += "<label>"+plot_obj.label+"</label>\n"
            for (plot_data_label,plot_data_path) in plot_obj.datas:
                plot_str += "<p class='m-0'><small>Data: <a href='"+data_path+plot_data_path+"'>" + plot_data_label + "</a></small></p>\n";
            plot_str += "</div></div>\n";
            html_str += plot_str

        html_str += """
</div>
</div>
</div>
     </body>
</html>
"""
        with open(report_file,'w') as fo:
            fo.write(html_str)
        logging.info('Wrote ' + report_file)

main()
