import argparse
from collections import defaultdict
import gzip
import io
import logging
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import subprocess
import sys

#todo: validate with R2


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

    primer_chr = ""
    primer_loc = -1
    if 'primer_site' in settings:
        primer_info = settings['primer_site'].split("_")
        primer_chr = primer_info[0]
        primer_loc = int(primer_info[1])
    primer_seq = None
    if 'primer_seq' in settings:
        primer_seq = settings['primer_seq']

    target_padding = settings['alignment_extension']
    target_length = av_read_length
    if target_padding < 0:
        target_length += target_padding
        target_padding = 0
        target_padding += av_read_length

    (target_names,target_info) = make_artificial_targets(
            cuts=settings['cuts'],
            genome=settings['genome'],
            target_length=target_length,
            target_padding=target_padding,
            primer_chr=primer_chr,
            primer_loc=primer_loc,
            primer_seq=primer_seq,
            add_non_primer_cut_targets=settings['add_non_primer_cut_targets'],
            samtools_command=settings['samtools_command'])

    custom_index_fasta = make_target_index(
                    root = settings['root'],
                    target_names = target_names,
                    target_info = target_info,
                    bowtie2_command=settings['bowtie2_command'],
                    bowtie2_threads = settings['bowtie2_threads'],
                    )

    reads_to_align_r1 = settings['fastq_r1'] #if alignment to genome happens first, the input for artificial target mapping will be reads that don't align to the genome
    reads_to_align_r2 = settings['fastq_r2']

    #if umis are provided, add them to the fastqs
    if settings['fastq_umi']:
        umi_r1, umi_r2 = add_umi_from_umi_file(
            root = settings['root']+'.addUMI',
            fastq_r1 = reads_to_align_r1,
            fastq_r2 = reads_to_align_r2,
            fastq_umi = settings['fastq_umi']
        )
        reads_to_align_r1 = umi_r1
        reads_to_align_r2 = umi_r2

    if settings['dedup_input']:
        dedup_r1, dedup_r2 = dedup_file(
            root = settings['root']+'.dedup',
            fastq_r1 = reads_to_align_r1,
            fastq_r2 = reads_to_align_r2,
            umi_regex = settings['umi_regex']
        )
        reads_to_align_r1 = dedup_r1
        reads_to_align_r2 = dedup_r2

    custom_aligned_count = 0
    crispresso_commands = [] #list of crispresso commands to run -- from custom targets and the genomic alignments
    crispresso_infos = [] #meta info about the crispresso runs
            #tuple of: name, chr, start, end, readCount, amplicon

    #keep track of where every read is assigned
    #these lists keep track of the files containing assignments from the separate analyses
    r1_assignment_files = []
    r2_assignment_files = []

    (genome_r1_assignments,genome_r2_assignments,genome_unmapped_r1, genome_unmapped_r2, genome_aligned_count, genome_mapped_bam_file,genome_chr_aln_plot_obj,genome_tlen_plot_object
        ) = align_reads(
                root = settings['root']+'.genomeAlignment',
                fastq_r1 = reads_to_align_r1,
                fastq_r2 = reads_to_align_r2,
                bowtie2_reference = settings['bowtie2_genome'],
                reference_name = 'genome',
                target_info = target_info,
                bowtie2_command = settings['bowtie2_command'],
                bowtie2_threads = settings['bowtie2_threads'],
                samtools_command = settings['samtools_command'],
                keep_intermediate = settings['keep_intermediate']
                )
    r1_assignment_files.append(genome_r1_assignments)
    if genome_r2_assignments is not None:
        r2_assignment_files.append(genome_r2_assignments)

    if len(target_names) > 0:
        #first, align the unaligned r1s
        logging.info('Processing unaligned R1 reads')
        (custom_r1_assignments,none_file,custom_unmapped_r1, none_file2, custom_r1_aligned_count, custom_r1_mapped_bam_file,custom_r1_chr_aln_plot_obj,none_custom_tlen_plot_object
            ) = align_reads(
                    root = settings['root']+'.r1.customTargetAlignment',
                    fastq_r1 = genome_unmapped_r1,
                    fastq_r2 = None,
                    bowtie2_reference=custom_index_fasta,
                    reference_name='custom targets',
                    target_info = target_info,
                    bowtie2_command=settings['bowtie2_command'],
                    bowtie2_threads=settings['bowtie2_threads'],
                    samtools_command=settings['samtools_command'],
                    keep_intermediate = settings['keep_intermediate']
                    )
        r1_assignment_files.append(custom_r1_assignments)
        #if input is paired, align second reads to the custom targets as welll
        if genome_unmapped_r2 is not None:
            logging.info('Processing unaligned R2 reads')
            (custom_r2_assignments,none_file,custom_unmapped_r2, none_file2, custom_r2_aligned_count, custom_r2_mapped_bam_file,custom_r2_chr_aln_plot_obj,none_custom_tlen_plot_object
                ) = align_reads(
                        root = settings['root']+'.r2.customTargetAlignment',
                        fastq_r1 = genome_unmapped_r2,
                        fastq_r2 = None,
                        bowtie2_reference=custom_index_fasta,
                        reference_name='custom targets',
                        target_info = target_info,
                        bowtie2_command=settings['bowtie2_command'],
                        bowtie2_threads=settings['bowtie2_threads'],
                        samtools_command=settings['samtools_command'],
                        keep_intermediate = settings['keep_intermediate']
                        )
            r2_assignment_files.append(custom_r2_assignments)

    #chop reads
    (frag_r1_assignments,frag_r2_assignments,translocation_count,large_deletion_count,unidentified_count,frags_plot_obj
        ) = chop_reads(
                root = settings['root']+".frags",
                unmapped_reads_fastq_r1 = custom_unmapped_r1,
                unmapped_reads_fastq_r2 = custom_unmapped_r2,
                bowtie2_genome = settings['bowtie2_genome'],
                fragment_size = settings['fragment_size'],
                fragment_step_size = settings['fragment_step_size'],
                bowtie2_command = settings['bowtie2_command'],
                bowtie2_threads = settings['bowtie2_threads'],
                samtools_command = settings['samtools_command'],
                )
    r1_assignment_files.append(frag_r1_assignments)
    if frag_r2_assignments is not None:
        r2_assignment_files.append(frag_r2_assignments)

    final_assignment_file,r1_read_ids_for_crispresso = make_final_read_assignments(
            root = settings['root']+'.final',
            r1_assignment_files = r1_assignment_files,
            r2_assignment_files = r2_assignment_files,
            cut_sites = cut_sites,
            target_info = target_info
            )

    quit()

    # Prepare CRISPResso runs
    (genome_crispresso_infos, genome_crispresso_commands,genome_read_count_at_cuts
        ) = prep_crispresso2_global(
                root = settings['root'],
                cuts = settings['cuts'],
                genome = settings['genome'],
                genome_len_file = settings['genome']+'.fai',
                crispresso_cutoff = settings['crispresso_cutoff'],
                r1_read_ids_for_crispresso = r1_read_ids_for_crispresso,
                aligned_locs = genome_aligned_locs,
                av_read_length = av_read_length,
                genome_mapped_bam_file = genome_mapped_bam_file,
                run_crispresso_genome_sites = settings['run_crispresso_genome_sites'],
                crispresso_min_aln_score = settings['crispresso_min_aln_score'],
                samtools_command=settings['samtools_command'],
                crispresso_command=settings['crispresso_command'],
                query_bp_around_cut=10
                )
    crispresso_commands.extend(genome_crispresso_commands)
    crispresso_infos.extend(genome_crispresso_infos)

    if len(target_names) > 0:
        (targets_crispresso_infos, targets_crispresso_commands
            ) = prep_crispresso2_artificial_targets(
                    root = settings['root'],
                    genome_len_file=settings['genome']+'.fai',
                    crispresso_cutoff=settings['crispresso_cutoff'],
                    custom_index_fasta = custom_index_fasta,
                    custom_mapped_bam_file = custom_mapped_bam_file,
                    target_names = target_names,
                    target_info = target_info,
                    genome_read_count_at_cuts = genome_read_count_at_cuts,
                    crispresso_min_aln_score = settings['crispresso_min_aln_score'],
                    samtools_command=settings['samtools_command'],
                    crispresso_command=settings['crispresso_command'],
                    )
        crispresso_commands.extend(targets_crispresso_commands)
        crispresso_infos.extend(targets_crispresso_infos)



    crispresso_results = run_and_aggregate_crispresso(
                root = settings['root'],
                crispresso_infos = crispresso_infos,
                crispresso_commands = crispresso_commands
                )

    labels = ["Input Reads","Aligned Templates","Aligned Genome","Chopped Translocations","Chopped Large Deletions","Chopped Unidentified"]
    values = [num_reads_input,custom_aligned_count,genome_aligned_count,translocation_count,large_deletion_count,unidentified_count]
    alignment_summary_root = settings['root']+".alignmentSummary"
    with open(alignment_summary_root+".txt",'w') as summary:
        summary.write("\t".join(labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    ax.pie(values[1:],labels=[labels[idx]+"\n("+str(values[idx])+")" for idx in range(1,len(labels))],autopct="%1.2f%%")
    plot_name = alignment_summary_root
    plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')


    summary_plot_objects.append(
            PlotObject(plot_name = plot_name,
                    plot_title = 'Alignment Summary',
                    plot_label = 'Pie chart showing distribution of reads. Total reads in input: ' + str(num_reads_input),
                    plot_datas = [('Alignment summary',alignment_summary_root + ".txt")]
                    ))

    summary_plot_objects.append(global_aln_plot_obj)
    summary_plot_objects.append(frags_plot_obj)

    make_report(report_file=settings['root']+".html",
            report_name = 'Report',
            crisprlungo_folder = '',
            crispresso_run_names = crispresso_results['run_names'],
            crispresso_sub_html_files = crispresso_results['run_sub_htmls'],
            summary_plot_objects = summary_plot_objects,
            )


    cleanup(settings['root'])

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
    if len(args) < 2:
        raise Exception('Error: Settings file is missing\nUsage: ' +
                    args[0] + ' {settings file}')
    settings_file = args[1]
    settings = {}

    logging_level = logging.INFO
    if len(args) > 2 and 'debug' in args[2].lower():
        logging_level=logging.DEBUG
        settings['debug'] = True

    log_formatter = logging.Formatter("%(asctime)s:%(levelname)s: %(message)s")
    logging.basicConfig(
            level=logging_level,
            format="%(asctime)s:%(levelname)s: %(message)s",
            filename=settings_file+".CRISPRlungo.log",
            filemode='w'
            )
    ch = logging.StreamHandler()
    ch.setFormatter(log_formatter)
    logging.getLogger().addHandler(ch)

    logging.info('Parsing settings file')

    with open(settings_file, 'r') as fin:
        for line in fin:
            line_els = line.split("#")[0].strip().split("\t")
            if (len(line_els) < 2):
                raise Exception('Cannot parse line ' + line + '\nA tab must separate the key and value in the settings file')
            key = line_els[0].strip()
            val = line_els[1].strip()
            settings[key] = val

    settings['root'] = settings_file + ".CRISPRlungo"

    settings['keep_intermediate'] = (settings['keep_intermediate'] == 'True') if 'keep_intermediate' in settings else False

    # two parameters control the size and stride of fragment creation
    # size of fragment
    settings['fragment_size'] = int(settings['fragment_size']) if 'fragment_size' in settings else 20
    # step size between fragments
    settings['fragment_step_size'] = int(settings['fragment_step_size']) if 'fragment_step_size' in settings else 10

    # number of bp to extend beyond av read length around cut site for custom index
    settings['alignment_extension'] = int(settings['alignment_extension']) if 'alignment_extension' in settings else 50

    #min number of reads to analyze using crispresso
    settings['crispresso_cutoff'] = settings['crispresso_cutoff'] if 'crispresso_cutoff' in settings else 50

    #space-delimited list of cut sites in the form chr1:234 chr2:234
    settings['cuts'] = settings['cut_sites'].split(" ") if 'cut_sites' in settings else []

    #for finding offtargets with casoffinder
    settings['PAM'] = settings['PAM'] if 'PAM' in settings else None
    settings['on-target_guides'] = settings['on-target_guides'].split(" ") if 'on-target_guides' in settings else None
    settings['casoffinder_num_mismatches'] = int(settings['casoffinder_num_mismatches']) if 'casoffinder_num_mismatches' in settings else 5
    settings['cleavage_offset'] = int(settings['cleavage_offset']) if 'cleavage_offset' in settings else -3

    if 'add_non_primer_cut_targets' in settings and settings['add_non_primer_cut_targets'] == 'True':
        settings['add_non_primer_cut_targets'] = True
    elif 'primer_site' not in settings or settings['primer_site'] is not None:
        settings['add_non_primer_cut_targets'] = True
    elif 'primer_seq' not in settings or settings['primer_seq'] is not None:
        settings['add_non_primer_cut_targets'] = True
    else:
        settings['add_non_primer_cut_targets'] = False


    #boolean for whether to run crispresso on highly-aligned genomic locations (if false, crispresso will only be run at cut sites and artificial targets)
    if 'run_crispresso_genome_sites' in settings and settings['run_crispresso_genome_sites'] == 'True':
        settings['run_crispresso_genome_sites'] = True
    else:
        settings['run_crispresso_genome_sites'] = False

    settings['samtools_command'] = settings['samtools_command'] if 'samtools_command' in settings else 'samtools'
    settings['bowtie2_command'] = settings['bowtie2_command'] if 'bowtie2_command' in settings else 'bowtie2'
    settings['bowtie2_threads'] = int(settings['bowtie2_threads']) if 'bowtie2_threads' in settings else 1
    settings['crispresso_command'] = settings['crispresso_command'] if 'crispresso_command' in settings else 'CRISPResso'
    settings['casoffinder_command'] = settings['casoffinder_command'] if 'casoffinder_command' in settings else 'cas-offinder'

    settings['umi_regex'] = settings['umi_regex'] if 'umi_regex' in settings else 'NNWNNWNNN'
    settings['dedup_input'] = (settings['dedup_input'] == 'True') if 'dedup_input' in settings else False


    settings['crispresso_min_aln_score'] = int(settings['crispresso_min_aln_score']) if 'crispresso_min_aln_score' in settings else 20

    if 'genome' not in settings:
        raise Exception('Error: genome is required in settings file, pointing to a .fa file. The accompanying .fai file should also exist.')
    if not os.path.isfile(settings['genome']):
        raise Exception('Error: The genome file %s does not exist'%settings['genome'])
    genome_len_file=settings['genome']+'.fai'
    if not os.path.isfile(genome_len_file):
        raise Exception('Error: The genome length file %s does not exist'%genome_len_file)

    if 'bowtie2_genome' not in settings:
        potential_bowtie2_path = re.sub('.fa$','',settings['genome'])
        if os.path.isfile(potential_bowtie2_path+'.1.bt2'):
            settings['bowtie2_genome']= potential_bowtie2_path
        else:
            raise Exception('Error: bowtie2_genome is required in settings file, pointing to a bowtie2 genome (minus trailing .X.bt2)\nAlternately, set genome to the .fa file in a bowtie2 directory.')

    if 'fastq' in settings:
        settings['fastq_r1'] = settings['fastq']

    if 'fastq_r1' not in settings:
        raise Exception('Error: fastq_r1 file is required in settings file')

    if not os.path.isfile(settings['fastq_r1']):
        raise Exception('Error: fastq_r1 file %s does not exist',settings['fastq_r1'])

    if 'fastq_r2' in settings and settings['fastq_r2'].rstrip() != "":
        if not os.path.isfile(settings['fastq_r2']):
            raise Exception('Error: fastq_r2 file %s does not exist',settings['fastq_r2'])
    else:
        settings['fastq_r2'] = None

    if 'fastq_umi' in settings and settings['fastq_umi'].rstrip() != "":
        if not os.path.isfile(settings['fastq_umi']):
            raise Exception('Error: fastq_umi file %s does not exist',settings['fastq_umi'])
    else:
        settings['fastq_umi'] = None

    with open (settings['root']+".settingsUsed.txt",'w') as fout:
        for setting in settings:
            fout.write("%s: %s\n"%(str(setting),str(settings[setting])))
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
    except Exception:
        raise Exception('Error: CRISPResso2 is required')

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
    av_read_len = float(subprocess.check_output(cmd,shell=True).decode('utf-8').strip())

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

    res = int(subprocess.check_output(cmd,shell=True).decode('utf-8').split(" ")[0])

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

def make_artificial_targets(cuts,genome,target_length,target_padding,primer_chr ="",primer_loc=-1,primer_seq=None,add_non_primer_cut_targets=False,samtools_command='samtools'):
    """
    Generates fasta sequences surrounding cuts for alignment
    At each cut point, sequence of length target_length is generated
    Combinations of sequences at each cut site are produced
    If a primer is given, only cut sites that include the primer are used

    params:
        cuts: array of cut locations
        genome: location of fasta genome
        target_length: how long the query fragment should be (on one side of the cut) (not including padding). If the primer_seq is given, targets with primer_seq may be sorter than read_length.
        taget_padding: sequence (bp) padding around target (no padding for primer_seq).
        primer_chr: chromosome of primer (if any)
        primer_loc: location of primer (if any)
        primer_seq: sequence of primer binding (if any) (e.g. for dsODN integration and priming). Normally either primer_seq or (primer_chr and primer_loc) are given
        add_non_primer_cut_targets: boolean for whether to add targets for cuts without a primer. If false, only primer-cut1 pairings will be added. If true, cut1-cut2 pairings will be added.
        samtools_command: location of samtools to run

    returns:
        target_names: array of target names (corresponding to targets)
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targets)
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
            target_info[target_name]['query_pos']: genomic start of query (bp)
            target_info[target_name]['query_start']: bp after which query starts (in case of padding)
            target_info[target_name]['query_end']
            target_info[target_name]['target_cut_idx']: bp of cut in target
    """
    logging.info('Making artificial targets')
    target_names = []
    target_info = {}
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
                'cut1_site':'Primer',
                'cut2_chr':'Primer',
                'cut2_site':'Primer',
                'query_pos':0,
                'query_start':0,
                'query_end':len(primer_seq)*2,
                'target_cut_idx':len(primer_seq)
                }


        LALAc = left_bit_A + reverse_complement(left_bit_A)
        target_name = 'CRISPRlungo_PPc' #primer primer complement
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': LALAc,
                'class': 'Primed',
                'cut1_chr':'Primer',
                'cut1_site':'Primer',
                'cut2_chr':'Primer',
                'cut2_site':'Primer',
                'query_pos':0,
                'query_start':0,
                'query_end':len(primer_seq)*2,
                'target_cut_idx':len(primer_seq)
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
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,cut_start_B,cut_start_B_stop),shell=True).decode('utf-8').strip()
            left_bit_B = fasta_cache[left_bit_B_key]

            right_bit_B_key = '%s %s %d %d'%(genome,chr_B,site_B,cut_end_B)
            if right_bit_B_key not in fasta_cache:
                fasta_cache[right_bit_B_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,site_B,cut_end_B),shell=True).decode('utf-8').strip()
            right_bit_B = fasta_cache[right_bit_B_key]

            LARB = left_bit_A + right_bit_B
            target_name = 'CRISPRlungo_P' + 'R' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LARB,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'query_pos':0,
                    'query_start':0,
                    'query_end':len(primer_seq)+target_length,
                    'target_cut_idx':len(primer_seq)
                    }

            LARBc = left_bit_A + complement(right_bit_B)
            target_name = 'CRISPRlungo_P' + 'Rc' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LARBc,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'query_pos':0,
                    'query_start':0,
                    'query_end':len(primer_seq)+target_length,
                    'target_cut_idx':len(primer_seq)
                    }

            LALB = left_bit_A + reverse(left_bit_B)
            target_name = 'CRISPRlungo_P' + 'L' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LALB,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'query_pos':0,
                    'query_start':0,
                    'query_end':len(primer_seq)+target_length,
                    'target_cut_idx':len(primer_seq)
                    }

            LALBc = left_bit_A + reverse_complement(left_bit_B)
            target_name = 'CRISPRlungo_P' + 'Lc' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LALBc,
                    'class': 'Primed',
                    'cut1_chr':'Primer',
                    'cut1_site':'Primer',
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'query_pos':0,
                    'query_start':0,
                    'query_end':len(primer_seq)+target_length,
                    'target_cut_idx':len(primer_seq)
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
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_A,cut_start_A,cut_start_A_stop),shell=True).decode('utf-8').strip()
            left_bit_A = fasta_cache[left_bit_A_key]

            right_bit_A_key = '%s %s %d %d'%(genome,chr_A,site_A,cut_end_A)
            if right_bit_A_key not in fasta_cache:
                fasta_cache[right_bit_A_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_A,site_A,cut_end_A),shell=True).decode('utf-8').strip()
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
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'query_pos':cut_start_A,
                        'query_start':target_padding,
                        'query_end':target_padding + target_length*2,
                        'target_cut_idx':target_padding + target_length
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
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'query_pos':cut_start_A,
                        'query_start':target_padding,
                        'query_end':target_padding + target_length*2,
                        'target_cut_idx':target_padding + target_length
                        }


                LALAc = left_bit_A + reverse_complement(left_bit_A)
                target_name = 'CRISPRlungo_L' + str(i) + 'Lc' + str(i)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': LALAc,
                        'class': 'Chimera',
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'query_pos':cut_start_A,
                        'query_start':target_padding,
                        'query_end':target_padding + target_length*2,
                        'target_cut_idx':target_padding + target_length
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
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'query_pos':cut_start_A,
                        'query_start':target_padding,
                        'query_end':target_padding + target_length*2,
                        'target_cut_idx':target_padding + target_length
                        }

                RARAc = right_bit_A + reverse_complement(right_bit_A)
                target_name = 'CRISPRlungo_R' + str(i) + 'Rc' + str(i)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': RARAc,
                        'class': 'Chimera',
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut2_chr':chr_A,
                        'cut2_site':site_A,
                        'query_pos':cut_start_A,
                        'query_start':target_padding,
                        'query_end':target_padding + target_length*2,
                        'target_cut_idx':target_padding + target_length
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
                        '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,cut_start_B,cut_start_B_stop),shell=True).decode('utf-8').strip()
                left_bit_B = fasta_cache[left_bit_B_key]

                right_bit_B_key = '%s %s %d %d'%(genome,chr_B,site_B,cut_end_B)
                if right_bit_B_key not in fasta_cache:
                    fasta_cache[right_bit_B_key] = subprocess.check_output(
                        '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,site_B,cut_end_B),shell=True).decode('utf-8').strip()
                right_bit_B = fasta_cache[right_bit_B_key]


                if primer_chr == "" or primer_is_in_left_bit_A or primer_is_in_right_bit_B:
                    LARB = left_bit_A + right_bit_B
                    target_name = 'CRISPRlungo_L' + str(i) + 'R' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LARB,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
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
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
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
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
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
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
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
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                    LALBc = left_bit_A + reverse_complement(left_bit_B)
                    target_name = 'CRISPRlungo_L' + str(i) + 'Lc' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': LALBc,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                if primer_chr == "" or primer_is_in_right_bit_A or primer_is_in_right_bit_B:
                    RARB = right_bit_A + reverse(right_bit_B)
                    target_name = 'CRISPRlungo_R' + str(i) + 'R' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': RARB,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'

                    RARBc = right_bit_A + reverse_complement(right_bit_B)
                    target_name = 'CRISPRlungo_R' + str(i) + 'Rc' + str(j)
                    target_names.append(target_name)
                    target_info[target_name] = {
                            'sequence': RARBc,
                            'cut1_chr':chr_A,
                            'cut1_site':site_A,
                            'cut2_chr':chr_B,
                            'cut2_site':site_B,
                            'query_pos':cut_start_A,
                            'query_start':target_padding,
                            'query_end':target_padding + target_length*2,
                            'target_cut_idx':target_padding + target_length
                            }
                    if chr_A == chr_B:
                        target_info[target_name]['class'] = 'Large deletion'
                    else:
                        target_info[target_name]['class'] = 'Translocation'


    logging.info('Created ' + str(len(target_names)) + ' targets')
    return(target_names,target_info)


def make_target_index(root,target_names,target_info,bowtie2_command,bowtie2_threads=1):
    """
    Creates a bowtie index for artificial targets
    params:
        root: root for written files
        target_names: array of target names
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targets)
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
            target_info[target_name]['query_pos']: genomic start of query (bp)
            target_info[target_name]['query_start']: bp after which query starts (in case of padding)
            target_info[target_name]['query_end']
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
    returns:
        custom_index_fasta: fasta of artificial targets
    """
    custom_index_fasta = root + '.customIndex.fa'
    logging.info('Printing ' + str(len(target_names)) + ' targets to custom index (' + custom_index_fasta + ')')
    with open(custom_index_fasta,'w') as fout:
        for i in range(len(target_names)):
            fout.write('>'+target_names[i]+'\n'+target_info[target_names[i]]['sequence']+'\n')

    if not os.path.isfile(custom_index_fasta + '.1.bt2'):
        logging.info('Indexing custom targets using ' + bowtie2_command + '-build (' + custom_index_fasta + ')')
        index_result = subprocess.check_output(bowtie2_command + '-build --offrate 3 --threads ' + str(bowtie2_threads) + ' ' + custom_index_fasta + ' ' + custom_index_fasta, shell=True,stderr=subprocess.STDOUT)

    return custom_index_fasta

def dedup_file(root,fastq_r1,fastq_r2,umi_regex,min_umi_seen_to_keep_read=0,write_UMI_counts=False):
    """
    Deduplicates fastq files based on UMI (UMI is assumed to be the last part of the read ID after the ':'

    params:
        root: root for written files
        fastq_r1: R1 reads to dedup
        fastq_r2: R2 reads to dedup
        umi_regex: string specifying regex that UMI must match
        min_umi_seen_to_keep_read: min number of times a umi must be seen to keep the umi and read (e.g. if set to 2, a read-UMI pair that is only seen once will be discarded)
        write_UMI_counts: if True, writes a file with the UMI counts
    returns:
        fastq_r1_dedup: fastq_r1 file of deduplicated reads
        fastq_r2_dedup: fastq_r2 file of deduplicated reads
        tot_read_count: number of total reads read
        count_with_regex: number of reads with umi matching regex
        post_dedup_count: number of reads post deduplication
        post_dedup_read_count: number of reads that contributed to the deduplicated reads (post_dedup_read_count reads were seen that were deduplicated to post_dedup_count reads. If post_dedup_count == post_dedup_read_count then every UMI-read pair was seen once.)
    """

    dedup_stats_file = root + ".info"
    fastq_r1_dedup = "NA"
    fastq_r2_dedup = "NA"
    tot_read_count = -1
    count_with_regex = -1
    post_dedup_count = -1
    post_dedup_read_count = -1
    #if already finished, attempt to read in stats
    if os.path.isfile(dedup_stats_file):
        with open(dedup_stats_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().strip().split("\t")
            if len(line_els) > 3:
                (fastq_r1_dedup,fastq_r2_dedup,tot_read_count_str,count_with_regex_str,post_dedup_count_str,post_dedup_read_count_str) = line_els
                tot_read_count = int(tot_read_count_str)
                count_with_regex = int(count_with_regex_str)
                post_dedup_count = int(post_dedup_count_str)
                post_dedup_read_count = int(post_dedup_read_count_str)
                logging.info('Using previously-deduplicated fastqs with %d reads'%post_dedup_count)

    #otherwise perform deduplication (check if we were able to read in stats as well -- if we couldn't read them in, tot_read_count will be -1
    if tot_read_count == -1:
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
            f1_in = io.BufferedReader(gzip.open(fastq_r1,'rb'))
        else:
            f1_in = open(fastq_r1,'r')

        if fastq_r2.endswith('.gz'):
            f2_in = io.BufferedReader(gzip.open(fastq_r2,'rb'))
        else:
            f2_in = open(fastq_r2,'r')

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
        f1_out = gzip.open(fastq_r1_dedup, 'wb')
        fastq_r2_dedup = root + '.r2.gz'
        f2_out = gzip.open(fastq_r2_dedup, 'wb')

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
            loggin.info('Wrote UMI counts to ' + root+".umiCounts.txt")


        logging.info('Wrote %d deduplicated reads'%post_dedup_count)
        with open(dedup_stats_file,'w') as fout:
            fout.write("\t".join(["fastq_r1_dedup","fastq_r2_dedup","tot_read_count","count_with_regex","post_dedup_count","post_dedup_read_count"])+"\n")
            fout.write("\t".join([str(x) for x in [fastq_r1_dedup,fastq_r2_dedup,tot_read_count,count_with_regex,post_dedup_count,post_dedup_read_count]])+"\n")

    #done processing just plotting now
    return(fastq_r1_dedup,fastq_r2_dedup)

def add_umi_from_umi_file(root,fastq_r1,fastq_r2,fastq_umi):
    """
    Adds the UMI to the read ID from a UMI file (a third file with the UMI sequences per read)

    params:
        root: root for written files
        fastq_r1: R1 reads to dedup
        fastq_r2: R2 reads to dedup
        fastq_umi: UMI fastq to dedup on
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
    if os.path.isfile(umi_stats_file):
        with open(umi_stats_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().strip().split("\t")
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

def align_reads(root,fastq_r1,fastq_r2,bowtie2_reference,reference_name,target_info,bowtie2_command='bowtie2',bowtie2_threads=1,samtools_command='samtools',keep_intermediate=False):
    """
    Aligns reads to the provided reference (either artificial targets or genome)

    params:
        root: root for written files
        fastq_r1: fastq_r1 to align
        fastq_r2: fastq_r2 to align
        bowtie2_reference: bowtie2 reference to align to (either artificial targets or reference)
        reference_name: Name displayed to user for updates (e.g. 'Genome' or 'Artificial Targets')
        target_info: hash of information for each target_name
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        samtools_command: location of samtools to run
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)

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

    logging.info('Aligning to %s'%reference_name)

    mapped_bam_file = root + ".bam"

    bowtie_log = root + '.bowtie2Log'
    if fastq_r2 is not None: #paired-end reads
        logging.info('Aligning paired reads to %s using %s'%(reference_name,bowtie2_command))
        aln_command = '{bowtie2_command} --sam-no-qname-trunc --end-to-end --threads {bowtie2_threads} -x {bowtie2_reference} -1 {fastq_r1} -2 {fastq_r2} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logging.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT)
        logging.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning paired reads to %s\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(reference_name,aln_command,aln_result))
    #unpaired reads
    else:
        logging.info('Aligning single-end reads to %s using %s'%(reference_name,bowtie2_command))
        aln_command = '{bowtie2_command} --sam-no-qname-trunc --end-to-end --threads {bowtie2_threads} -x {bowtie2_reference} -U {fastq_r1} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logging.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT)
        logging.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning single-end reads to %s\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(reference_name,aln_command,aln_result))

    logging.info('Analyzing reads aligned to %s'%(reference_name))

    #analyze alignments
    # - if read aligned, store read id and loc in dict r1_assignments
    # - if read unaligned, write read to unaligned fastq
    mapped_tlens = defaultdict(int) # observed fragment lengths
    mapped_chrs = {} #stores the counts aligning to chrs or to templates
    aligned_locs = {} # aligned reads will not be chopped, but keep track of where they aligned for CRISPResso output
    aligned_chr_counts = {}

    unmapped_fastq_r1_file = root + '.unmapped.fastq'
    unmapped_fastq_r2_file = None
    r1_assignments_file = root + '.assignments.txt'
    r2_assignments_file = None
    if fastq_r2 is not None: #paired-end reads
        r1_assignments_file = root + '.assignments_r1.txt'
        unmapped_fastq_r1_file = root + '.unmapped_r1.fastq'
        r2_assignments_file = root + '.assignments_r2.txt'
        unmapped_fastq_r2_file = root + '.unmapped_r2.fastq'
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

        line_unmapped = int(line_els[1]) & 0x4
        if line_unmapped:
            seq = line_els[9]
            qual = line_els[10]
            if read_is_r1:
                uf1.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                unmapped_r1_count += 1
            else:
                uf2.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                unmapped_r2_count += 1
            continue


        line_chr = line_els[2]
        if not line_chr in mapped_chrs:
            mapped_chrs[line_chr] = 0
        mapped_chrs[line_chr] += 1

        line_mapq = line_els[5]
        line_start = int(line_els[3])-1

        if line_chr not in aligned_locs:
            aligned_locs[line_chr] = {}
        if line_start not in aligned_locs[line_chr]:
            aligned_locs[line_chr][line_start] = 0

        aligned_locs[line_chr][line_start] += 1
        if line_chr not in aligned_chr_counts:
            aligned_chr_counts[line_chr] = 0
        aligned_chr_counts[line_chr] += 1

        curr_classification = 'Linear'
        curr_position = line_chr + ":" + str(line_start)
        curr_annotation = ''
        curr_cut = line_chr + ":" + str(line_start) + "_" + str(line_start + len(line_els[9]))
        if line_chr in target_info: #if this aligned to a custom chromosome
            curr_classification = target_info[line_chr]['class']
            curr_annotation = line_chr
            custom_start = target_info[line_chr]['query_pos'] + line_start
            curr_position = target_info[line_chr]['cut1_chr'] +":" + str(custom_start)
            #if custom ref is not linear, set cuts and genome aln locations
            if curr_classification != 'Linear':

                cut1 =  target_info[line_chr]['cut1_chr'] + ":" +  str(target_info[line_chr]['cut1_site'])
                cut2 =  target_info[line_chr]['cut2_chr'] + ":" +  str(target_info[line_chr]['cut2_site'])
                curr_cut = cut1 + "-" + cut2
                if cut1 == cut2: #cut1 could == cut2 if it's a chimeric read
                    curr_cut = cut1
                else:
                    curr_position += "-"+ target_info[line_chr]['cut2_chr'] + ":" +  str(target_info[line_chr]['cut2_site'])

        if read_is_r1:
            af1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(read_id,reference_name,curr_classification,curr_annotation,curr_position,curr_cut))
        else:
            af2.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(read_id,reference_name,curr_classification,curr_annotation,curr_position,curr_cut))

        if read_is_paired_read_1: #only set if paired
            insert_size = int(line_els[8])
            mapped_tlens[insert_size] += 1

    uf1.close()
    af1.close()
    if unmapped_fastq_r2_file is not None:
        uf2.close()
        af2.close()

    chr_aln_plot_root = root + ".chrs"
    keys = sorted(aligned_chr_counts.keys())
    vals = [str(aligned_chr_counts[key]) for key in keys]
    with open(chr_aln_plot_root+".txt","w") as chrs:
        chrs.write('chr\tnumReads\n')
        for key in sorted(aligned_chr_counts.keys()):
            chrs.write(key + '\t' + str(aligned_chr_counts[key]) + '\n')

    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    ax.bar(range(len(keys)),vals,tick_label=keys)
    ax.set_ylabel('Number of Reads')
    ax.set_title('Location of reads aligned to the '+reference_name)
    plt.savefig(chr_aln_plot_root+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(chr_aln_plot_root+".png",pad_inches=1,bbox_inches='tight')

    chr_aln_plot_obj = PlotObject(
            plot_name = chr_aln_plot_root,
            plot_title = reference_name.capitalize() + ' Alignment Summary',
            plot_label = 'Bar plot showing alignment location of reads aligned to ' + reference_name,
            plot_datas = [(reference_name.capitalize() + 'Alignment Summary',chr_aln_plot_root + ".txt")]
            )

    tlen_plot_root = root + ".insertSizes"
    keys = sorted(mapped_tlens.keys())
    vals = [mapped_tlens[key] for key in keys]
    with open(tlen_plot_root+".txt","w") as fout:
        fout.write('insertSize\tnumReads\n')
        for key in keys:
            fout.write(str(key) + '\t' + str(mapped_tlens[key]) + '\n')

    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    ax.bar(keys,vals)
    ax.set_ylabel('Number of Reads')
    ax.set_title('Insert size')
    plt.savefig(tlen_plot_root+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(tlen_plot_root+".png",pad_inches=1,bbox_inches='tight')

    tlen_plot_obj = PlotObject(
            plot_name = tlen_plot_root,
            plot_title = reference_name.capitalize() + ' Alignment Insert Size Summary',
            plot_label = 'Bar plot showing insert size of reads aligned to ' + reference_name,
            plot_datas = [(reference_name.capitalize() + ' Alignment Insert Size Summary',tlen_plot_root + ".txt")]
            )

    aligned_count = read_count - (unmapped_r1_count + unmapped_r2_count)

    with open (root+".info",'w') as fout:
        fout.write("\t".join([str(x) for x in ["read_count","r1_count","r2_count","unmapped_r1_count","unmapped_r2_count","aligned_count","unmapped_fastq_r1_file","unmapped_fastq_r2_file"]])+"\n")
        fout.write("\t".join([str(x) for x in [read_count,r1_count,r2_count,unmapped_r1_count,unmapped_r2_count,aligned_count,unmapped_fastq_r1_file,unmapped_fastq_r2_file]]))
    return(r1_assignments_file,r2_assignments_file,unmapped_fastq_r1_file,unmapped_fastq_r2_file,aligned_count,mapped_bam_file,chr_aln_plot_obj,tlen_plot_obj)

def make_final_read_assignments(root,r1_assignment_files,r2_assignment_files,cut_sites,target_info,cut_merge_dist=200,genome_map_resolution=1000000):
    """
    Makes final read assignments after:
        Pairing R1 and R2 if they were processed in separate steps
        Deduplicating based on UMI and alignment location

    params:
        root: root for written files
        r1_assignment_files: list of r1 assignment files
        r2_assignment_files: list of r2 assignment files
        cut_sites: array of known cut sites
        target_info: hash of information for each target_name
        cut_merge_dist: cuts within this distance (bp) will be merged
        genome_map_resolution: window size (bp) for reporting number of reads aligned

    returns:
        final_assignment_filename: filename of final assignments
        r1_read_ids_for_crispresso: dict of cut=>list of read ids for each cut
    """

    logging.info('Making final read assignments')

    sorted_r1_file = root + '.r1.sorted'
    sorted_r2_file = root + '.r2.sorted'
    final_file = root + '.final_assignments'
    final_file_tmp = root + '.final_assignments.tmp'
    cat_and_sort_cmd = 'cat ' + " ".join(r1_assignment_files) + ' | sort > ' + sorted_r1_file
    logging.debug('r1 cat command: ' + str(cat_and_sort_cmd))
    cat_and_sort_result = subprocess.check_output(cat_and_sort_cmd, shell=True,stderr=subprocess.STDOUT)

    cat_and_sort_cmd = 'cat ' + " ".join(r2_assignment_files) + ' | sort > ' + sorted_r2_file
    logging.debug('r2 cat command: ' + str(cat_and_sort_cmd))
    cat_and_sort_result = subprocess.check_output(cat_and_sort_cmd, shell=True,stderr=subprocess.STDOUT)

    #first pass through file:
    # -merge r1 and r2
    # -deduplicate by UMI
    # -aggregate cut sites that will be later cleaned

    id_ind = 0
    source_ind = 1
    classification_ind = 2
    annotation_ind = 3
    alignment_ind = 4
    cut_point_ind = 5

    seen_reads = {} #put into this dict if it's been seen for deduplicating
    total_reads_processed = 0
    dups_seen = 0
    r1_seen_custom_targets = defaultdict(int) #count of each custom target
    cut_points_by_chr = {} # dict of arrays for each chr contianing observedcuts
    with open(sorted_r1_file,'r') as sf1, open(sorted_r2_file,'r') as sf2, open(final_file_tmp,'w') as ff:
        for r1_line in sf1:
            r1_line = r1_line.strip()
            r2_line = sf2.readline().strip()
            total_reads_processed += 1

            r1_line_els = r1_line.split("\t")
            r1_id = r1_line_els[0]
            r1_umi = r1_line_els[0].split(":")[-1]

            r2_line_els = r2_line.split("\t")
            r2_id = r2_line_els[0]
            r2_umi = r2_line_els[0].split(":")[-1]


            if r1_id.split(" ")[0] != r2_id.split(" ")[0]:
                raise Exception("Some reads were lost -- read ids don't match up:\n1:"+r1_line+"\n2:"+r2_line)

            this_key = r1_umi + " " + r1_line_els[alignment_ind] + " " + r2_line_els[alignment_ind]
            #if this is a duplicate, skip it
            if this_key in seen_reads:
                dups_seen += 1
                ff.write(r1_line + "\t" + r2_line + "\tduplicate\t"+seen_reads[this_key]+"\n")
                continue

            #if not a duplicate:
            seen_reads[this_key] = r1_id
            ff.write(r1_line + "\t" + r2_line + "\tnot_duplicate\tNA\n")

            r1_seen_custom_targets[r1_line_els[annotation_ind]] += 1

            #parse out cut positions and record them
            r1_cut_points = r1_line_els[cut_point_ind]
            r2_cut_points = r2_line_els[cut_point_ind]
            if r1_line_els[classification_ind] != 'Linear':
                for cut_point in r1_cut_points.split("-"):
                    (cut_chr,cut_pos) = cut_point.split(":")
                    if cut_chr != "*":
                        if cut_chr not in cut_points_by_chr:
                            cut_points_by_chr[cut_chr] = []
                        cut_points_by_chr[cut_chr].append(int(cut_pos))
            if r2_line_els[classification_ind] != 'Linear':
                for cut_point in r2_cut_points.split("-"):
                    (cut_chr,cut_pos) = cut_point.split(":")
                    if cut_chr != "*":
                        if cut_chr not in cut_points_by_chr:
                            cut_points_by_chr[cut_chr] = []
                        cut_points_by_chr[cut_chr].append(int(cut_pos))

    #now that we've aggregated all cut sites, find likely cut positions and record assignments to those positions in final_cut_point_lookup
    #cut sites are input by user
    known_cut_points_by_chr = {}
    for cut_site in cut_sites:
        (cut_chr,cut_pos) = cut_site.split(":")
        if cut_chr not in known_cut_points_by_chr:
            known_cut_points_by_chr[cut_chr] = []
        known_cut_points_by_chr[cut_chr].append(int(cut_pos))

    with open(root+'.seen_cut_points','w') as fout:
        for cut_chr in cut_points_by_chr:
            these_cut_points = sorted(cut_points_by_chr[cut_chr])
            for cut_point in these_cut_points:
                fout.write(cut_chr + ':' + str(cut_point)+"\n")


    final_cut_points_by_chr = {} #dict of final cut points by chromozome
    final_cut_point_lookup = {} #dict from old(fuzzy/imprecise) position to new

    observed_cut_points = defaultdict(int) # hash of cut point -> count
    for cut_chr in cut_points_by_chr:
        these_cut_points = sorted(cut_points_by_chr[cut_chr])
        last_seen_point = these_cut_points[0]
        curr_points = [last_seen_point]
        for i in range(1,len(these_cut_points)+1):
            #if this point is within cut_merge_dist, merge!
            if i == len(these_cut_points) or abs(these_cut_points[i] - these_cut_points[i-1]) > cut_merge_dist:
                this_mean = int(sum(curr_points)/float(len(curr_points)))
                this_pos = this_mean
                min_dist = cut_merge_dist
                if cut_chr in known_cut_points_by_chr:
                    for known_cut_point in known_cut_points_by_chr[cut_chr]:
                        this_dist = abs(known_cut_point - this_mean)
                        if this_dist <= min_dist:
                            this_pos = known_cut_point
                            min_dist = this_dist
                observed_cut_points[cut_chr+":"+str(this_pos)] += len(these_cut_points)
                if cut_chr not in final_cut_points_by_chr:
                    final_cut_points_by_chr[cut_chr] = []
                final_cut_points_by_chr[cut_chr].append(this_pos)

                #add this assignment to the lookup
                for curr_point in curr_points:
                    final_cut_point_lookup[cut_chr+":"+str(curr_point)] = cut_chr+":"+str(this_pos)

                curr_points = []

            if i < len(these_cut_points):
                curr_points.append(these_cut_points[i])

    with open(root + ".final_cut_points.txt",'w') as fout:
        for chrom in sorted(final_cut_points_by_chr.keys()):
            for pos in sorted(final_cut_points_by_chr[chrom]):
                fout.write("%s\t%d\n"%(chrom,pos))
    with open(root+".final_cut_points_lookup.txt",'w') as fout:
        for key in sorted(final_cut_point_lookup.keys()):
            fout.write("%s\t%s\n"%(key,final_cut_point_lookup[key]))


    final_total_count = 0
    final_duplicate_count = 0


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

    #keep track of translocations
    r1_translocation_cut_counts = {}

    #keep track of read ids aligning to cuts
    #reads with these ids will be pulled out of the fastq for analysis by CRISPResso
    r1_read_ids_for_crispresso = {}

    #column indices in the final file
    ind_offset = cut_point_ind + 1
    r1_id_ind = id_ind
    r1_source_ind = source_ind
    r1_classification_ind = classification_ind
    r1_alignment_ind = alignment_ind
    r1_cut_point_ind = cut_point_ind
    r2_source_ind = source_ind + ind_offset
    r2_classification_ind = classification_ind + ind_offset
    r2_alignment_ind = alignment_ind + ind_offset
    r2_cut_point_ind = cut_point_ind + ind_offset
    duplicate_ind = 12
    with open(final_file_tmp,'r') as fin, open(final_file,'w') as fout:
        for line in fin:
            line = line.strip()
            final_total_count += 1
            line_els = line.split("\t")
            if line_els[duplicate_ind] == "duplicate":
                final_duplicate_count += 1
                continue
            r1_classification = line_els[r1_classification_ind]
            r1_cut_points = line_els[r1_cut_point_ind]
            r1_source = line_els[r1_source_ind]

            r1_associated_cut_points = []

            if r1_classification == 'Linear':
                (r1_chr,r1_range) = r1_cut_points.split(":")
                (r1_start,r1_stop) = r1_range.split("_")
                #check to see if this read overlaps with any known cutpoints
                if r1_chr in final_cut_points_by_chr:
                    for cut_point in final_cut_points_by_chr[r1_chr]:
                        if cut_point > int(r1_start) and cut_point < int(r1_stop):
                            key = r1_chr+":"+str(cut_point)
                            r1_associated_cut_points.append(key)
                            r1_uncut_counts[key] += 1
                            if r1_source == 'genome':
                                r1_uncut_genome_counts[key] += 1
                            elif r1_source == 'fragmented':
                                r1_uncut_frag_counts[key] += 1
                            elif r1_source == 'custom targets':
                                r1_uncut_custom_aln[key] += 1
            else:
                #if it's not Linear, we've already annotated the cut points
                for cut_point in r1_cut_points.split("-"):
                    if "*" not in cut_point:
                        final_cut_point = final_cut_point_lookup[cut_point]
                        r1_associated_cut_points.append(final_cut_point)
                        r1_cut_counts[final_cut_point] += 1
                        if r1_source == 'fragmented':
                            r1_cut_frag_counts[final_cut_point] += 1
                        elif r1_source == 'custom targets':
                            r1_cut_custom_aln_counts[final_cut_point] += 1

            #map translocation cuts
            if len(r1_associated_cut_points) > 0:
                cut1 = r1_associated_cut_points[0]
                cut2 = r1_associated_cut_points[-1] #this may be the same as [0] for linear ones
                if cut1 not in r1_translocation_cut_counts:
                    r1_translocation_cut_counts[cut1] = defaultdict(int)
                if cut2 not in r1_translocation_cut_counts:
                    r1_translocation_cut_counts[cut2] = defaultdict(int)
                r1_translocation_cut_counts[cut1][cut2] += 1
                if cut1 != cut2:
                    r1_translocation_cut_counts[cut2][cut1] += 1

                #also keep track of read ids covering each cut for analysis by CRISPResso
                cut_key = cut1 + " " + cut2
                if cut_key not in r1_read_ids_for_crispresso:
                    r1_read_ids_for_crispresso[cut_key] = []
                r1_read_ids_for_crispresso[cut_key].append(line_els[r1_id_ind])

            r1_position = line_els[r1_alignment_ind]
            for r1_aln_pos in r1_position.split("-"):
                (aln_chr,aln_pos) = r1_aln_pos.split(":")
                aln_pos_window = int(aln_pos)/genome_map_resolution
                if aln_chr not in r1_aln_pos_counts_by_chr:
                    r1_aln_pos_counts_by_chr[aln_chr] = defaultdict(int)
                    r1_aln_pos_genome_counts_by_chr[aln_chr] = defaultdict(int)
                    r1_aln_pos_custom_aln_counts_by_chr[aln_chr] = defaultdict(int)
                    r1_aln_pos_frag_counts_by_chr[aln_chr] = defaultdict(int)
                #assign each read to neighboring windows
                r1_aln_pos_counts_by_chr[aln_chr][aln_pos_window] += 1
                if r1_source == 'genome':
                    r1_aln_pos_genome_counts_by_chr[aln_chr][aln_pos_window] += 1
                elif r1_source == 'fragmented':
                    r1_aln_pos_frag_counts_by_chr[aln_chr][aln_pos_window] += 1
                elif r1_source == 'custom targets':
                    r1_aln_pos_custom_aln_counts_by_chr[aln_chr][aln_pos_window] += 1

            fout.write("%s\t%s"%(line,",".join(r1_associated_cut_points))+"\n")

    logging.info('Deduplicated %d/%d reads (kept %d)'%(final_duplicate_count,final_total_count,final_total_count-final_duplicate_count))


    cut_report = root + ".r1_cut_report.txt"
    with open (cut_report,'w') as fout:
        fout.write('chr\tstart\tcut_tot\tuncut_tot\tuncut_genome\tcut_custom\tuncut_custom\tcut_frags\tuncut_frags\n')
        for chrom in sorted(final_cut_points_by_chr.keys()):
            if chrom == "*":
                continue
            for pos in sorted(final_cut_points_by_chr[chrom]):
                cut_key = chrom + ":" + str(pos)
                known_str = 'Novel'
                if chrom in known_cut_points_by_chr and pos in known_cut_points_by_chr[chrom]:
                    known_str = 'Known'
                fout.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"%(chrom,pos,
                    r1_cut_counts[cut_key],
                    r1_uncut_counts[cut_key],
                    r1_uncut_genome_counts[cut_key],
                    r1_cut_custom_aln_counts[cut_key],
                    r1_uncut_custom_aln_counts[cut_key],
                    r1_cut_frag_counts[cut_key],
                    r1_uncut_frag_counts[cut_key],
                    known_str))

    logging.info('Wrote cut report ' + cut_report)
    tx_keys = sorted(r1_translocation_cut_counts.keys())
    tx_report = root + ".r1_translocation_table.txt"
    with open (tx_report,'w') as fout:
        fout.write('Translocations\t'+"\t".join(tx_keys)+"\n")
        for tx_key in tx_keys:
            fout.write(tx_key+"\t"+"\t".join([str(r1_translocation_cut_counts[tx_key][x]) for x in tx_keys])+"\n")
    logging.info('Wrote translocation table ' + tx_report)

    align_report = root + ".r1_align_report.txt"
    with open (align_report,'w') as fout:
        fout.write('chr\tstart\tread_tot\tread_genome\tread_custom\tread_frags\n')
        for chrom in sorted(r1_aln_pos_counts_by_chr.keys()):
            for pos in sorted(r1_aln_pos_counts_by_chr[chrom].keys()):

                #divide counts by two because we counted each read in two windows above
                fout.write("%s\t%d\t%d\t%d\t%d\t%d\n"%(chrom,pos*genome_map_resolution,
                    r1_aln_pos_counts_by_chr[chrom][pos],
                    r1_aln_pos_genome_counts_by_chr[chrom][pos],
                    r1_aln_pos_custom_aln_counts_by_chr[chrom][pos],
                    r1_aln_pos_frag_counts_by_chr[chrom][pos]))


    print('Finished final assignments')
    return final_file,r1_read_ids_for_crispresso




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

def chop_reads(root,unmapped_reads_fastq_r1,unmapped_reads_fastq_r2,bowtie2_genome,fragment_size,fragment_step_size,bowtie2_command='bowtie2',bowtie2_threads=1,samtools_command='samtools'):
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
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        samtools_command: location of samtools to run


    returns:
        r1_assignments_file: file of location assignments for each r1 read
        r2_assignments_file: file of location assignments for each r2 read
        translocation_count: number of reads that were identified as translocations
        large_deletion_count: number of reads that were identified as large deletion
        unidentified_count: number of reads that couldn't be identified as translocations
        frags_plot_obj: plot object summarizing fragments
    """
    logging.info('Creating read fragments')

    unmapped_ids = {}
    unmapped_id = 0
    max_frags = 0 #max frags produced for a read
    frags_per_read = {} #number of fragments created per read

    read_files_to_frag = [unmapped_reads_fastq_r1]
    if unmapped_reads_fastq_r2 is not None:
        read_files_to_frag.append(unmapped_reads_fastq_r2)
    unmapped_frag_file = root + ".unmapped.fastq"
    with open(unmapped_frag_file,"w") as unmapped_fastq:
        for (read_idx,unmapped_reads_fastq) in enumerate(read_files_to_frag):
            with open(unmapped_reads_fastq,'r') as unmapped_reads:
                line_id = unmapped_reads.readline()
                while(line_id):
                    line_id = line_id.strip()
                    line_seq = unmapped_reads.readline().strip()
                    line_plus = unmapped_reads.readline().strip()
                    line_qual = unmapped_reads.readline().strip()

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
                        unmapped_fastq.write("%s\n%s\n+\n%s\n"%(new_id,new_seq,new_qual))
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
    vals = [str(frags_per_read[key]) for key in keys]
    with open(frags_per_unaligned_read_root + ".txt","w") as frags:
        frags.write('numFragments\tnumReads\n')
        for key in keys:
            frags.write(str(key) + '\t' + str(frags_per_read[key]) + '\n')

    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    ax.bar(range(len(keys)),vals,tick_label=keys)
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


    mapped_chopped_sam_file = root + ".mapped.sam"
    chopped_bowtie_log = root + ".mapped.bowtie2Log"
    chopped_aln_command = '%s --sam-no-qname-trunc --reorder --no-hd --threads %d --end-to-end -x %s -U %s -S %s' %(bowtie2_command,bowtie2_threads,bowtie2_genome,unmapped_frag_file,mapped_chopped_sam_file)

    logging.debug(chopped_aln_command)
    aln_result = subprocess.check_output(chopped_aln_command,shell=True,stderr=subprocess.STDOUT).decode('utf-8')
    logging.debug('Alignment of chopped fragments to genome: ' + aln_result)
    with open (chopped_bowtie_log,'w') as lout:
        lout.write('Alignment of chopped fragments to genome\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(chopped_aln_command,aln_result))



    logging.info('Analyzing chopped reads')

    def analyze_curr_id_vals(curr_id_vals,fragment_size,fragment_step_size,max_frags,min_seen_count=2):
        """
        Helper function to analyze a list of curr_id_vals
        A midpoint is found that splits the aligned fragments such that the number of locations in each split is minimized. 
        A string describing this read and the aligned fragments is also created
        I pulled it out into a separate function because it's used in two spots and may need to be changed

        params:
            curr_id_vals: dict of fragNum:alnPos for each fragment in a specific read
            fragment_size: size of resulting chopped fragments.
            fragment_step_size: stride between chopped fragments. 5 results in a new fragment every 5bp
            max_frags: the max number of frags to print (for constructing the outline to be printed)
            min_seen_count: the minimum number of fragments aligned to a location for it to be counted (e.g. if a read only has 3 fragments aligning to location A and 1 fragment aligning to location B, a min_seen_count of 2 would not allow position A)
        returns:
            outline: String to be printed to output with locationa and max_frags number of columns describing where each fragment aligned
            curr_chr_count: how many chrs fragments from this read aligned to 
            left_chr: Chr location where the left part of this read was determined to come from
            left_cut_pos: Base Location where the left part could have resulted from a cut at
            left_pos: Base location where the left part of the read aligned
            right_chr: Location where the right part of this read was determined to come from
            right_cut_pos
            right_pos
            left_potential_cut: Location of genomic cut producing left part of this read
            right_potential_cut: Location of genomic cut producing right part of read
        """
        # first create arrays from the dict of positions
        curr_id_chrs = {} #dict of where fragments align chr->count
        pos_list = [None]*max_frags #chr and pos where each frag mapped
        val_list = [] #for printing the outline
        mapped_frag_len = 0 # length of array for which all frags are mapped frag (some reads are shorter, and don't have any fragments on the left side). I used the length so I can perform range and indexing operations (as opposed to keeping the index of the last mapped frag)

        for i in range(max_frags):
            val = ""
            if i in curr_id_vals:
                mapped_frag_len = i+1
                val = curr_id_vals[i]
                val_chr,val_pos = curr_id_vals[i].split(" ")
                if val_chr != "*":
                    if val_chr not in curr_id_chrs:
                        curr_id_chrs[val_chr] = 0
                    curr_id_chrs[val_chr] += 1
                    pos_list[i] = val

            val_list.append(val)

        min_pos_unique_sum = max_frags
        max_majority_sum = 0
        min_i = None
        #create partitions of the fragments, choose the one that:
        # - minimizes the number of unique locations seen in each partition
        # - maximizes the number of fragments that are the majority in their respective partitions
#        print('pos list: ' + str(pos_list))
#        print('last frag ind: ' + str(mapped_frag_len))
        for i in range(1,mapped_frag_len):
            left_list = filter(None,pos_list[0:i])
            num_pos_left = len(set(left_list))

            right_list = filter(None,pos_list[i:mapped_frag_len])
            len_right_list = len(right_list)
            num_pos_right = len(set(filter(None,pos_list[i:mapped_frag_len])))

            #if one of the partitions has no pos (all are None), continue
            if num_pos_right == 0 or num_pos_left == 0:
                continue

#            if val_list[0] == "chr1 223763905" and val_list[1] == "chr2 60498399":
#                print("i:"+str(i)+"/"+str(mapped_frag_len))
#                print('left bit: ' + str(pos_list[0:i]))
#                print('right bit: ' + str(pos_list[i:mapped_frag_len]))
#                print('num_pos_left: ' + str(num_pos_left))
#                print('num_pos_right: ' + str(num_pos_right))
#                print('ratio left: ' + str(ratio_left))
#                print('ratio right: ' + str(ratio_right))
#                if i > 5:
#                    toggle = 1
#                if toggle == 1 and i < 5:
#                    asdf()

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

#        print('min_i:' + str(min_i))
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
                if left_counts[left_val] >= min_seen_count:
                    left_chr,left_pos = left_val.split(" ")
                else: #an alignment was seen, but below the cutoff
                    left_chr = "*"
                    left_pos = "1"
            else: #no alignments were seen
                left_chr = "*"
                left_pos = "2"

            right_counts = {}
            right_val = None;
            right_count = 0
            for val in pos_list[min_i:mapped_frag_len]:
                right_counts[val] = right_counts.get(val,0) + 1
                if right_counts[val] > right_count and not val is None:
                    right_count = right_counts[val]
                    right_val = val
            if right_count > 0:
                if right_counts[right_val] >= min_seen_count:
                    right_chr,right_pos = right_val.split(" ")
                else: #an alignment was seen, but below the cutoff
                    right_chr = "*"
                    right_pos = "1"
            else: #no alignments were seen
                right_chr = "*"
                right_pos = "2"
        else: #no alignments were seen, min_i is None
            left_chr = "*"
            left_pos = "3"
            right_chr = "*"
            right_pos = "3"
            min_i = -1

        # estimate position where cuts may have occurred to create this read
        left_cut_pos = "*"
        if left_chr != "*":
            left_cut_pos = int(left_pos) + int((min_i*fragment_step_size)*fragment_size/2)
        right_cut_pos = "*"
        if right_chr != "*":
            right_cut_pos = int(right_pos) + int((min_i*fragment_step_size)*fragment_size/2)

        outline = '%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n'%(left_chr,str(left_cut_pos),left_pos,right_chr,str(right_cut_pos),right_pos,min_i,"\t".join(val_list))

        #sort so we don't count mix up chr1>chr2 with chr2>chr1 when they are really the same
        if right_chr < left_chr or (right_chr == left_chr and int(right_pos) < int(left_pos)):
            tmp_chr = left_chr
            tmp_pos = left_pos
            tmp_cut_pos = left_cut_pos
            left_chr = right_chr
            left_pos = right_pos
            left_cut_pos = right_cut_pos
            right_chr = tmp_chr
            right_pos = tmp_pos
            right_cut_pos = tmp_cut_pos

        curr_chr_count = len(curr_id_chrs.keys())
        return outline,curr_chr_count,left_chr,left_cut_pos,left_pos,right_chr,right_cut_pos,right_pos

    frag_meta_file = root + ".meta.txt"

    r1_assignments_file = root + '.assignments_r1.txt'
    r2_assignments_file = root + '.assignments_r2.txt'
    frags_mapped_chrs = defaultdict(int)
    with open(frag_meta_file,"w") as frag_file, open(r1_assignments_file,'w') as af1, open(r2_assignments_file,'w') as af2:
        head = "ID\tleft_chr\tleft_cut\tleft_pos\tright_chr\tright_cut\tright_pos\tbreak_index"
        for i in range(max_frags):
            head += "\t"+str(i)
        frag_file.write(head+"\n")

        with open(mapped_chopped_sam_file,'r') as aln_frags:
            curr_id = "" #whole id of the read (alignments will be aggregated based on this id)
            curr_idx = "" #read index (either 1 for r1 or 2 for r2
            curr_id_chrs = {} # keep track of which chrs frags from this id map to
            curr_id_vals = {} # keep track of where each frag maps to (this is a hash but it's like an array, so we can check for uninitialized values)
            translocations = defaultdict(int) # hash of all translocation locations
            translocation_count = 0
            large_deletion_count = 0
            unidentified_count = 0
            chroms_per_frag_read_count = defaultdict(int) # keep track for all read how many chrs the fragments mapped to
            frags_mapped_count = 0
            for line in aln_frags:
                line_els = line.split("\t")
                line_chr = line_els[2]
                frags_mapped_chrs[line_chr] += 1
                frags_mapped_count += 1

                line_mapq = line_els[5]
                line_unmapped = int(line_els[1]) & 0x4
                line_reverse = int(line_els[1]) & 0x10
                line_start = int(line_els[3])-1

                match = re.match('(.*)\.CLR(\d)\.CLID(\d+)\.CLO(-?\d+)\.CLF(-?\d+)$',line_els[0])
                if not match:
                    raise Exception('Cannot parse id %s from line %s in %s. Perhaps line was trimmed?\n'%(line_els[0],line,mapped_chopped_sam_file))
                (next_id,read_idx,lungo_id,lungo_offset,lungo_frag) = match.groups()

                #write results from last id if this line is for a new id
                if curr_id != next_id and curr_id != "":
                    (outline,curr_chr_count,first_chr,first_cut_pos,first_pos,last_chr,last_cut_pos,last_pos) = analyze_curr_id_vals(curr_id_vals,fragment_size,fragment_step_size,max_frags)
                    frag_file.write(curr_id+"\t"+outline)

                    curr_classification = "NA"
                    if first_chr == "*" or last_chr == "*":
                        unidentified_count += 1
                        curr_classification = "Unidentified"
                    else:
                        key = first_chr + " " + first_pos + " " + last_chr + " " + last_pos
                        translocations[key] += 1

                        #can't be linear otherwise it would have mapped in genomic alignment previuosly
                        if first_chr == last_chr:
                            large_deletion_count += 1
                            curr_classification = "Large deletion"
                        else:
                            translocation_count += 1
                            curr_classification = "Translocation"

                    chroms_per_frag_read_count[curr_chr_count] += 1
                    curr_loc = first_chr + ":" + str(first_pos) + "-" + last_chr + ":" + str(last_pos)
                    curr_cuts = first_chr + ":" + str(first_cut_pos) + "-" + last_chr + ":" + str(last_cut_pos)
                    curr_annotation = ''
                    if int(curr_idx) == 1:
                        af1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))
                    else:
                        af2.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))

                    curr_id_chrs = {}
                    curr_id_vals = {}

                curr_id = next_id
                curr_idx = read_idx

                inferred_start = line_start - int(lungo_offset)
                if line_reverse:
                    inferred_start = line_start + int(lungo_offset)
                curr_id_vals[int(lungo_frag)] = line_chr + " " + str(inferred_start)

            #finish last item
            if frags_mapped_count > 0:
                (outline,curr_chr_count,first_chr,first_cut,first_pos,last_chr,last_cut,last_pos) = analyze_curr_id_vals(curr_id_vals,fragment_size,fragment_step_size,max_frags)
                frag_file.write(curr_id+"\t"+outline)

                if first_chr == "*" or last_chr == "*":
                    unidentified_count += 1
                else:
                    key = first_chr + " " + first_pos + " " + last_chr + " " + last_pos
                    if key not in translocations:
                        translocations[key] = 0
                    translocations[key] += 1

                    if first_chr == last_chr:
                        large_deletion_count += 1
                    else:
                        translocation_count += 1

                chroms_per_frag_read_count[curr_chr_count] += 1
                curr_loc = first_chr + ":" + str(first_pos) + "-" + last_chr + ":" + str(last_pos)
                curr_cuts = first_chr + ":" + str(first_cut_pos) + "-" + last_chr + ":" + str(last_cut_pos)
                curr_annotation = ''
                if int(curr_idx) == 1:
                    af1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))
                else:
                    af2.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(curr_id,'fragmented',curr_classification,curr_annotation,curr_loc,curr_cuts))
            #done with last one

    logging.info("debug read frags mapped count: " + str(frags_mapped_count))
    logging.info("Found %d translocations, %d large deletions, and %d unidentified reads"%(translocation_count,large_deletion_count,unidentified_count))

    frags_aligned_chrs_root = root + ".chrs"
    keys = sorted(frags_mapped_chrs.keys())
    vals = [str(frags_mapped_chrs[key]) for key in keys]
    with open(frags_aligned_chrs_root+".txt","w") as fout:
        fout.write("chr\tnumReads\n")
        for key in keys:
            fout.write(str(key) + '\t' + str(frags_mapped_chrs[key]) + '\n')
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    ax.bar(range(len(keys)),vals,tick_label=keys)
    ax.set_ylabel('Number of Reads')
    ax.set_title('Location of fragments from unaligned reads')

    plt.savefig(frags_aligned_chrs_root+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(frags_aligned_chrs_root+".png",pad_inches=1,bbox_inches='tight')

    frag_chroms_per_read_root = root + ".alignedChromsPerRead"
    keys = sorted(chroms_per_frag_read_count.keys())
    vals = [str(chroms_per_frag_read_count[key]) for key in keys]
    with open(frag_chroms_per_read_root+".txt","w") as fout:
        fout.write("numChroms\tnumReads\n")
        for key in keys:
            fout.write(str(key) + '\t' + str(chroms_per_frag_read_count[key]) + '\n')

    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    ax.bar(range(len(keys)),vals,tick_label=keys)
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

    return (r1_assignments_file,r2_assignments_file,translocation_count,large_deletion_count,unidentified_count,frags_plot_obj)


def prep_crispresso2_global(root,cuts,genome,genome_len_file,crispresso_cutoff,aligned_locs,av_read_length,genome_mapped_bam_file,run_crispresso_genome_sites,crispresso_min_aln_score,samtools_command,crispresso_command,query_bp_around_cut):
    """
    Prepares globally-aligned data for crispresso2
    Frequently-aligned locations with a min number of reads (crispresso_cutoff) are identified
    Reads from these locations are extracted and prepared for analysis by CRISPResso2

    params:
        root: root for written files
        final_assignments_file: file with final assignments for each read
        cuts: array of cut locations -- these will be pulled out for crispresso analysis
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targets)
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
            target_info[target_name]['query_pos']: genomic start of query (bp)
            target_info[target_name]['query_start']: bp after which query starts (in case of padding)
            target_info[target_name]['query_end']
            target_info[target_name]['target_cut_idx']: bp of cut in target

        genome: path to genome fa file
        genome_len_file: path to tab-sep file with lengths of chrs in genome (ends with .fai)
        crispresso_cutoff: min number of reads at site for crispresso2 processing
        aligned_locs: hash of locations aligned to on each chromosome aligned_locs[chr][position] = count
        av_read_length: length of average sequencing read in sample
        genome_mapped_bam_file: bam of reads mapped to genome
        run_crispresso_genome_sites: boolean of whether to run crispresso on highly-aligned genomic locations (if false, crispresso will be run only at cut sites)
        crispresso_min_aln_score: minimum score for reads to align to amplicons
        samtools_command: location of samtools to run
        crispresso_command: location of crispresso to run
        query_bp_around_cut: extract reads that are within this bp of the cut site

    returns:
        cut_infos
        crispresso_infos: metadata about each crispresso run
            tuple of: name, chr, start, end, readCount, amplicon
        crispresso_commands: list of crispresso commands to run
        genome_read_count_at_cuts: dict for each target name of the number of reads discovered at that cut site (to be added to Linear class)
    """
    logging.info('Preparing global alignments for CRISPResso2')

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
    crispresso_commands = []

    if not os.path.isdir(root+'.CRISPResso_data'):
        os.mkdir(root+'.CRISPResso_data')

    genome_read_count_at_cuts = {}
    total_genome_read_count_at_cuts = 0
    for i,cut in enumerate(cuts):
        cut_els = cut.split(":")
        cut_chrom = cut_els[0]
        cut_loc = int(cut_els[1])
        cut_query_start = cut_loc - query_bp_around_cut #pull out reads in this region
        cut_query_end = cut_loc + query_bp_around_cut
        cut_start = cut_loc - av_read_length #define amplicon
        cut_end = cut_loc + av_read_length

        target_name = 'CRISPRlungo_WT'+str(i)
        name = target_name + "_" + cut_chrom + "_" + str(cut_loc)
        logging.debug('Analyzing cut target in genomic alignment ' + str(i) + ': ' + name)

        # -F 260 - filter read unmapped (4) or not primary alignment (240)
        read_count_command = '%s view -F 260 -c %s %s:%d-%d'%(samtools_command,genome_mapped_bam_file,cut_chrom,cut_query_start,cut_query_end)
        logging.debug(read_count_command)
        read_count = int(subprocess.check_output(read_count_command,shell=True).strip())
        logging.debug('got ' + str(read_count) + ' reads')
        total_genome_read_count_at_cuts += read_count
        genome_read_count_at_cuts[target_name] = read_count

        if read_count > crispresso_cutoff:
            amp_seq = subprocess.check_output(
                '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,cut_chrom,cut_start,cut_end),shell=True).decode('utf-8').strip()
            crispresso_infos.append((name,cut_chrom,cut_start,cut_end,amp_seq))

            read_count = 0
            reads_file = root + ".CRISPResso_data/"+name+".fastq"
            with open(reads_file,'w') as reads_out:
                for line in read_command_output('%s view -F 260 %s %s:%d-%d'%(samtools_command,genome_mapped_bam_file,cut_chrom,cut_query_start,cut_query_end)):
                    if line.strip() == "": break
                    line_els = line.split("\t")
                    seq = line_els[9]
                    qual = line_els[10]
                    read_id = line_els[0]
                    reads_out.write("%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                    read_count += 1
            crispresso_cmd = "%s -o %s -n %s --default_min_aln_score %d -a %s -r1 %s &> %s.log"%(crispresso_command,root + '.CRISPResso_runs',name,crispresso_min_aln_score,amp_seq,reads_file,reads_file)
            crispresso_commands.append(crispresso_cmd)

    logging.info('Got %d reads aligned in genome-wide alignment to cut sites.'%total_genome_read_count_at_cuts)

    #if we should run crispresso on highly-aligned genomic sites
    if run_crispresso_genome_sites:
        for chrom in aligned_locs.keys():
            starts = aligned_locs[chrom].keys()
            if len(starts) == 0:
                continue
            starts.sort()
            big_start_locs = [] #these regions can get pretty big if multiple sites get them
            big_end_locs = []
            last_ind = 0
            for start in starts:
                count = aligned_locs[chrom][start]
                if count > crispresso_cutoff:
                    if len(big_start_locs) == 0:
                        big_start_locs = [start]
                        big_end_locs = [start + av_read_length]
                    elif big_end_locs[last_ind] > start:
                        big_end_locs[last_ind] = start + av_read_length
                    else:
                        last_ind += 1
                        big_start_locs.append(start)
                        big_end_locs.append(start + av_read_length)
            if len(big_start_locs) == 0:
                continue

            #split big regions into smaller av_read_len-sized regions for analysis
            region_buffer = 30 #how much longer than av_read_len to allow amplicons to be
            start_locs = []
            end_locs = []
            for i in range(len(big_start_locs)):
                if big_end_locs[i] - big_start_locs[i] > av_read_length+region_buffer:
                    for j in range(big_start_locs[i],big_end_locs[i]-(av_read_length+region_buffer),int(av_read_length-region_buffer)):
                        start_locs.append(j)
                        end_locs.append(j+av_read_length+region_buffer)
                    start_locs.append(big_end_locs[i] - (av_read_length+region_buffer))
                    end_locs.append(big_end_locs[i])
                else:
                    start_locs.append(big_start_locs[i])
                    end_locs.append(big_end_locs[i])


            for i in range(len(start_locs)):
                name = chrom + "_" + str(start_locs[i])
                amp_seq = subprocess.check_output(
                        '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chrom,start_locs[i],end_locs[i]),shell=True).decode('utf-8').strip()

                read_count = 0
                reads_file = root + ".CRISPResso_data/"+name+".fastq"
                with open(reads_file,'w') as reads_out:
                    for line in read_command_output('%s view -F 260 %s %s:%d-%d'%(samtools_command,genome_mapped_bam_file,chrom,start_locs[i],end_locs[i])):
                        if line.strip() == "": break
                        line_els = line.split("\t")
                        seq = line_els[9]
                        qual = line_els[10]
                        read_id = line_els[0]
                        reads_out.write("%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                        read_count += 1
                crispresso_infos.append((name,chrom,start_locs[i],end_locs[i],read_count,amp_seq))
                crispresso_cmd = "%s -o %s -n %s --default_min_aln_score %d -a %s -r1 %s &> %s.log"%(crispresso_command,root + '.CRISPResso_runs',name,crispresso_min_aln_score,amp_seq,reads_file,reads_file)
                crispresso_commands.append(crispresso_cmd)
    logging.info('Created ' + str(len(crispresso_commands)) + ' CRISPResso commands from global alignment')
    return (crispresso_infos,crispresso_commands,genome_read_count_at_cuts)

def prep_crispresso2_artificial_targets(root,genome_len_file,crispresso_cutoff,custom_index_fasta,custom_mapped_bam_file,target_names,target_info,genome_read_count_at_cuts,crispresso_min_aln_score,samtools_command='samtools',crispresso_command='CRISPResso'):
    """
    Prepares data aligned to artificial targets for analysis with crispresso2
    Frequently-aligned targets with a min number of reads (crispresso_cutoff) are identified
    Reads from these locations are extracted and prepared for analysis by CRISPResso2

    params:
        root: root for written files
        genome_len_file: path to tab-sep file with lengths of chrs in genome (ends with .fai)
        crispresso_cutoff: min number of reads at site for crispresso2 processing
        custom_index_fasta: fasta of artificial targets
        custom_mapped_bam_file: aligned reads aligning to artificial targets
        target_names: list of artificial target names
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targets)
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
            target_info[target_name]['query_pos']: genomic start of query (bp)
            target_info[target_name]['query_start']: bp after which query starts (in case of padding)
            target_info[target_name]['query_end']
            target_info[target_name]['target_cut_idx']: bp of cut in target
        genome_read_count_at_cuts: number of reads discovered at cut sites in genomic alignment (to be added to 'Linear' class category)
        crispresso_min_aln_score: minimum score for reads to align to amplicons
        samtools_command: location of samtools to run
        crispresso_command: location of crispresso to run

    returns:
        crispresso_infos: array of metadata information for CRISPResso
            tuple of: name, chr, start, end, readCount, amplicon
        crispresso_commands: array of commands to run CRISPResso
    """
    logging.info('Preparing artificial targets for CRISPResso2')

    #first, add the counts from the genome alignment
    #these should all probably be 'Linear', but I'm leaving this generalized here in case things change
    class_counts = {}
    for target_name in genome_read_count_at_cuts:
        custom_class = target_info[target_name]['class']
        if custom_class not in class_counts:
            class_counts[custom_class] = 0
        class_counts[custom_class] += genome_read_count_at_cuts[target_name]

    crispresso_commands = []
    crispresso_infos = []

    target_summary_file = root + '.customTargetAlignment.summary.txt'
    with open (target_summary_file,'w') as out:
        out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('target_name','target_class','total_aligned','aligned_to_artificial_targets','aligned_to_genome','cut1_chr','cut1_site','cut2_chr','cut2_site'))
        for i,target_name in enumerate(target_names):
            start_ind = target_info[target_names[i]]['query_start']
            end_ind = target_info[target_names[i]]['query_end']

            logging.debug('Analyzing target ' + str(i) + ': ' + target_name)
            logging.debug('%s view -F 260 -c %s %s:%d-%d'%(samtools_command,custom_mapped_bam_file,target_name,start_ind,end_ind))
            read_count = int(subprocess.check_output('%s view -F 260 -c %s %s:%d-%d'%(samtools_command,custom_mapped_bam_file,target_name,start_ind,end_ind),shell=True).strip())
            logging.debug('got ' + str(read_count) + ' reads')

            custom_chr_class = target_info[target_name]['class']
            if custom_chr_class not in class_counts:
                class_counts[custom_chr_class] = 0
            class_counts[custom_chr_class] += read_count

            genome_aligned = 0
            if target_name in genome_read_count_at_cuts:
                genome_aligned = genome_read_count_at_cuts[target_name]
            total_aligned = genome_aligned + read_count
            out.write('%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n'%(target_name,target_info[target_name]['class'],total_aligned,read_count,genome_aligned,target_info[target_name]['cut1_chr'],str(target_info[target_name]['cut1_site']),target_info[target_name]['cut2_chr'],str(target_info[target_name]['cut2_site'])))

            if read_count > crispresso_cutoff:
                amp_seq = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,custom_index_fasta,target_name,start_ind,end_ind),shell=True).decode('utf-8').strip()

                reads_file = root + ".crispresso."+target_name+".fastq"
                printed_read_count = 0
                with open(reads_file,'w') as reads_out:
                    for line in read_command_output('%s view -F 260 %s %s:%d-%d'%(samtools_command,custom_mapped_bam_file,target_name,start_ind,end_ind)):
                        if line.strip() == "": break
                        line_els = line.split("\t")
                        seq = line_els[9]
                        qual = line_els[10]
                        read_id = line_els[0]
                        reads_out.write("%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                        printed_read_count += 1

                loc1 = target_info[target_name]['cut1_chr'] + "_" + str(target_info[target_name]['cut1_site'])
                loc2 = target_info[target_name]['cut2_chr'] + "_" + str(target_info[target_name]['cut2_site'])
                crispresso_infos.append((target_name,target_name,loc1,loc2,read_count,amp_seq))
                crispresso_commands.append("%s -o %s -n %s --default_min_aln_score %d -a %s -r1 %s &> %s.log"%(crispresso_command,root+'.CRISPResso_runs',target_name,crispresso_min_aln_score,amp_seq,reads_file,reads_file))

    class_counts_root = root + ".class_counts"
    keys = sorted(class_counts.keys())
    vals = [str(class_counts[key]) for key in keys]
    with open(class_counts_root+".txt",'w') as class_out:
        class_out.write("\t".join(keys)+"\n")
        class_out.write("\t".join(vals)+"\n")


    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    ax.pie(vals,labels=[keys[idx]+"\n("+str(vals[idx])+")" for idx in range(len(keys))],autopct="%1.2f%%")
    plt.savefig(class_counts_root+".pdf",pad_inches=1,bbox_inches='tight')

    logging.info('Created ' + str(len(crispresso_commands)) + ' CRISPResso commands from artificial targets')
    return (crispresso_infos,crispresso_commands)

def run_and_aggregate_crispresso(root,crispresso_infos,crispresso_commands):
    """
    Runs CRISPResso2 commands and aggregates output

    params:
        root: root for written files
        crispresso_infos: array of metadata information for CRISPResso
            tuple of: name, chr, start, end, readCount, amplicon
        crispresso_commands: array of commands to run CRISPResso

    returns:
        crispresso_infos: array of metadata information for CRISPResso
        crispresso_commands: array of commands to run CRISPResso

    """
    logging.info('Analyzing alignments using CRISPResso2')


    # getting rid of pickle in crispresso... so I'm leaving this import here in the meantime
    running_python3 = False
    if sys.version_info > (3, 0):
        running_python3 = True

    if running_python3:
            import pickle as cp #python 3
    else:
            import cPickle as cp #python 2.7
    # end nasty pickle part



    crispresso_results = {}
    crispresso_results['run_names'] = []
    crispresso_results['run_sub_htmls'] = {}
    crispresso_command_file = root + '.CRISPResso.commands.txt'
    crispresso_info_file = root + '.CRISPResso.info.txt'
    with open(crispresso_info_file,'w') as crispresso_info_file, open(crispresso_command_file,'w') as crispresso_command_file:
        crispresso_info_file.write('name\tchr\tstart\tend\treadCount\tamplicon\tn_total\treads_aligned\treads_unmod\treads_mod\treads_discarded\treads_insertion\treads_deletion\treads_substitution\treads_only_insertion\treads_only_deletion\treads_only_substitution\treads_insertion_and_deletion\treads_insertion_and_substitution\treads_deletion_and_substitution\treads_insertion_and_deletion_and_substitution\n')
        for idx,command in enumerate(crispresso_commands):
            crispresso_command_file.write(command)
            logging.debug('Running CRISPResso command '+str(idx))
            crispresso_output = subprocess.check_output(
                    command,shell=True).decode('utf-8').strip()
            name = crispresso_infos[idx][0]
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

            run_file = os.path.join(root + '.CRISPResso_runs','CRISPResso_on_'+name,'CRISPResso2_info.pickle')
            if os.path.isfile(run_file):
                run_data = cp.load(open(run_file,'rb'))
                report_filename = run_data['report_filename']
                report_file_loc = os.path.join(root+'.CRISPResso_runs',report_filename)
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
            new_vals = [n_total,n_aligned,n_unmod,n_mod,n_discarded,n_insertion,n_deletion,n_substitution,n_only_insertion,n_only_deletion,n_only_substitution,n_insertion_and_deletion,n_insertion_and_substitution,n_deletion_and_substitution,n_insertion_and_deletion_and_substitution]
            crispresso_info_file.write("\t".join([str(x) for x in crispresso_infos[idx]])+"\t"+"\t".join([str(x) for x in new_vals])+"\n")
        return crispresso_results

def cleanup(root):
    """
    Deletes intermediate files

    params:
        root: root for written files
    """
    delete_result = subprocess.check_output('rm -rf ' + root + '.customIndex.fa.*', stderr=subprocess.STDOUT,shell=True)
    logging.debug('Deleted bowtie indexes ' + delete_result)

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
