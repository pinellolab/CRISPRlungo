import os
import sys
import subprocess
import logging
import re

def main():
    settings = parse_settings(sys.argv)
    assert_dependencies()
    av_read_length = get_av_read_len(settings['fastq'])


    primer_chr = ""
    primer_loc = -1
    if 'primer_site' in settings:
        primer_info = settings['primer_site'].split("_")
        primer_chr = primer_info[0]
        primer_loc = int(primer_info[1])

    (target_names,target_info) = make_artificial_targets(
            cuts=settings['cuts'],
            genome=settings['genome'],
            target_length=av_read_length+settings['alignment_extension'],
            primer_chr=primer_chr,
            primer_loc=primer_loc)

    input_for_genome_mapping = settings['fastq'] #if targets are provided, the input for genome mapping will be the reads that don't align to targets
    custom_aligned_count = 0
    crispresso_commands = [] #list of crispresso commands to run -- first, we add the commands from the custom targets, then the genome aligned reads
    crispresso_infos = [] #meta info about the crispresso runs

    if len(target_names) > 0:
        (input_for_genome_mapping,custom_aligned_count, custom_mapped_bam_file, custom_index_fasta
            ) = align_to_artificial_targets(
                    root = settings['root'],
                    fastq = settings['fastq'],
                    target_names = target_names,
                    target_info = target_info
                    )

        (crispresso_infos, crispresso_commands
            ) = prep_crispresso2_artificial_targets(
                    root = settings['root'],
                    genome_len_file=re.sub('.fa$','.dict',settings['genome']),
                    crispresso_cutoff=settings['crispresso_cutoff'],
                    av_read_length = av_read_length,
                    alignment_extension=settings['alignment_extension'],
                    custom_index_fasta = custom_index_fasta,
                    custom_mapped_bam_file = custom_mapped_bam_file,
                    target_names = target_names,
                    target_info = target_info,
                    query_bp_around_cut=50
                    )

    (genome_mapped_bam_file, genome_aligned_count
        ) = align_to_genome(
                root = settings['root'],
                fastq = input_for_genome_mapping,
                bowtie2_genome = settings['bowtie2_genome']
                )

    (mapped_chopped_sam_file, aligned_locs, mapped_chrs, max_frags
        ) = analyze_global_aln(
                root = settings['root'],
                genome_mapped_bam_file = genome_mapped_bam_file,
                bowtie2_genome = settings['bowtie2_genome'],
                fragment_size=settings['fragment_size'],
                fragment_step_size=settings['fragment_step_size']
                )

    analyze_chopped_reads(
            root = settings['root'],
            mapped_chopped_sam_file=mapped_chopped_sam_file,
            mapped_chrs = mapped_chrs,
            max_frags = max_frags)

    (genome_crispresso_infos, genome_crispresso_commands
            ) = prep_crispresso2_global(
                root = settings['root'],
                genome_len_file = re.sub('.fa$','.dict',settings['genome']),
                crispresso_cutoff = settings['crispresso_cutoff'],
                aligned_locs = aligned_locs,
                av_read_length = av_read_length,
                genome_mapped_bam_file = genome_mapped_bam_file
                )

    crispresso_commands.extend(genome_crispresso_commands)
    crispresso_infos.extend(genome_crispresso_infos)

    crispresso_info_file = settings['root'] + '.CRISPResso.info'
    with open(crispresso_info_file,'w') as crispresso_info_file:
        crispresso_info_file.write('name\tchr\tstart\tend\treadCount\tamplicon\n')
        for info in crispresso_infos:
            crispresso_info_file.write(info)

    crispresso_command_file = settings['root'] + '.CRISPResso.commands'
    with open(crispresso_command_file,'w') as crispresso_command_file:
        for idx,command in enumerate(crispresso_commands):
            crispresso_command_file.write(command)
            logging.debug('Running CRISPResso command '+str(idx))
            crispresso_output = subprocess.check_output(
                    command,shell=True).decode('utf-8').strip()



def parse_settings(args):
    if len(args) < 2:
        raise Exception('Error: Settings file is missing\nUsage: ' +
                    args[0] + ' {settings file}')
    settings_file = args[1]

    logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s:%(levelname)s: %(message)s",
            handlers=[
                logging.FileHandler(settings_file+".log"),
                logging.StreamHandler()
                ]
            )

    settings = {}
    with open(settings_file, 'r') as fin:
        for line in fin.readlines():
            line_els = line.split("#")[0].strip().split("\t")
            key = line_els[0].strip()
            val = line_els[1].strip()
            settings[key] = val


    # two parameters control the size and stride of fragment creation
    # size of fragment
    settings['fragment_size'] = int(settings['fragment_size']) if 'fragment_size' in settings else 20
    # step size between fragments
    settings['fragment_step_size'] = int(settings['fragment_step_size']) if 'fragment_step_size' in settings else 10

    # number of bp to extend beyond av read length around cut site for custom index
    settings['alignment_extension'] = settings['alignment_extension'] if 'alignment_extension' in settings else 50

    settings['crispresso_cutoff'] = settings['crispresso_cutoff'] if 'crispresso_cutoff' in settings else 50

    settings['cuts'] = settings['cut_sites'].split(" ") if 'cut_sites' in settings else []

    fastq = settings['fastq']
    if not os.path.isfile(fastq):
        raise Exception('Error: Fastq file %s does not exist',fastq)

    settings['root'] = settings_file + ".CRISPRlungo"

    return settings


def assert_dependencies():
    """
    Asserts the presence of required software (faidx, bowtie2)
    """
    # check faidx
    try:
        faidx_result = subprocess.check_output('samtools faidx', stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: samtools faidx is required')
    if not 'Usage: samtools faidx' in str(faidx_result):
        raise Exception('Error: samtools faidx is required')
    #check bowtie2
    try:
        bowtie_result = subprocess.check_output('bowtie2 --version', stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: bowtie2 is required')

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
        number_reads_to_check: the number of reads to read in

    returns:
        av_read_len: average read length
    """
    sum_len = 0
    with open(fastq,'r') as fin:
        for i in range(number_reads_to_check):
            fastq_id = fin.readline()
            seq = fin.readline()
            read_len = len(seq)
            sum_len += read_len
            fastq_plus = fin.readline()
            fastq_qual = fin.readline()
    av_read_len = int(sum_len/float(number_reads_to_check))
    return(av_read_len)

def make_artificial_targets(cuts,genome,target_length,primer_chr ="",primer_loc=-1):
    """
    Generates fasta sequences surrounding cuts for alignment
    At each cut point, sequence of length target_length is generated
    Combinations of sequences at each cut site are produced
    If a primer is given, only cut sites that include the primer are used

    params:
        cuts: array of cut locations
        genome: location of fastq genome
        target_length: how long the resulting fragment should be (on one side fothe cut)
        primer_chr: chromosome of primer (if any)
        primer_loc: location of primer (if any)

    returns:
        target_names: array of target names (corresponding to targets)
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targest
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
    """
    target_names = []
    target_info = {}

    for i in range(len(cuts)):
        chr_els = cuts[i].split("_")
        chr_A = chr_els[0]
        site_A = int(chr_els[1])
        cut_start_A = site_A - target_length
        cut_start_A_stop = site_A - 1
        cut_end_A = site_A + target_length

        primer_is_in_left_bit_A = (primer_chr != "" and primer_chr == chr_A and primer_loc >= cut_start_A and primer_loc <= cut_start_A_stop)
        primer_is_in_right_bit_A = (primer_chr != "" and primer_chr == chr_A and primer_loc >= site_A and primer_loc <= cut_end_A)

        left_bit_A = subprocess.check_output(
                'samtools faidx -n 10000 %s %s:%d-%d | tail -n 1'%(genome,chr_A,cut_start_A,cut_start_A_stop),shell=True).decode('utf-8').strip()

        right_bit_A = subprocess.check_output(
                'samtools faidx -n 10000 %s %s:%d-%d | tail -n 1'%(genome,chr_A,site_A,cut_end_A),shell=True).decode('utf-8').strip()

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
                    'cut2_site':site_A
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
                    'cut2_site':site_A
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
                    'cut2_site':site_A
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
                    'cut2_site':site_A
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
                    'cut2_site':site_A
                    }


        for j in range(i+1,len(cuts)):

            chr_els = cuts[j].split("_")
            chr_B = chr_els[0]
            site_B = int(chr_els[1])
            cut_start_B = site_B - target_length
            cut_start_B_stop = site_B - 1
            cut_end_B = site_B + target_length

            primer_is_in_left_bit_B = (primer_chr != "" and primer_chr == chr_B and primer_loc >= cut_start_B and primer_loc <= cut_start_B_stop)
            primer_is_in_right_bit_B = (primer_chr != "" and primer_chr == chr_B and primer_loc >= site_B and primer_loc <= cut_end_B)

            left_bit_B = subprocess.check_output(
                    'samtools faidx -n 10000 %s %s:%d-%d | tail -n 1'%(genome,chr_B,cut_start_B,cut_start_B_stop),shell=True).decode('utf-8').strip()

            right_bit_B = subprocess.check_output(
                    'samtools faidx -n 10000 %s %s:%d-%d | tail -n 1'%(genome,chr_B,site_B,cut_end_B),shell=True).decode('utf-8').strip()


            if primer_chr == "" or primer_is_in_left_bit_A or primer_is_in_right_bit_B:
                LARB = left_bit_A + right_bit_B
                target_name = 'CRISPRlungo_L' + str(i) + 'R' + str(j)
                target_names.append(target_name)
                target_info[target_name] = {
                        'sequence': LARB,
                        'cut1_chr':chr_A,
                        'cut1_site':site_A,
                        'cut2_chr':chr_B,
                        'cut2_site':site_B
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
                        'cut2_site':site_B
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
                        'cut2_site':site_B
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
                        'cut2_site':site_B
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
                        'cut2_site':site_B
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
                        'cut2_site':site_B
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
                        'cut2_site':site_B
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
                        'cut2_site':site_B
                        }
                if chr_A == chr_B:
                    target_info[target_name]['class'] = 'Large deletion'
                else:
                    target_info[target_name]['class'] = 'Translocation'


    return(target_names,target_info)

def align_to_artificial_targets(root,fastq,target_names,target_info):
    """
    Aligns reads to artificial targets

    params:
        root: root for written files
        fastq: fastq to align
        target_names: array of target names
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targest
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']

    returns:
        custom_unmapped_fastq_file: fastq file of reads not aligning to artificial targets
        custom_aligned_count: number of reads aligned to artificial targets
        custom_mapped_bam_file: aligned reads aligning to artificial targets
        custom_index_fasta: fasta of artificial targets

    """

    custom_mapped_bam_file = root + ".customMapped.bam" # file only generated if cut sites provided
    custom_index_fasta = root + ".customIndex.fa"
    logging.info('Printing ' + str(len(target_names)) + ' targets to custom index (' + custom_index_fasta + ')')
    with open(custom_index_fasta,'w') as fout:
        for i in range(len(target_names)):
            fout.write('>'+target_names[i]+'\n'+target_info[target_names[i]]['sequence']+'\n')

    if not os.path.isfile(custom_index_fasta + '.1.bt2'):
        logging.info('Indexing custom targets using bowtie2-build (' + custom_index_fasta + ')')
        index_result = subprocess.check_output('bowtie2-build ' + custom_index_fasta + ' ' + custom_index_fasta, shell=True,stderr=subprocess.STDOUT)

    custom_unmapped_fastq_file = root + '.customUnmapped.fastq'
    custom_bowtie_log = root + '.customBowtie2Log'
    if not os.path.isfile(custom_mapped_bam_file):
        logging.info('Aligning reads to custom targets using bowtie2')
        custom_aln_command = 'bowtie2 -x %s -U %s --un-gz %s 2> %s | samtools view -Shu - | samtools sort -o %s - && samtools index %s' % (custom_index_fasta,fastq,custom_unmapped_fastq_file,custom_bowtie_log,custom_mapped_bam_file,custom_mapped_bam_file)
        aln_result = subprocess.check_output(custom_aln_command,shell=True,stderr=subprocess.STDOUT) 

    custom_aligned_count = int(subprocess.check_output('samtools view -F 4 -c %s'%custom_mapped_bam_file,shell=True).strip())

    logging.info('Aligned ' + str(custom_aligned_count) + ' reads to custom targets')
    return(custom_unmapped_fastq_file,custom_aligned_count,custom_mapped_bam_file,custom_index_fasta)



def align_to_genome(root,fastq,bowtie2_genome):
    """
    Aligns reads to entire genome

    params:
        root: root for written files
        fastq: fastq to align
        genome: location of bowtie2 genome files (root of files)

    returns:
        genome_aligned_bam: bam of reads after alignment to genome
        genome_aligned_count: number of reads aligning to genome
    """

    genome_mapped_bam_file = root + ".genomeMapped.bam"
    bowtie_log = root + '.bowtie2Log'
    aln_command = 'bowtie2 -x %s -U %s 2> %s | samtools view -Shu - | samtools sort -o %s - && samtools index %s'%(bowtie2_genome,fastq,bowtie_log,genome_mapped_bam_file,genome_mapped_bam_file)

    aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT,encoding='utf-8') 
    logging.debug('Alignment to genome: ' + aln_result)

    genome_aligned_count = int(subprocess.check_output('samtools view -F 4 -c %s'%genome_mapped_bam_file,shell=True).strip())

    logging.info('Aligned ' + str(genome_aligned_count) + ' reads to genome')
    return(genome_mapped_bam_file,genome_aligned_count)


def run_command(command):
    """
    Runs a shell command and returns an iter to read the output

    param: 
        command: shell command to run

    returns:
        iter to read the output
    """

    p = subprocess.Popen(command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,shell=True)#,
#            encoding='utf-8',universal_newlines=True)
    return iter(p.stdout.readline, b'')

def analyze_global_aln(root,genome_mapped_bam_file,bowtie2_genome,fragment_size,fragment_step_size):
    """
    Analyzes reads aligned globally
    Keeps track of locations where reads map
    For reads that do not map, creats chopped fragments which are be mapped

    params:
        root: root for written files
        genome_mapped_bam_file: bam of reads after alignment to genome (some will be unaligned)
        bowtie2_genome: location of bowtie2 indexed reference genome
        fragment_size: size of resulting chopped fragments. 
        fragment_step_size: stride between chopped fragments. 5 results in a new fragment every 5bp



    returns:
        mapped_chopped_sam_file: mapped chopped fragments
        aligned_locs: hash of locations aligned to on each chromosome aligned_locs[chr][position] = count
        mapped_chrs: hash of counts of reads aligned to each chr in the genome-alignment step
        max_frags: max number of fragments created per read
    """

    mapped_chrs = {}
    unmapped_ids = {}
    unmapped_id = 0
    aligned_locs = {} # aligned reads will not be chopped, but keep track of where they aligned for CRISPResso output
    aligned_chr_counts = {}
    max_frags = 0 #max frags produced for a read
    frags_per_read = {} #number of fragments created per read
    unmapped_frag_file = root + ".unmapped.fastq"
    with open(unmapped_frag_file,"w") as unmapped_fastq:
        for line in run_command('samtools view %s'%genome_mapped_bam_file):
            if line.strip() == "": break
            line_els = line.split("\t")
            line_chr = line_els[2]
            if not line_chr in mapped_chrs:
                mapped_chrs[line_chr] = 0
            mapped_chrs[line_chr] += 1

            line_mapq = line_els[5]
            line_unmapped = int(line_els[1]) & 0x4
            line_start = int(line_els[3])-1

            if line_unmapped:
                line_seq = line_els[9]
                line_qual = line_els[10]
                line_id = line_els[0]
                frag_num = 0
                offset = 0
                while offset + fragment_size < len(line_seq):
                    new_id = "@%s.CLID%s.CLO%s.CLF%s"%(line_id,unmapped_id,offset,frag_num)
                    new_seq =  line_seq[offset:offset + fragment_size]
                    new_qual = line_qual[offset:offset + fragment_size]
                    unmapped_fastq.write("%s\n%s\n+\n%s\n"%(new_id,new_seq,new_qual))
                    frag_num += 1
                    offset += fragment_step_size
                
                offset = len(line_seq)-fragment_size
                new_id = "@%s.CLID%s.CLO%s.CLF%s"%(line_id,unmapped_id,offset,frag_num)
                new_seq =  line_seq[offset:offset + fragment_size]
                new_qual = line_qual[offset:offset + fragment_size]
                unmapped_fastq.write("%s\n%s\n+\n%s\n"%(new_id,new_seq,new_qual))
                if frag_num > max_frags:
                    max_frags = frag_num
                if frag_num not in frags_per_read:
                    frags_per_read[frag_num] = 0
                frags_per_read[frag_num] += 1

                unmapped_id += 1
            else:
                if line_chr not in aligned_locs:
                    aligned_locs[line_chr] = {}
                if line_start not in aligned_locs[line_chr]:
                    aligned_locs[line_chr][line_start] = 0

                aligned_locs[line_chr][line_start] += 1
                if line_chr not in aligned_chr_counts:
                    aligned_chr_counts[line_chr] = 0
                aligned_chr_counts[line_chr] += 1

    number_unmapped_reads_chopped = unmapped_id
    logging.info('Created chopped reads for ' + str(number_unmapped_reads_chopped) + ' reads')

    frags_per_unaligned_read_file = root + ".fragsPerUnalignedRead"
    with open(frags_per_unaligned_read_file,"w") as frags:
        frags.write('numFragments\tnumReads')
        for key in sorted(frags_per_read.keys()):
            frags.write(str(key) + '\t' + str(frags_per_read[key]) + '\n')

    global_aligned_chrs_file = root + ".globalAlignedChrs"
    with open(global_aligned_chrs_file,"w") as chrs:
        chrs.write('chr\tnumReads\n')
        for key in sorted(aligned_chr_counts.keys()):
            chrs.write(key + '\t' + str(aligned_chr_counts[key]) + '\n')

    mapped_chopped_sam_file = root + ".fragMapped.sam"
    chopped_bowtie_log = root + ".fragMapped.bowtie2Log"
    chopped_aln_command = 'bowtie2 -x %s -U %s -S %s --end-to-end 2> %s' %(bowtie2_genome,unmapped_frag_file,mapped_chopped_sam_file,chopped_bowtie_log)

    aln_result = subprocess.check_output(chopped_aln_command,shell=True,stderr=subprocess.STDOUT,encoding='utf-8')
    logging.debug('Alignment of chopped fragments to genome: ' + aln_result)


    return(mapped_chopped_sam_file,aligned_locs,mapped_chrs,max_frags)

def analyze_chopped_reads(root,mapped_chopped_sam_file,mapped_chrs,max_frags):
    """
    Analyzes reads aligned globally
    Keeps track of locations where reads map
    For reads that do not map, creats chopped fragments which are be mapped

    params:
        root: root for written files
        mapped_chopped_sam_file: file with mapped chopped fragments
        mapped_chrs: counts of reads aligned to each chr in the genome-alignment step
        max_frags: max number of fragments a read was split into


    returns:
    TODO
    """
    frag_meta_file = root + ".fragMeta"
    with open(frag_meta_file,"w") as frag_file:
        head = "ID"
        for i in range(max_frags):
            head += "\t"+str(i)
        frag_file.write(head+"\n")

        with open(mapped_chopped_sam_file,'r') as aln_frags:
            #read header
            line1 = aln_frags.readline()
            while line1 and not line1.startswith('@PG'):
                line1 = aln_frags.readline()

            curr_id = ""
            curr_id_chrs = {} # keep track of which chrs frags from this id map to
            curr_id_vals = {} # keep track of where each frag maps to (this is a hash but it's like an array, so we can check for uninitialized values)
            translocations = {} # hash of all translocation locations
            translocation_count = 0
            unidentified_count = 0
            chroms_per_frag_read_count = {} # keep track for all read how many chrs the fragments mapped to
            frags_mapped_count = 0
            for line in aln_frags.readlines():
                line_els = line.split("\t")
                line_chr = line_els[2]
                if line_chr not in mapped_chrs:
                    mapped_chrs[line_chr] = 0
                mapped_chrs[line_chr] += 1
                frags_mapped_count += 1

                line_mapq = line_els[5]
                line_unmapped = int(line_els[1]) & 0x4
                line_start = int(line_els[3])-1

                match = re.match('(.*)\.CLID(\d+)\.CLO(\d+)\.CLF(\d+)$',line_els[0])
                if not match:
                    raise Exception('Cannot parse id %s from line %s. Perhaps line was trimmed?\n'%(line_els[0],line))
                (orig_id,lungo_id,lungo_offset,lungo_frag) = match.groups()

                if curr_id != orig_id and curr_id != "":
                    outline = curr_id;
                    for i in range(max_frags):
                        val = ""
                        if i in curr_id_vals:
                            val = curr_id_vals[i]
                            val_chr,val_pos = curr_id_vals[i].split(" ")
                            if val_chr != "*":
                                if val_chr not in curr_id_chrs:
                                    curr_id_chrs[val_chr] = 0
                                curr_id_chrs[val_chr] += 1
                        outline += "\t"+val
                    frag_file.write(outline)

                    if '0' in curr_id_vals and curr_id_vals['0'] == "*":
                        unidentified_count += 1
                    elif 'last_chr' in curr_id_vals and curr_id_vals['last_chr'] == "*":
                        unidentified_count += 1
                    else:
                        (start_chr,start_pos) = curr_id_vals['0'].split(" ")
                        key = start_chr + " " + curr_id_vals['last_chr']
                        if key not in translocations:
                            translocations[key] = 0
                        translocations[key] += 1
                        translocation_count += 1

                    curr_chr_count = len(curr_id_chrs.keys())
                    if curr_chr_count not in chroms_per_frag_read_count:
                        chroms_per_frag_read_count[curr_chr_count] = 0
                    chroms_per_frag_read_count[curr_chr_count] += 1

                    curr_id_chrs = {}
                    curr_id_vals = {}

                curr_id = orig_id

                inferred_start = line_start - int(lungo_offset)
                curr_id_vals[lungo_frag] = line_chr + " " + str(inferred_start)
                curr_id_vals['last'] = line_chr + " " + str(inferred_start)
                curr_id_vals['last_chr'] = line_chr

            #finish last item
            if frags_mapped_count > 0:
                outline = curr_id;
                for i in range(max_frags):
                    val = ""
                    if i in curr_id_vals:
                        val = curr_id_vals[i]
                        val_chr,val_pos = curr_id_vals[i].split(" ")
                        if val_chr != "*":
                            if val_chr not in curr_id_chrs:
                                curr_id_chrs[val_chr] = 0
                            curr_id_chrs[val_chr] += 1
                    outline += "\t"+val
                frag_file.write(outline)

                if '0' in curr_id_vals and curr_id_vals['0'] == "*":
                    unidentified_count += 1
                elif 'last_chr' in curr_id_vals and curr_id_vals['last_chr'] == "*":
                    unidentified_count += 1
                else:
                    (start_chr,start_pos) = curr_id_vals['0'].split(" ")
                    key = start_chr + " " + curr_id_vals['last_chr']
                    if key not in translocations:
                        translocations[key] = 0
                    translocations[key] += 1
                    translocation_count += 1

                curr_chr_count = len(curr_id_chrs.keys())
                if curr_chr_count not in chroms_per_frag_read_count:
                    chroms_per_frag_read_count[curr_chr_count] = 0
                chroms_per_frag_read_count[curr_chr_count] += 1
            #done with last one

    logging.info("Found %d translocations and %d unidentified reads"%(translocation_count,unidentified_count))

    frags_aligned_chrs_file = root + ".fragsAlignedChrs"
    with open(frags_aligned_chrs_file,"w") as fout:
        fout.write("chr\tnumReads\n")
        for key in sorted(mapped_chrs.keys()):
            fout.write(str(key) + '\t' + str(mapped_chrs[key]) + '\n')

    frag_chroms_per_read_file = root + ".fragsChromsPerRead"
    with open(frag_chroms_per_read_file,"w") as fout:
        fout.write("numChroms\tnumReads\n")
        for key in sorted(chroms_per_frag_read_count.keys()):
            fout.write(str(key) + '\t' + str(chroms_per_frag_read_count[key]) + '\n')

    #make translocation table
    translocation_table = {}
    translocation_report_file = root + ".translocationReport"
    with open(translocation_report_file,"w") as fout:
        fout.write("from\tto\tcount\n")
        for key in sorted(translocations.keys()):
            (from_chr,to_chr) = key.split(" ")
            fout.write(from_chr+'\t'+to_chr+ '\t' + str(translocations[key]) + '\n')
            if from_chr not in translocation_table:
                translocation_table[from_chr] = {}
            if to_chr not in translocation_table[from_chr]:
                translocation_table[from_chr][to_chr] = 0
            translocation_table[from_chr][to_chr] += translocations[key]

    translocation_report_table_file = root + ".translocationReport.table"
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



def prep_crispresso2_global(root,genome_len_file,crispresso_cutoff,aligned_locs,av_read_length,genome_mapped_bam_file):
    """
    Prepares globally-aligned data for crispresso2
    Frequently-aligned locations with a min number of reads (crispresso_cutoff) are identified
    Reads from these locations are extracted and prepared for analysis by CRISPResso2

    params:
        root: root for written files
        genome_len_file: path to tab-sep file with lengths of chrs in genome (ends with .dict)
        crispresso_cutoff: min number of reads at site for crispresso2 processing
        aligned_locs: hash of locations aligned to on each chromosome aligned_locs[chr][position] = count
        av_read_length: length of average sequencing read in sample
        genome_mapped_bam_file: bam of reads mapped to genome

    returns:
        crispresso_infos: metadata about each crispresso run
        crispresso_commands: list of crispresso commands to run
    """
    logging.info('Preparing for CRISPResso2')

    chrom_lens = {}
    chroms = []
    with open(genome_len_file,'r') as gfile:
        for line in gfile:
            line_els = line.split("\t")
            line_chrom = line_els[0]
            if not re.match('^SN:chr',line_chrom):
                continue
            line_chrom = re.sub('SN:','',line_chrom)
            line_len = line_els[2]
            line_len = re.sub('LN:','',line_len)
            chrom_lens[line_chrom] = line_len
            chroms.append(line_chrom)

        crispresso_names = []
        crispresso_commands = []
        crispresso_infos = []

        for chrom in chroms:
            start_locs = [0]
            end_locs = [0]
            last_ind = 0
            starts = aligned_locs[chrom].keys()
            starts.sort()
            for start in starts:
                count = aligned_locs[chrom][start]
                if count > crispresso_cutoff:
                    if end_locs[last_ind] > start:
                        end_locs[last_ind] = start + av_read_length
                    else:
                        last_ind += 1
                        start_locs[last_ind] = start
                        end_locs[last_ind] = start + av_read_length
            for i in range(len(start_locs)):
                name = chrom + "_" + start_locs[i]
                crispresso_names.append(name)
                amp_seq = subprocess.check_output(
                        'samtools faidx -n 10000 %s %s:%d-%d | tail -n 1'%(genome,chrom,start_locs[i], end_locs[i]),shell=True).decode('utf-8').strip()

                read_count = 0
                reads_file = root + ".crispresso."+name+".fastq"
                with open(reads_file,'w') as reads_out:
                    for line in run_command('samtools view -F 4 %s %s:%d-%d'%(genome_mapped_bam_file,chrom,start_locs[i],end_locs[i])):
                        if line.strip() == "": break
                        line_els = line.split("\t")
                        seq = line_els[9]
                        qual = line_els[10]
                        read_id = line_els[0]
                        reads_out.write("%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                        read_count += 1
                crispresso_infos.append("%s\t%s\t%d\t%d\t$d\t$s\n"%(name,chrom,start_locs[i],end_locs[i],read_count,amp_seq))
                crispresso_command = "CRISPResso -n %s -a %s -r1 %s &> %s.log"%(name,amp_seq,reads_file,reads_file)
                crispresso_commands.append(crispresso_command)
    return (crispresso_infos,crispresso_commands)




def prep_crispresso2_artificial_targets(root,genome_len_file,crispresso_cutoff,av_read_length,alignment_extension,custom_index_fasta,custom_mapped_bam_file,target_names,target_info,query_bp_around_cut=50):
    """
    Prepares data aligned to artificial targets for analysis with crispresso2
    Frequently-aligned targets with a min number of reads (crispresso_cutoff) are identified
    Reads from these locations are extracted and prepared for analysis by CRISPResso2

    params:
        root: root for written files
        genome_len_file: path to tab-sep file with lengths of chrs in genome (ends with .dict)
        crispresso_cutoff: min number of reads at site for crispresso2 processing
        av_read_length: length of average sequencing read in sample
        alignment_extension: window around cut sites for artificial targets
        custom_index_fasta: fasta of artificial targets
        custom_mapped_bam_file: aligned reads aligning to artificial targets
        target_names: list of artificial target names
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targest
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
        query_bp_around_cut: extract reads that are within this bp of the cut site

    returns:
        crispresso_infos: array of metadata information for CRISPresso
        crispresso_commands: array of commands to run CRISPresso

    """
    logging.info('Preparing artificial targets for CRISPResso2')

    class_counts = {}

    #because the custom amplicons are the same length, the start and stop ind will be the same for all amplicons
    start_ind = (av_read_length + alignment_extension) - query_bp_around_cut
    end_ind = (av_read_length + alignment_extension) + query_bp_around_cut

    crispresso_commands = []
    crispresso_infos = []

    for i in range(len(target_names)):
        custom_chr_name = target_names[i]
        logging.debug('Analyzing target ' + str(i) + ': ' + custom_chr_name)
        logging.debug('samtools view -F 4 -c %s %s:%d-%d'%(custom_mapped_bam_file,custom_chr_name,start_ind,end_ind))
        read_count = int(subprocess.check_output('samtools view -F 4 -c %s %s:%d-%d'%(custom_mapped_bam_file,custom_chr_name,start_ind,end_ind),shell=True).strip())

        custom_chr_class = target_info[custom_chr_name]['class']
        if custom_chr_class not in class_counts:
            class_counts[custom_chr_class] = 0
        class_counts[custom_chr_class] += read_count


        if read_count > crispresso_cutoff:
            amp_seq = subprocess.check_output(
                'samtools faidx -n 10000 %s %s:%d-%d | tail -n 1'%(custom_index_fasta,custom_chr_name,start_ind,end_ind),shell=True).decode('utf-8').strip()

            reads_file = root + ".crispresso."+custom_chr_name+".fastq"
            printed_read_count = 0
            with open(reads_file,'w') as reads_out:
                for line in run_command('samtools view -F 4 %s %s:%d-%d'%(custom_mapped_bam_file,custom_chr_name,start_ind,end_ind)):
                    if line.strip() == "": break
                    line_els = line.split("\t")
                    seq = line_els[9]
                    qual = line_els[10]
                    read_id = line_els[0]
                    reads_out.write("%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                    printed_read_count += 1

            loc1 = target_info[custom_chr_name]['cut1_chr'] + "_" + str(target_info[custom_chr_name]['cut1_site'])
            loc2 = target_info[custom_chr_name]['cut2_chr'] + "_" + str(target_info[custom_chr_name]['cut2_site'])
            crispresso_infos.append("%s\t%s\t%s\t%s\t%d\t%s\n"%(custom_chr_name,custom_chr_name,loc1,loc2,read_count,amp_seq))
            crispresso_commands.append("CRISPResso -n %s -a %s -r1 %s &> %s.log"%(custom_chr_name,amp_seq,reads_file,reads_file))

    class_log_file = root + ".class_counts"
    with open(class_log_file,'w') as class_out:
        head_line = ""
        first_line = ""
        keys = sorted(class_counts.keys())
        vals = [str(class_counts[key]) for key in keys]
        class_out.write("\t".join(keys)+"\n")
        class_out.write("\t".join(vals)+"\n")

    return (crispresso_infos,crispresso_commands)




main()
