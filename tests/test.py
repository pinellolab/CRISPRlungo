import sys
import os
import glob
import unittest

test_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.dirname(test_dir))
from CRISPRlungo import CRISPRlungoCore as cl

this_bowtie2_genome = os.path.join(test_dir,'genome.chr11_5225364_5225863/genome')
this_genome = this_bowtie2_genome + '.fa'
test_inputs_folder = os.path.join(test_dir,'testInputs')
expected_outputs_folder = os.path.join(test_dir,'expectedResults')
temp_outputs_folder = os.path.join(test_dir,'tmpOutputs')

if not os.path.exists(this_genome):
    raise Exception(f'Missing genome file at {this_genome}')
if not os.path.exists(test_inputs_folder):
    raise Exception(f'Missing test inputs folder at {test_inputs_folder}')
if not os.path.exists(expected_outputs_folder):
    raise Exception(f'Missing expected outputs folder at {expected_outputs_folder}')
if not os.path.exists(temp_outputs_folder):
    os.makedirs(temp_outputs_folder)

verbose = True

class FunctionTests(unittest.TestCase):
    def test_prep_input(self):
        assert_prep_input(self,this_genome,this_bowtie2_genome,test_inputs_folder,verbose)
    def test_collapse_cuts(self):
        assert_collapse_cuts(self,verbose)

class RunTests(unittest.TestCase):
    def test_primer(self):
        assert_primer(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose)
    
    def test_dedup(self):
        assert_dedup(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose)

    def test_cuts(self):
        assert_cuts(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose)

    def test_no_results(self):
        assert_no_results(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose)

    def test_r1r2_agreement(self):
        assert_r1r2_agreement(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose)

    def test_runs(self):
        assert_runs(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose)

def assert_prep_input(self, this_genome, this_bowtie2_genome, test_inputs_folder, verbose):
    #test primer fw, guide fw
    print('Testing Primer FW Guide FW')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, av_read_length, num_reads_input = cl.prep_input(
                    root = 'test.py.primer_fw_guide_fw',
                    primer_seq = 'TTGCAATGAAAATAAATGTTT',
                    min_primer_length = 0,
                    guide_seqs = ['CTCAAGGCCCTTCATAATAT'],
                    cleavage_offset = -3,
                    fastq_r1 = os.path.join(test_inputs_folder,'primer_fw.fq'),
                    samtools_command = 'samtools',
                    genome = this_genome,
                    bowtie2_command = 'bowtie2',
                    bowtie2_genome = this_bowtie2_genome,
                    can_use_previous_analysis = False,
                    suppress_file_output = True,
                    )
    if verbose:
        print(f"{origin_seq=}")
        print(f"{cut_sites=}")
        print(f"{cut_annotations=}")
        print(f"{primer_chr=}")
        print(f"{primer_loc=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    self.assertEqual(origin_seq,'TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAA')
    self.assertEqual(primer_chr,'chr11_sm')
    self.assertEqual(primer_loc,101)
    self.assertEqual(av_read_length,99)
    self.assertEqual(num_reads_input,6)

    print('Testing Primer FW Guide RV')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, av_read_length, num_reads_input = cl.prep_input(
                    root = 'test.py.primer_fw_guide_rv',
                    primer_seq = 'TTGCAATGAAAATAAATGTTT',
                    min_primer_length = 0,
                    guide_seqs = ['TACTAAACTGGGGGATATTA'],
                    cleavage_offset = -3,
                    fastq_r1 = os.path.join(test_inputs_folder,'primer_fw.fq'),
                    samtools_command = 'samtools',
                    genome = this_genome,
                    bowtie2_command = 'bowtie2',
                    bowtie2_genome = this_bowtie2_genome,
                    can_use_previous_analysis = False,
                    suppress_file_output = True,
                    )
    if verbose:
        print(f"{origin_seq=}")
        print(f"{cut_sites=}")
        print(f"{cut_annotations=}")
        print(f"{primer_chr=}")
        print(f"{primer_loc=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    self.assertEqual(origin_seq, 'TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAA')
    self.assertEqual(primer_chr, 'chr11_sm')
    self.assertEqual(primer_loc, 101)
    self.assertEqual(av_read_length, 99)
    self.assertEqual(num_reads_input, 6)

    print('Testing Primer RV Guide RW')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, av_read_length, num_reads_input = cl.prep_input(
                    root = 'test.py.primer_rv_guide_fw',
                    primer_seq = 'TTCCTTTGTTCCCTAAGTCCA',
                    min_primer_length = 0,
                    guide_seqs = ['CTCAAGGCCCTTCATAATAT'],
                    cleavage_offset = -3,
                    fastq_r1 = os.path.join(test_inputs_folder,'primer_fw.fq'),
                    samtools_command = 'samtools',
                    genome = this_genome,
                    bowtie2_command = 'bowtie2',
                    bowtie2_genome = this_bowtie2_genome,
                    can_use_previous_analysis = False,
                    suppress_file_output = True,
                    )
    if verbose:
        print(f"{origin_seq=}")
        print(f"{cut_sites=}")
        print(f"{cut_annotations=}")
        print(f"{primer_chr=}")
        print(f"{primer_loc=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    self.assertEqual(origin_seq, 'TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATA')
    self.assertEqual(primer_chr, 'chr11_sm')
    self.assertEqual(primer_loc, 201)
    self.assertEqual(av_read_length, 99)
    self.assertEqual(num_reads_input, 6)

    print('Testing Primer RV Guide RV')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, av_read_length, num_reads_input = cl.prep_input(
                    root = 'test.py.primer_rv_guide_rv',
                    primer_seq = 'TTCCTTTGTTCCCTAAGTCCA',
                    min_primer_length = 0,
                    guide_seqs = ['TACTAAACTGGGGGATATTA'],
                    cleavage_offset = -3,
                    fastq_r1 = os.path.join(test_inputs_folder,'primer_fw.fq'),
                    samtools_command = 'samtools',
                    genome = this_genome,
                    bowtie2_command = 'bowtie2',
                    bowtie2_genome = this_bowtie2_genome,
                    can_use_previous_analysis = False,
                    suppress_file_output = True,
                    )
    if verbose:
        print(f"{origin_seq=}")
        print(f"{cut_sites=}")
        print(f"{cut_annotations=}")
        print(f"{primer_chr=}")
        print(f"{primer_loc=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    self.assertEqual(origin_seq, 'TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATA')
    self.assertEqual(primer_chr, 'chr11_sm')
    self.assertEqual(primer_loc, 201)
    self.assertEqual(av_read_length, 99)
    self.assertEqual(num_reads_input, 6)

    print('Testing additional cut site file')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, av_read_length, num_reads_input = cl.prep_input(
                    root = 'test.py.primer_fw_guide_fw',
                    primer_seq = 'TTGCAATGAAAATAAATGTTT',
                    min_primer_length = 0,
                    guide_seqs = ['CTCAAGGCCCTTCATAATAT'],
                    cleavage_offset = -3,
                    fastq_r1 = os.path.join(test_inputs_folder,'primer_fw.fq'),
                    samtools_command = 'samtools',
                    genome = this_genome,
                    bowtie2_command = 'bowtie2',
                    bowtie2_genome = this_bowtie2_genome,
                    additional_cut_site_file=os.path.join(test_inputs_folder,'sample_cuts.txt'),
                    can_use_previous_analysis = False,
                    suppress_file_output = True,
                    )
    if verbose:
        print(f"{origin_seq=}")
        print(f"{cut_sites=}")
        print(f"{cut_annotations=}")
        for cut in cut_sites:
            print(f'{cut}: {cut_annotations[cut]}')
        print(f"{primer_chr=}")
        print(f"{primer_loc=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    #here, we inserted a new cut site at chr11_sm:150. The old closest cut site was chr11_sm:161
    self.assertEqual(origin_seq, 'TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAG')
    self.assertEqual(primer_chr, 'chr11_sm')
    self.assertEqual(primer_loc, 101)
    self.assertTrue(cut in cut_sites for cut in ['chr11_sm:500', 'chr11_sm:700', 'chr11_sm:150', 'chr1:120', 'chr2:300', 'chr11_sm:161'])
    self.assertEqual(len(cut_sites), 6)
    # assert that new site is now the origin cut site
    self.assertEqual(cut_annotations['chr11_sm:150'], ['newprimer', 'Origin:left', 'RC', 'ATTCCG'])
    # assert that old site is no longer the origin cut site
    self.assertEqual(cut_annotations['chr11_sm:161'], ['Programmed', 'Not-origin', 'FW', 'CTCAAGGCCCTTCATAATAT'])

    print('TESTS PASSED')

def assert_collapse_cuts(self, verbose):
    print('Testing cut collapsing')
    novel_cut_merge_distance = 10
    known_cut_merge_distance = 15
    origin_cut_merge_distance = 50
    cut_points_by_chr = {
            'chr1':{100:1,
                    101:3,
                    102:5,
                    112:5,
                    123:5,
                    133:5,
                    144:5
                    },
            'chr2':{49:1,
                    50:1,
                    80:1,
                    90:1,
                    100:1, #origin
                    101:3,
                    102:5,
                    110:1,
                    111:1,
                    115:1,
                    120:1, #programmed not-origin
                    125:1,
                    135:1,
                    136:1,
                    140:1,
                    150:1,
                    280:1,
                    290:1,
                    300:1, #programmed not-origin
                    310:1,
                    320:1,
                    330:1
                    },
            'chr3':{100:10,
                    110:3,
                    180:1,
                    190:1,
                    200:1, #homology
                    210:1,
                    220:1},
            'chr4':{100:10},
            'chr5':{100:10},#homology
            }
    cut_annotations = {
            'chr2:100':['Programmed','Origin:left','FW','AATTCG'],
            'chr2:120':['Programmed','Not-origin','FW','AATTCG'],
            'chr2:300':['Programmed','Not-origin','FW','AATTCG']
            }
    cut_point_homology_info = {
            'chr3:200':['Programmed','Not-origin','FW','AATTCG'],
            'chr5:100':['Programmed','Not-origin','FW','AATTCG']
            }
    origin_chr, origin_cut_pos, origin_direction, final_cut_point_lookup, final_cut_points_by_chr, found_cut_point_total = cl.collapse_cut_points(
        novel_cut_merge_distance = novel_cut_merge_distance,
        known_cut_merge_distance = known_cut_merge_distance,
        origin_cut_merge_distance = origin_cut_merge_distance,
        cut_points_by_chr = cut_points_by_chr,
        cut_annotations = cut_annotations,
        cut_point_homology_info = cut_point_homology_info
        )

    if verbose:
        print(f'{origin_chr=} {origin_cut_pos=} {origin_direction=}')
        print('Final cut point lookup: ' + str(final_cut_point_lookup))
        print('Final cut points by chr: ' + str(final_cut_points_by_chr))
        print(f'{found_cut_point_total=}')

    def assert_pair(chr1,pos1,pos2):
        self.assertEqual(final_cut_point_lookup[chr1+":"+str(pos1)], chr1+":"+str(pos2))

    assert_pair('chr1',100,101) # |
    assert_pair('chr1',101,101) # | collapsed novel (10bp)
    assert_pair('chr1',102,101) # |
    assert_pair('chr1',112,112) #   not collapsed novel (because even though it's within 10bp from 102, the mean for 100-102 is 101)
    assert_pair('chr1',123,128) # | collapsed novel (10bp)
    assert_pair('chr1',133,128) # | collapsed novel (10bp)
    assert_pair('chr1',144,144) #   not collapsed novel, tests for novel not-collapsed at end of chr

    assert_pair('chr2',49,49)   # not collapsed novel
    assert_pair('chr2',50,100)  # |
    assert_pair('chr2',80,100)  # |
    assert_pair('chr2',90,100)  # |
    assert_pair('chr2',100,100) # | origin - collapsed to origin (50bp)
    assert_pair('chr2',101,100) # |
    assert_pair('chr2',102,100) # |
    assert_pair('chr2',110,100) # |   site is closer to origin than 120
    assert_pair('chr2',111,120) # |  |
    assert_pair('chr2',115,120) # |  |
    assert_pair('chr2',120,120) # |  |  programmed not-origin -collapsed to known (15bp)
    assert_pair('chr2',125,120) # |  |
    assert_pair('chr2',135,120) # |  |
    assert_pair('chr2',136,100) # |
    assert_pair('chr2',140,100) # |
    assert_pair('chr2',150,100) # |
    assert_pair('chr2',280,280) # not collapsed novel
    assert_pair('chr2',290,300) # |
    assert_pair('chr2',300,300) # |   programmed not-origin -collapsed to known (15bp)
    assert_pair('chr2',310,300) # |
    assert_pair('chr2',320,325) #   | collapsed novel, #tests for collapsing at end of chr
    assert_pair('chr2',330,325) #   | collapsed novel

    assert_pair('chr3',100,102) # | collapsed, weighted novel
    assert_pair('chr3',110,102) # |
    assert_pair('chr3',180,180) # not collapsed novel
    assert_pair('chr3',190,200) # |
    assert_pair('chr3',200,200) # | homology - collapsed to known (15bp)
    assert_pair('chr3',210,200) # |
    assert_pair('chr3',220,220) # not collapsed novel

    assert_pair('chr4',100,100) # not collapsed novel, tests for single novel site
    assert_pair('chr5',100,100) # homology, tests for single homology site

    print('TESTS PASSED')

def assert_no_results(self, this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose):
    file_roots_to_compare = [
            '.final.classifications.txt',
            '.final.final_assignments.txt',
            '.final.final_cut_points.txt',
            '.final.fragment_translocation_list.txt',
            '.final.discarded_reads.txt',
            ]
    print('Running no results (tests program when no results)')
    this_root = 'no_results'
    delete_files_with_prefix(f'{temp_outputs_folder}/{this_root}')
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/badReads.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('Running no results from previous cached run')
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/badReads.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate'.split(" "))
    self.assertTrue(settings['can_use_previous_analysis'])
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('TESTS PASSED')

def assert_dedup(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose):
    file_roots_to_compare = [
            '.final.classifications.txt',
            '.final.final_assignments.txt',
            '.final.final_cut_points.txt',
            '.final.fragment_translocation_list.txt',
            ]

    print('Checking initial deduplication')
    this_root = 'initial_dedup'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/dedup.r1.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate --fastq_umi {test_inputs_folder}/dedup_UMIs.fq'.split(" "))

    umi_r1, umi_r2 = cl.add_umi_from_umi_file(
            root=settings['root']+'.addUMI',
            fastq_r1=settings['fastq_r1'],
            fastq_r2=settings['fastq_r2'],
            fastq_umi=settings['fastq_umi'],
            can_use_previous_analysis=False,
        )

    (experiment_had_UMIs, analyze_UMI_r1, analyze_UMI_r2, post_UMI_regex_count, post_initial_dedup_count
            ) = cl.analyze_UMIs_initial(
                    root=settings['root']+'.analyzeInputUMI',
                    fastq_r1=umi_r1,
                    fastq_r2=umi_r2,
                    umi_regex=None,
                    dedup_input_on_UMI=False,
                    can_use_previous_analysis=False,
        )
    (analyze_UMI_tot_read_count, post_analyze_UMI_read_count) = parse_UMI_stats_file(self, settings['root']+'.analyzeInputUMI')
    if verbose:
        print(f'{analyze_UMI_tot_read_count=}')
        print(f'{post_UMI_regex_count=}')
        print(f'{post_initial_dedup_count=}')
        print(f'{post_analyze_UMI_read_count=}')
    self.assertEqual(analyze_UMI_tot_read_count, 8)
    self.assertEqual(post_UMI_regex_count, 8)
    self.assertEqual(post_initial_dedup_count, 8)
    self.assertEqual(post_analyze_UMI_read_count, 8)


    print('Testing with UMI regex')
    (experiment_had_UMIs, analyze_UMI_r1, analyze_UMI_r2, post_UMI_regex_count, post_initial_dedup_count
            ) = cl.analyze_UMIs_initial(
                    root=settings['root']+'.analyzeInputUMI',
                    fastq_r1=umi_r1,
                    fastq_r2=umi_r2,
                    umi_regex='AAAAA',
                    dedup_input_on_UMI=False,
                    can_use_previous_analysis=False,
        )
    (analyze_UMI_tot_read_count, post_analyze_UMI_read_count) = parse_UMI_stats_file(self,settings['root']+'.analyzeInputUMI')
    if verbose:
        print(f'{analyze_UMI_tot_read_count=}')
        print(f'{post_UMI_regex_count=}')
        print(f'{post_initial_dedup_count=}')
        print(f'{post_analyze_UMI_read_count=}')
    self.assertEqual(analyze_UMI_tot_read_count, 8)
    self.assertEqual(post_UMI_regex_count, 7)
    self.assertEqual(post_initial_dedup_count, 7)
    self.assertEqual(post_analyze_UMI_read_count, 7)

    print('Testing dedup')
    (experiment_had_UMIs, analyze_UMI_r1, analyze_UMI_r2, post_UMI_regex_count, post_initial_dedup_count
            ) = cl.analyze_UMIs_initial(
                    root=settings['root']+'.analyzeInputUMI',
                    fastq_r1=umi_r1,
                    fastq_r2=umi_r2,
                    umi_regex=None,
                    dedup_input_on_UMI=True,
                    can_use_previous_analysis=False,
        )
    (analyze_UMI_tot_read_count, post_analyze_UMI_read_count) = parse_UMI_stats_file(self,settings['root']+'.analyzeInputUMI')
    if verbose:
        print(f'{analyze_UMI_tot_read_count=}')
        print(f'{post_UMI_regex_count=}')
        print(f'{post_initial_dedup_count=}')
        print(f'{post_analyze_UMI_read_count=}')
    self.assertEqual(analyze_UMI_tot_read_count, 8)
    self.assertEqual(post_UMI_regex_count, 8)
    self.assertEqual(post_initial_dedup_count, 2) #2 UMIs seen (WT; one bad, one good)
    self.assertEqual(post_analyze_UMI_read_count, 3) #2 of WT and 1 bad UMI = 3 reads

    print('Testing UMI regex and dedup')
    (experiment_had_UMIs, analyze_UMI_r1, analyze_UMI_r2, post_UMI_regex_count, post_initial_dedup_count
            ) = cl.analyze_UMIs_initial(
                    root=settings['root']+'.analyzeInputUMI',
                    fastq_r1=umi_r1,
                    fastq_r2=umi_r2,
                    umi_regex='AAAAA',
                    dedup_input_on_UMI=True,
                    can_use_previous_analysis=False,
        )
    (analyze_UMI_tot_read_count, post_analyze_UMI_read_count) = parse_UMI_stats_file(self,settings['root']+'.analyzeInputUMI')
    if verbose:
        print(f'{analyze_UMI_tot_read_count=}')
        print(f'{post_UMI_regex_count=}')
        print(f'{post_initial_dedup_count=}')
        print(f'{post_analyze_UMI_read_count=}')
    self.assertEqual(analyze_UMI_tot_read_count, 8)
    self.assertEqual(post_UMI_regex_count, 7)
    self.assertEqual(post_initial_dedup_count, 1) # 1 UMI/WT read combination with
    self.assertEqual(post_analyze_UMI_read_count, 2) # 2 reads

    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/dedup.r1.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate --fastq_umi {test_inputs_folder}/dedup_UMIs.fq --umi_regex AAAAA --dedup_input_on_UMI'.split(" "))
    cl.processCRISPRlungo(settings)

    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('Testing UMI regex and dedup for paired sample')
    this_root = 'initial_dedup_paired'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/dedup.r1.fq --fastq_r2 {test_inputs_folder}/dedup.r2.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate --fastq_umi {test_inputs_folder}/dedup_UMIs.fq --umi_regex AAAAA --dedup_input_on_UMI'.split(" "))
    cl.processCRISPRlungo(settings)

    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    #note here that @origin_ins1 is deduplicated with WT because the origin has a 1bp insertion like this:
    #@WT1_highQual
    #TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGA
    #TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAA <origin
    #@origin_ins1
    #TTGCAATGAAAATAAATGTTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGG
    #TTGCAATGAAAATAAATGTTTTTT ATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAA <origin
    # and because we already trimmed the origin it doesn't know that there was any difference in this read. I suppose this is ok because any errors in the origin region don't have 'biological' significance i.e. they don't arise from CRISPR editing like indels at the end of the primer sequence (which would count as different reads for deduplication)
    print('Testing post dedup for paired sample')
    this_root = 'post_dedup_paired'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/dedup.r1.fq --fastq_r2 {test_inputs_folder}/dedup.r2.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate --fastq_umi {test_inputs_folder}/dedup_UMIs.fq --umi_regex AAAAA --write_discarded_read_info'.split(" "))
    cl.processCRISPRlungo(settings)

    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    reads_per_umi_by_classification_root = ".final.deduplication_by_UMI_and_aln_pos.reads_per_umi_by_classification.txt"
    self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{reads_per_umi_by_classification_root}',f'{expected_outputs_folder}/{this_root}{reads_per_umi_by_classification_root}',debug=verbose))

    print('Testing post dedup on umi and final cut for paired sample')
    this_root = 'post_dedup_cut_paired'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/dedup.r1.fq --fastq_r2 {test_inputs_folder}/dedup.r2.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate --fastq_umi {test_inputs_folder}/dedup_UMIs.fq --umi_regex AAAAA --dedup_by_final_cut_assignment_and_UMI --write_discarded_read_info'.split(" "))
    cl.processCRISPRlungo(settings)

    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    reads_per_umi_by_classification_root = ".final.deduplication_by_UMI_and_final_cut_assignment.reads_per_umi_by_classification.txt"
    self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{reads_per_umi_by_classification_root}',f'{expected_outputs_folder}/{this_root}{reads_per_umi_by_classification_root}',debug=verbose))

def parse_UMI_stats_file(self,root):
    analyze_UMI_stats_file = root+".info"
    with open(analyze_UMI_stats_file,'r') as fin:
        head_line = fin.readline()
        while head_line.startswith('#'):
            head_line = fin.readline()
        line_els = fin.readline().rstrip('\n').split("\t")
        if len(line_els) > 3:
            self.assertEqual(len(line_els), 9)
            (experiment_had_UMIs_str,fastq_analyze_UMI_r1_str,fastq_analyze_UMI_r2_str,tot_read_count_str,count_with_regex_str,post_analyze_UMI_count_str,post_analyze_UMI_read_count_str,UMI_GINI_str,UMI_count) = line_els
            experiment_had_UMIs = True if experiment_had_UMIs_str == "True" else False
            if not experiment_had_UMIs:
                return(None,None)
            fastq_analyze_UMI_r1 = None if fastq_analyze_UMI_r1_str == "None" else fastq_analyze_UMI_r1_str
            fastq_analyze_UMI_r2 = None if fastq_analyze_UMI_r2_str == "None" else fastq_analyze_UMI_r2_str
            tot_read_count = int(tot_read_count_str)
            count_with_regex = int(count_with_regex_str)
            post_analyze_UMI_count = int(post_analyze_UMI_count_str)
            post_analyze_UMI_read_count = int(post_analyze_UMI_read_count_str)
            return(tot_read_count, post_analyze_UMI_read_count)
    raise Exception('Could not parse UMI stats file ' + analyze_UMI_stats_file)

def assert_r1r2_agreement(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose):
    file_roots_to_compare = [
            '.final.classifications.txt',
            '.final.final_assignments.txt',
            '.final.final_cut_points.txt',
            '.final.fragment_translocation_list.txt',
            ]
    #Primer forward, guide forward
#primer>>>>>>>>>>>>>>>                                                              guide>          ><cut
#TTGCAATGAAAATAAATGTTT                                                              CTTAGGGAACAAAGGAACCT
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAACCTTTAATAGAAATTGGACAGCAAGAAAGCGAGCTTAGTG
#
#chr21_sm
#                      guide>          ><cut
#guide:                CTTAGGGAACAAAGGAACCT
#ob3:                  CTCAGGGAAAAAAGGAACCC
#TTGGAAGCTCTGAAATTCACATCTCAGGGAAAAAAGGAACCCCGGAAAGTGTATTCCTAGAAAGGCATCCATAGTGAACAATGAGCTCATT

#

    print('Testing R1/R2 agreement')
    this_root = 'r1r2_agreement'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/r1r2_agreement.r1.fq --fastq_r2 {test_inputs_folder}/r1r2_agreement.r2.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTTAGGGAACAAAGGAACCT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate --write_discarded_read_info'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('Testing R1/R2 agreement with smaller max distance')
    this_root = 'r1r2_agreement_small'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/r1r2_agreement.r1.fq --fastq_r2 {test_inputs_folder}/r1r2_agreement.r2.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTTAGGGAACAAAGGAACCT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate --write_discarded_read_info --r1_r2_support_max_distance 100'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

def assert_runs(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose):
    file_roots_to_compare = [
            '.final.classifications.txt',
            '.final.final_assignments.txt',
            '.final.final_cut_points.txt',
            '.final.fragment_translocation_list.txt',
            ]
    #Primer forward, guide forward
#primer>>>>>>>>>>>>>>>                      guide>          ><cut
#TTGCAATGAAAATAAATGTTT                      CTCAAGGCCCTTCATAATAT
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
#reads:
#WT
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGA
#1DelLeft
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
#1DelRight
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAAATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA

##
#----|---- wt
#-  -|---- origin del1
#-+--|---- origin ins1
#--- |----- del 1 left
#----| ---- del 1 right
#-----------ins 1
#

    print('Running Primer FW Guide FW')
    this_root = 'primer_fw_guide_fw'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_fw.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('Running Primer FW Guide FW on cached run')
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_fw.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --debug --keep_intermediate'.split(" "))
    self.assertEqual(settings['can_use_previous_analysis'], True)
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))


    #Primer forward, guide reverse
#FW
#primer>>>>>>>>>>>>>>>                                      ><cut    <guide
#TTGCAATGAAAATAAATGTTT                                    TAATATCCCCCAGTTTAGTA
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
#reads:
#WT
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGA
#1DelLeft
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
#1DelRight
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAAATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
    print('Running Primer FW Guide RV')
    this_root = 'primer_fw_guide_rv'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_fw.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences TACTAAACTGGGGGATATTA --root {temp_outputs_folder}/{this_root} --keep_intermediate'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    #Primer reverse, guide forward
#(FW)
#                                           guide>          ><cut               <<<<<<<<<<<<<<<primer
#                                           CTCAAGGCCCTTCATAATAT                TGGACTTAGGGAACAAAGGAA
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
#(RC)
#primer>>>>>>>>>>>>>>>                  ><cut                                  
#TTCCTTTGTTCCCTAAGTCCA                ATATTATGAAGGGCCTTGAG                      
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA
#reads:
#WT
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCA
#1DelLeft (right bc rc)
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA
#1DelRight (left bc rc)
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATTTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA
    print('Running Primer RV Guide FW')
    this_root = 'primer_rv_guide_fw'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_rv.fq --primer_seq TTCCTTTGTTCCCTAAGTCCA --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root} --keep_intermediate --suppress_plots'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    #Primer reverse, guide forward
#(FW)
#                                                           ><cut   <guide      <<<<<<<<<<<<<<<primer
#                                                         TAATATCCCCCAGTTTAGTA  TGGACTTAGGGAACAAAGGAA
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
#(RC)
#primer>>>>>>>>>>>>>>>  guide>          ><cut                                  
#TTCCTTTGTTCCCTAAGTCCA  TACTAAACTGGGGGATATTA
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA
#reads:
#WT
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCA
#1DelLeft (right bc rc)
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA
#1DelRight (left bc rc)
#TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATTTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA

    print('Running Primer RV Guide RV')
    this_root = 'primer_rv_guide_rv'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_rv.fq --primer_seq TTCCTTTGTTCCCTAAGTCCA --guide_sequences TACTAAACTGGGGGATATTA --root {temp_outputs_folder}/{this_root} --keep_intermediate --suppress_plots'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('TESTS PASSED')

def test_lines_equal(file1,file2,debug):
    """
    Asserts that the lines in two files are present in each other - the order may be different

    Parameters:
        file1 (str): File1 to check
        file2 (str): File2 to check
        debug (boolean): Whether to print debug statements

    Returns:
        are_same: Whether two files are the same
    """
    if debug:
        print(f'Comparing {file1} to {file2}')

    count_in_1 = 0
    count_in_2 = 0
    count_in_both = 0
    lines1 = {}
    with open(file1,'r') as f1:
        for line in f1:
            count_in_1 += 1
            lines1[line] = 1

    with open(file2,'r') as f2:
        for line in f2:
            count_in_2 += 1
            if line in lines1:
                count_in_both += 1
                del lines1[line]
            else:
                if debug:
                    print ('file 2 contained line:\n===\n' + line + '===\nnot found in file1')
                return False
    if len(lines1.keys()) > 0:
        if debug:
            print ('file 1 contained line:\n===\n' + lines1.keys()[0] + '===\nnot found in file2')
        return False

    if debug:
        print(f'Read {count_in_1} lines in {file1}, {count_in_2} lines in {file2}, found {count_in_both} in both')

    return True

def assert_cuts(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose):
    file_roots_to_compare = [
            '.final.classifications.txt',
            '.final.final_assignments.txt',
            '.final.final_cut_points.txt',
            '.final.fragment_translocation_list.txt',
            ]
#>hg38_dna range=chr11:5225464-5225564 5'pad=0 3'pad=0 strand=+ repeatMasking=none
#TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAA
#
#>hg38_dna range=chr11:5225564-5225664 5'pad=0 3'pad=0 strand=+ repeatMasking=none
#CCTTTAATAGAAATTGGACAGCAAGAAAGCGAGCTTAGTGATACTTGTGGGCCAGGGCATTAGCCACACCAGCCACCACTTTCTGATAGGCAGCCTGCAC
#
#>hg38_dna range=chr11:5225664-5225764 5'pad=0 3'pad=0 strand=+ repeatMasking=none
#TGGTGGGGTGAATTCTTTGCCAAAGTGATGGGCCAGCACACAGACCAGCACGTTGCCCAGGAGCTGTGGGAGGAAGATAAGAGGTATGAACATGATTAGC
#
#guide
#CTTAGGGAACAAAGGAA|CCT
#
#primer
#TTGCAATGAAAATAAATGTTT

#
#----|---------linear
#----|  -------2d
#--  |---------2d primer
#----|       --100d
#----|        -100d4del


    print('Running Cuts - testing for casoffinder functionality')
    this_root = 'cuts'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/cuts.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTTAGGGAACAAAGGAACCT --PAM NGG --casoffinder_num_mismatches 3 --root {temp_outputs_folder}/{this_root} --keep_intermediate --debug'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('Running Cuts - testing for collapsing cut distance')
    this_root = 'cuts_collapse'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/cuts.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTTAGGGAACAAAGGAACCT --PAM NGG --casoffinder_num_mismatches 3 --root {temp_outputs_folder}/{this_root} --keep_intermediate --debug --origin_cut_merge_distance 10 --known_cut_merge_distance 10 --novel_cut_merge_distance 10'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('Running Cuts - testing for cut annotations')
    this_root = 'cuts_annotation'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/cuts.fq --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTTAGGGAACAAAGGAACCT --PAM NGG --casoffinder_num_mismatches 3 --root {temp_outputs_folder}/{this_root} --keep_intermediate --debug --origin_cut_merge_distance 10 --known_cut_merge_distance 10 --novel_cut_merge_distance 10 --cut_classification_annotations chr11_sm:301:left:My__Annotation --cut_region_annotation_file {test_inputs_folder}/cut_region_annos.txt'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))
    print('test!')

    print('TESTS PASSED')

def assert_primer(self,this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose):
    file_roots_to_compare = [
            '.final.classifications.txt',
            '.final.final_assignments.txt',
            '.final.final_cut_points.txt',
            '.final.fragment_translocation_list.txt',
            ]
    print('Running Exogenous Primer - testing for exogenous primer functionality')
    this_root = 'exogenous'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/exogenous.fq --primer_seq aagagggttgactatgattgtttggggt --guide_sequences CTTAGGGAACAAAGGAACCT --PAM NGG --root {temp_outputs_folder}/{this_root} --keep_intermediate --min_primer_length 20 --known_cut_merge_distance 104 --debug'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('Running Primer - testing for primer functionality without guide provided')
    this_root = 'exogenous_noguide'
    settings = cl.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/exogenous.fq --primer_seq aagagggttgactatgattgtttggggt --root {temp_outputs_folder}/{this_root} --keep_intermediate --min_primer_length 20 --debug'.split(" "))
    cl.processCRISPRlungo(settings)
    for file_root in file_roots_to_compare:
        self.assertTrue(test_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose))

    print('TESTS PASSED')


def delete_files_with_prefix(prefix):
    files = glob.glob(prefix+'*')
    for f in files:
        try:
            os.remove(f)
        except OSError as e:
            print("Error: %s : %s" % (f, e.strerror))

if __name__ == "__main__":
    unittest.main()
