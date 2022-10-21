import sys
import os
sys.path.append(os.path.abspath("../"))
import CRISPRlungo

this_bowtie2_genome = 'genome.chr11_5225364_5225663/genome'
this_genome = this_bowtie2_genome + '.fa'
test_inputs_folder = 'testInputs'
expected_outputs_folder = 'expectedResults'
temp_outputs_folder = 'tmpOutputs'

if not os.path.exists(this_genome):
    raise Exception(f'Missing genome file at {this_genome}')
if not os.path.exists(test_inputs_folder):
    raise Exception(f'Missing test inputs folder at {test_inputs_folder}')
if not os.path.exists(expected_outputs_folder):
    raise Exception(f'Missing expected outputs folder at {expected_outputs_folder}')
if not os.path.exists(temp_outputs_folder):
    os.makedirs(temp_outputs_folder)

verbose = False

def main():
    assert_prep_input(this_genome,this_bowtie2_genome,verbose)
    assert_runs(this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose)

def assert_prep_input(this_genome, this_bowtie2_genome, verbose):
    #test primer fw, guide fw
    print('Testing Primer FW Guide FW')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, primer_is_genomic, av_read_length, num_reads_input = CRISPRlungo.prep_input(
                    root = 'test.py.primer_fw_guide_fw',
                    primer_seq = 'TTGCAATGAAAATAAATGTTT',
                    guide_seqs = ['CTCAAGGCCCTTCATAATAT'],
                    cleavage_offset = -3,
                    fastq_r1 = 'testOrientations/primer_fw.fa',
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
        print(f"{primer_is_genomic=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    assert(origin_seq == 'TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAA')
    assert(primer_chr == 'chr11')
    assert(primer_loc == 5225464)
    assert(primer_is_genomic == True)
    assert(av_read_length == 99)
    assert(num_reads_input == 4)

    print('Testing Primer FW Guide RV')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, primer_is_genomic, av_read_length, num_reads_input = CRISPRlungo.prep_input(
                    root = 'test.py.primer_fw_guide_rv',
                    primer_seq = 'TTGCAATGAAAATAAATGTTT',
                    guide_seqs = ['TACTAAACTGGGGGATATTA'],
                    cleavage_offset = -3,
                    fastq_r1 = 'testOrientations/primer_fw.fa',
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
        print(f"{primer_is_genomic=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    assert(origin_seq == 'TTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAA')
    assert(primer_chr == 'chr11')
    assert(primer_loc == 5225464)
    assert(primer_is_genomic == True)
    assert(av_read_length == 99)
    assert(num_reads_input == 4)

    print('Testing Primer RV Guide RW')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, primer_is_genomic, av_read_length, num_reads_input = CRISPRlungo.prep_input(
                    root = 'test.py.primer_rv_guide_fw',
                    primer_seq = 'TTCCTTTGTTCCCTAAGTCCA',
                    guide_seqs = ['CTCAAGGCCCTTCATAATAT'],
                    cleavage_offset = -3,
                    fastq_r1 = 'testOrientations/primer_fw.fa',
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
        print(f"{primer_is_genomic=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    assert(origin_seq == 'TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATA')
    assert(primer_chr == 'chr11')
    assert(primer_loc == 5225564)
    assert(primer_is_genomic == True)
    assert(av_read_length == 99)
    assert(num_reads_input == 4)

    print('Testing Primer RV Guide RV')
    origin_seq, cut_sites, cut_annotations, primer_chr, primer_loc, primer_is_genomic, av_read_length, num_reads_input = CRISPRlungo.prep_input(
                    root = 'test.py.primer_rv_guide_rv',
                    primer_seq = 'TTCCTTTGTTCCCTAAGTCCA',
                    guide_seqs = ['TACTAAACTGGGGGATATTA'],
                    cleavage_offset = -3,
                    fastq_r1 = 'testOrientations/primer_fw.fa',
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
        print(f"{primer_is_genomic=}")
        print(f"{av_read_length=}")
        print(f"{num_reads_input=}")
    assert(origin_seq == 'TTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATA')
    assert(primer_chr == 'chr11')
    assert(primer_loc == 5225564)
    assert(primer_is_genomic == True)
    assert(av_read_length == 99)
    assert(num_reads_input == 4)

def assert_runs(this_genome,this_bowtie2_genome,test_inputs_folder,temp_outputs_folder,expected_outputs_folder,verbose):
    file_roots_to_compare = [
            '.final.classifications.txt',
            '.final.final_assignments',
            '.final.final_cut_points.txt',
            '.fragment_translocation_list.txt',
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
    print('Running Primer FW Guide FW')
    this_root = 'primer_fw_guide_fw'
    settings,logger = CRISPRlungo.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_fw.fa --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root}'.split(" "))
    CRISPRlungo.processCRISPRlungo(settings,logger)
    for file_root in file_roots_to_compare:
        assert_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose)


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
    settings,logger = CRISPRlungo.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_fw.fa --primer_seq TTGCAATGAAAATAAATGTTT --guide_sequences TACTAAACTGGGGGATATTA --root {temp_outputs_folder}/{this_root}'.split(" "))
    CRISPRlungo.processCRISPRlungo(settings,logger)
    for file_root in file_roots_to_compare:
        assert_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose)

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
    settings,logger = CRISPRlungo.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_rv.fa --primer_seq TTCCTTTGTTCCCTAAGTCCA --guide_sequences CTCAAGGCCCTTCATAATAT --root {temp_outputs_folder}/{this_root}'.split(" "))
    CRISPRlungo.processCRISPRlungo(settings,logger)
    for file_root in file_roots_to_compare:
        assert_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose)

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
    settings,logger = CRISPRlungo.parse_settings(f'scriptname --genome {this_genome} --fastq_r1 {test_inputs_folder}/primer_rv.fa --primer_seq TTCCTTTGTTCCCTAAGTCCA --guide_sequences TACTAAACTGGGGGATATTA --root {temp_outputs_folder}/{this_root}'.split(" "))
    CRISPRlungo.processCRISPRlungo(settings,logger)
    for file_root in file_roots_to_compare:
        assert_lines_equal(f'{temp_outputs_folder}/{this_root}{file_root}',f'{expected_outputs_folder}/{this_root}{file_root}',debug=verbose)

def assert_lines_equal(file1,file2,debug):
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

if __name__ == "__main__":
    main()
