"""Console script for crisprlungo."""
import argparse
import logging
import sys
from CRISPRlungo.CRISPRlungoCore import parse_settings, processCRISPRlungo
from CRISPRlungo.CRISPRlungoBatch import processBatch
from CRISPRlungo.CRISPRlungoCompare import compareCRISPRlungo


def main():
    if len(sys.argv) == 1:
        sys.argv[0] = "CRISPRlungo"
        crisprlungo()
    if sys.argv[1] == "CRISPRlungo":
        sys.argv.remove("CRISPRlungo")
        sys.argv[0] = "CRISPRlungo"
        crisprlungo()
    elif sys.argv[1] == "CRISPRlungoBatch":
        sys.argv.remove('CRISPRlungoBatch')
        sys.argv[0] = 'CRISPRlungoBatch'
        batch()
    elif sys.argv[1] == "CRISPRlungoCompare":
        sys.argv.remove('CRISPRlungoCompare')
        sys.argv[0] = 'CRISPRlungoCompare'
        compare()
    else:
        crisprlungo()



def crisprlungo():
    if len(sys.argv) == 1:
        print('CRISPRlungo: Analyzing unidirectional sequencing of genome editing\n' + \
            'usage: CRISPRlungo --fastq_r1 r1.fq --genome hg38.fa --guide_sequences ATTA --primer_seq GCTA\n' + \
            'commonly-used arguments:\n' + \
            '-h, --help            show the full list of arguments\n' + \
            '-v, --version         show program\'s version number and exit\n' + \
            '--fastq_r1            FASTQ file containing the reads to analyze (required)\n' + \
            '--genome              genome sequence file in FASTA format (required)\n' + \
            '--guide_sequences     guide sequences (may be specified multiple times)\n' + \
            '--primer_seq          primer sequence\n\n' + \
            'Alternately, these arguments may be specified in a settings file\n' + \
            '     which is a tab-separated file of setting names and values.\n'

        )
        sys.exit()
    try:
        settings = parse_settings(sys.argv)
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

def batch():
    try:
        processBatch(sys.argv)
    except Exception as e:
        if '--debug' in sys.argv:
            raise e
        else:
            print(str(e))
            sys.exit(1)

def compare():
    try:
        compareCRISPRlungo(sys.argv)
    except Exception as e:
        if '--debug' in sys.argv:
            raise e
        else:
            print(e)
            sys.exit(1)

if __name__ == "__main__":
    sys.exit(main())
