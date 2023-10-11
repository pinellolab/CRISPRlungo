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
