import sys
import os
from collections import defaultdict
import CRISPRlungo
import multiprocessing
import multiprocessing.pool
import logging

def processBatch(arg_els):
    usage_str = 'Usage: CRISPRlungoBatch.py {batch file} {settings file} {?--batch_processes int} {additional settings}'

    #batch_file #settings_file #list of other settings
    sep = "\t"
    if len(arg_els) < 2:
        raise Exception(usage_str)

    batch_process_count = 1
    for idx,arg in enumerate(arg_els):
        if arg == '--batch_processes':
            if len(arg_els) < idx+1:
                raise Exception('batch_processes parameter must be an integer')
            try:
                batch_process_count = int(arg_els[idx+1])
            except ValueError:
                raise Exception('batch_processes parameter must be an integer (got "%s")'%arg_els[idx+1])
            arg_els.pop(idx) #pop arg name
            arg_els.pop(idx) #pop value

    batch_file = arg_els[1]

    if not os.path.isfile(batch_file):
        raise Exception('Cannot find batch file "'+batch_file+'"\n'+usage_str)


    logger = logging.getLogger('CRISPRlungoBatch')
    logging_level = logging.INFO
    if '--debug' in arg_els:
        logging_level=logging.DEBUG

    log_formatter = logging.Formatter("%(asctime)s:%(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    #logger for Batch
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging_level)
    sh.setFormatter(log_formatter)
    logger.addHandler(sh)
    logger.setLevel(logging_level)

    #logger for sub runs
    sub_logger = logging.getLogger('CRISPRlungo')
    fh = logging.FileHandler(batch_file+'.CRISPRlungoBatch.log')
    fh.setFormatter(log_formatter)
    sub_logger.addHandler(fh)
    logger.addHandler(fh)
    if '--debug' in arg_els:
        sub_logger.addHandler(sh)

    logger.info('Printing sub-command logs to: ' + batch_file+'.CRISPRlungoBatch.log')

    sub_arg_els = []
    if len(arg_els) > 2:
        sub_arg_els = arg_els[2:]

    settings_arr = []
    names_arr = []
    with open(batch_file,'r') as fin:
        head = fin.readline().strip()
        head_els = head.split(sep)
        if 'name' not in head_els and 'root' not in head_els:
            raise Exception('Run name must be provided in batch file')
        line_idx = 1
        for line in fin:
            line_els = line.strip().split(sep)

            line_idx += 1
            if line_els == ['']:
                continue

            if len(line_els) != len(head_els):
                raise Exception('Incorrect number of items on line ' + str(line_idx) + ' (got ' + str(len(line_els)) + ' expected ' + str(len(head_els)) + ')')

            this_name = ''
            this_command_args = ['CRISPRlungo']
            for idx in range(len(line_els)):
                head_val = head_els[idx]
                this_val = line_els[idx]

                if head_val == 'name' or head_val == 'root':
                    this_name = this_val

                if this_val.lower() == 'true':
                    this_command_args.extend(['--'+head_val])
                elif this_val.lower() == 'false':
                    pass
                else:
                    this_command_args.extend(['--'+head_val,this_val])

            names_arr.append(this_name)

            this_command_args.extend(sub_arg_els)
            print('this command args: ' + str(this_command_args))

            try:
                logger.debug('Parsing settings for line %d'%line_idx)
                settings = CRISPRlungo.parse_settings(this_command_args)
                settings_arr.append(settings)
            except Exception as e:
                raise Exception('Error in parsing ' + batch_file + ' line ' + str(line_idx)+ ":\n" + str(e)) from e

    logger.info('Running CRISPRlungo with %d processes'%batch_process_count)
    result_summary_files = []
    with NonDaemonPool(processes=batch_process_count) as pool:
        result_summary_files = pool.map(runOneCRISPRlungo,settings_arr)

    if len(result_summary_files) != len(names_arr):
        raise Exception('Number of results and number of input batches does not match')

    data = defaultdict(lambda: defaultdict(int))
    all_head_els = {}
    for run_idx in range(len(result_summary_files)):
        this_name = names_arr[run_idx]
        with open(result_summary_files[run_idx],'r') as fin:
            head_line = fin.readline()
            val_line = fin.readline()
            head_els = head_line.strip().split("\t")
            val_els = val_line.strip().split("\t")
            for (head_el, val_el) in zip(head_els,val_els):
                all_head_els[head_el] = 1
                data[this_name][head_el] = val_el

    column_heads = ['total_input_reads','post_dedup_reads','post_primer_filter_reads','post_quality_filter_reads','discarded_quality_filter_reads','Not primary alignment','Duplicate read','Not supported by R2','Poor alignment','analyzed_read_count']
    for column_head in column_heads:
        if column_head in all_head_els:
            del all_head_els[column_head]

    indel_heads = [x for x in all_head_els.keys() if 'indels' in x]

    for indel_head in indel_heads:
        del all_head_els[indel_head]

    classification_heads = all_head_els.keys()

    column_heads.extend(classification_heads)
    column_heads.extend(indel_heads)


    output_summary_file = batch_file+".summary.txt"
    with open(output_summary_file,'w') as fout:
        fout.write("Name\t"+'\t'.join(column_heads)+"\n")
        for run_idx in range(len(result_summary_files)):
            this_name = names_arr[run_idx]
            row_str = this_name
            for head in column_heads:
                val = data[this_name][head]
                row_str += "\t"+str(val)
            fout.write(row_str+"\n")

    logger.info('Printed ' + output_summary_file)
    logger.info('Finished')

def runOneCRISPRlungo(settings):
    logger = logging.getLogger('CRISPRlungoBatch')
    logger.info('Running ' + settings['root'])
    final_file = CRISPRlungo.processCRISPRlungo(settings)
    logger.info('Finished ' + settings['root'])
    return final_file

class NonDaemonPool(multiprocessing.pool.Pool):
    """
    Pool that is not daemonic so it can have children

    https://stackoverflow.com/questions/52948447/error-group-argument-must-be-none-for-now-in-multiprocessing-pool
    """
    def Process(self, *args, **kwds):
        proc = super(NonDaemonPool, self).Process(*args, **kwds)

        class NonDaemonProcess(proc.__class__):
            """Monkey-patch process to ensure it is never daemonized"""
            @property
            def daemon(self):
                return False

            @daemon.setter
            def daemon(self, val):
                pass

        proc.__class__ = NonDaemonProcess
        return proc


if __name__ == "__main__":
    try:
        processBatch(sys.argv)
    except Exception as e:
        if '--debug' in sys.argv:
            raise e
        else:
            print(str(e))
            sys.exit(1)
