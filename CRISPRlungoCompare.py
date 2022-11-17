import sys
import os
import argparse
from collections import defaultdict
import CRISPRlungo
import logging
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle,Patch
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def compareCRISPRlungo(arg_els):

    logger = logging.getLogger('CRISPRlungoCompare')
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

    parser = argparse.ArgumentParser(description='CRISPRlungoCompare: Compares two runs from CRISPRlungo', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder1_root','-f1', type=str, help='Root of the first completed CRISPRlungo analysis', required=True)
    parser.add_argument('--folder1_name','-n1', type=str, help='Name of the first completed CRISPRlungo analysis')
    parser.add_argument('--folder2_root','-f2', type=str, help='Root of the second completed CRISPRlungo analysis', required=True)
    parser.add_argument('--folder2_name','-nl', type=str, help='Name of the second completed CRISPRlungo analysis')
    parser.add_argument('--output_root','-o', type=str, help='Root of the second completed CRISPRlungo analysis', required=True)
    parser.add_argument('--debug', help='Print debug messages', action='store_true')

    args = parser.parse_args(arg_els[1:])

    f1_tx_list_report = args.folder1_root + '.final.fragment_translocation_list.txt'
    check_file(f1_tx_list_report)
    f2_tx_list_report = args.folder2_root + '.final.fragment_translocation_list.txt'
    check_file(f2_tx_list_report)
    f1_final_cut_counts,f1_cut_annotations,f1_cut_classification_lookup = parse_frag_tx_list(f1_tx_list_report)
    f2_final_cut_counts,f2_cut_annotations,f2_cut_classification_lookup = parse_frag_tx_list(f2_tx_list_report)

    f1_name = "Folder 1" if args.folder1_name is None else args.folder1_name
    f1_name = f1_name.replace(" ","_").replace("/","_")
    f2_name = "Folder 2" if args.folder2_name is None else args.folder2_name
    f2_name = f2_name.replace(" ","_").replace("/","_")

    plot_root = args.output_root+"."+f1_name+"_vs_"+f2_name
    plot_tx_list(f1_final_cut_counts,f1_cut_annotations,f1_cut_classification_lookup, f2_final_cut_counts, f2_cut_annotations, f2_cut_classification_lookup, plot_root, f1_name=f1_name, f2_name=f2_name)


    plot_root = args.output_root+"."+f2_name+"_vs_"+f1_name
    plot_tx_list(f2_final_cut_counts,f2_cut_annotations,f2_cut_classification_lookup, f1_final_cut_counts, f1_cut_annotations, f1_cut_classification_lookup, plot_root, f1_name=f2_name, f2_name=f1_name)

    f1_summary_file = args.folder1_root + '.summary.txt'
    check_file(f1_summary_file)
    f2_summary_file = args.folder2_root + '.summary.txt'
    check_file(f2_summary_file)





    result_summary_files = [f1_summary_file, f2_summary_file]
    names_arr = [f1_name,f2_name]

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


    output_summary_file = args.output_root+".summary.txt"
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

def parse_frag_tx_list(filename):
    """
    Parses a fragment translocation list file from a completed CRISPRlungo run

    Args:
        filename: Filename from which to parse translocations. Should have headers "fragment, cut_annotation, fragment_annotation, count" and has the suffix ".final.fragment_translocation_list.txt"

    Returns:
        final_cut_counts: dict of cut>count (e.g. chr1:234:left > 5)
        cut_annotations: dict of cut>annotation (e.g. chr1:234:left > 'Novel') Note that in CRISPRlungo the key for this dict is only 'chr1:234' instead of 'chr1:234:left'
        cut_classification_lookup: dict of cutID=> type (e.g. Linear, Translocation, etc)
    """

    logger = logging.getLogger('CRISPRlungoCompare')
    tx_keys = []
    final_cut_counts = {}
    cut_annotations = {}
    cut_classification_lookup = {}
    read_count = 0
    with open(filename,'r') as fin:
        head = fin.readline()
        for line in fin:
            read_count += 1
            (tx_key,cut_annotation,tx_classification,tx_count) = line.split("\t")
            tx_cut = ':'.join(tx_key.split(":")[0:2])
            final_cut_counts[tx_key] = int(tx_count)
            cut_annotations[tx_key] = cut_annotation
            cut_classification_lookup[tx_key] = tx_classification

    logger.debug('Read %d lines from %s'%(read_count,filename))
    return(final_cut_counts,cut_annotations,cut_classification_lookup)


def check_file(filename):
    """
    Check to see if file exists. If not, raise exception.
    
    Args:
        filename: str, path of file to check
    """
    if not os.path.isfile(filename):
        raise Exception('File ' + str(filename) + ' does not exist')

def plot_tx_list(f1_final_cut_counts,f1_cut_annotations,f1_cut_classification_lookup, 
        f2_final_cut_counts, f2_cut_annotations,f2_cut_classification_lookup, 
        plot_name, f1_name=None, f2_name=None):
    """
    Plots comparison of translocation counts and writes summary

    requires two dicts for each finished CRISPRlungo run:
        final_cut_counts: dict of cut>count (e.g. chr1:234:left > 5)
        cut_annotations: dict of cut>annotation (e.g. chr1:234:left > 'Novel') Note that in CRISPRlungo the key for this dict is only 'chr1:234' instead of 'chr1:234:left'
        cut_classification_lookup: dict of cutID=> type (e.g. Linear, Translocation, etc) (e.g. chr1:234:left > 'Linear')

    Args:
        f1_final_cut_counts: final_cut_counts for folder 1
        f1_cut_annotations: cut_annotations for folder 1
        f1_cut_classification_lookup: cut_classification_lookup for folder 1
        f2_final_cut_counts: final_cut_counts for folder 2
        f2_cut_annotations: cut_annotations for folder 2
        f2_cut_classification_lookup: cut_classification_lookup for folder 2

    Returns:
        Nothin
"""
    logger = logging.getLogger('CRISPRlungoCompare')

    sorted_tx_list = sorted(f1_final_cut_counts, key=f1_final_cut_counts.get,reverse=True)
    top_number_to_plot = 20
    top_sorted_tx_list = sorted_tx_list[:min(top_number_to_plot,len(sorted_tx_list))]

    f1_other_tx_count = 0 #number of locations not plotted
    f1_other_tx_read_count = 0 # number of reads to other locations not plotted
    if len(sorted_tx_list) > top_number_to_plot:
        for i in range(top_number_to_plot,len(sorted_tx_list)):
            f1_other_tx_count += 1
            f1_other_tx_read_count += f1_final_cut_counts[sorted_tx_list[i]]

    f2_other_tx_count = 0 #number of locations not plotted
    f2_other_tx_read_count = 0 # number of reads to other locations not plotted
    f2_included_tx_count = 0
    f2_included_tx_read_count = 0
    for tx_order_obj in top_sorted_tx_list:
        if tx_order_obj in f2_final_cut_counts:
            f2_included_tx_count += 1
            f2_included_tx_read_count += f2_final_cut_counts[tx_order_obj]
        else:
            f2_other_tx_count += 1
            f2_other_tx_read_count += f2_final_cut_counts[tx_order_obj]

    if len(top_sorted_tx_list) > 0:
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
        counts_1 = []
        counts_2 = []
        fill_cols = []
        outline_cols = []
        for tx_order_obj in top_sorted_tx_list:
            left_labs.append(f1_cut_classification_lookup[tx_order_obj])
            right_labs.append(tx_order_obj)
            counts_1.append(f1_final_cut_counts[tx_order_obj])
            val = 0
            if tx_order_obj in f2_final_cut_counts:
                val = f2_final_cut_counts[tx_order_obj]
            counts_2.append(val)

            this_chr_pos = ':'.join(tx_order_obj.split(':')[0:2])
            this_fill_col = color_lookup[this_chr_pos]
            this_cut_anno='Novel'
            if this_chr_pos in f1_cut_annotations:
                this_cut_anno = f1_cut_annotations[this_chr_pos][0]

            fill_cols.append(this_fill_col)
            if 'Cas-OFFinder' in this_cut_anno: #casoffinder categories look like "Cas-OFFinder OB 3"
                this_cut_anno = 'Cas-OFFinder'
            this_outline_col = cut_category_lookup[this_cut_anno]
            outline_cols.append(this_outline_col)

        if f1_other_tx_count > 0 or f2_other_tx_count > 0:
            left_labs.append('Other')
            right_labs.append('('+str(f1_other_tx_count)+' vs ' + str(f2_other_tx_count) + ' locations)')
            counts_1.append(f1_other_tx_read_count)
            counts_2.append(f2_other_tx_read_count)

            fill_cols.append('0.8') #light gray
            outline_cols.append('None')

        legend_outline_cols=cut_colors.colors
        legend_outline_labs=cut_categories

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

        fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharey=True, figsize=(12,8))

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

        f1_rects = []
        f2_rects = []
        for ind in range(num_boxes):
            val = max(1,counts_1[ind]) #min value is 1 or matplotlib flips out
            f1_rects.append(Rectangle((1,ys[ind]+count_bar_ydiff),val,y_height-(count_bar_ydiff*2)))
            ax2.text(x=val,y=ys[ind]+y_height/2,s=" " + str(counts_1[ind]),ha='left',va='center')

            val = max(1,counts_2[ind]) #min value is 1 or matplotlib flips out
            f2_rects.append(Rectangle((1,ys[ind]+count_bar_ydiff),val,y_height-(count_bar_ydiff*2)))
            ax3.text(x=val,y=ys[ind]+y_height/2,s=" " + str(counts_2[ind]),ha='left',va='center')

        f1_pc = PatchCollection(f1_rects)
        ax2.add_collection(f1_pc)
        ax2.axes.get_yaxis().set_visible(False)
        ax2.set_frame_on(False)

        ax2.set_xscale("log")
        ax2.set_xlim(1,max(5,max(counts_1)))
        ax2.set_xlabel('Number of reads')

        if f1_name is not None:
            ax2.set_title(f1_name)


        f2_pc = PatchCollection(f1_rects)
        ax3.add_collection(f2_pc)
        ax3.axes.get_yaxis().set_visible(False)
        ax3.set_frame_on(False)

        ax3.set_xscale("log")
        ax3.set_xlim(1,max(5,max(counts_2)))
        ax3.set_xlabel('Number of reads')

        if f2_name is not None:
            ax3.set_title(f2_name)

        plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')

    tx_list_report = plot_name + '.txt'
    with open (tx_list_report,'w') as fout:
        fout.write('fragment\tcut_annotation\tfragment_annotation\tcount_%s\tcount_%s\n'%(f1_name,f2_name))
        for tx_key in f1_final_cut_counts:
            cut_annotation = 'Novel (BOGUS'
            if tx_key in f1_cut_annotations:
                cut_annotation = f1_cut_annotations[tx_key]
            tx_classification = 'Unknown (BOGUS'
            if tx_key in f1_cut_classification_lookup:
                tx_classification = f1_cut_classification_lookup[tx_key]
            count1 = f1_final_cut_counts[tx_key]
            count2 = 0 if tx_key not in f2_final_cut_counts else f2_final_cut_counts[tx_key]
            fout.write("%s\t%s\t%s\t%s\t%s\n"%(tx_key,cut_annotation,tx_classification,count1,count2))
        logger.debug('Wrote fragment translocation list ' + tx_list_report)


if __name__ == "__main__":
    try:
        compareCRISPRlungo(sys.argv)
    except Exception as e:
        if '--debug' in sys.argv:
            raise e
        else:
            print(str(e))
            sys.exit(1)
