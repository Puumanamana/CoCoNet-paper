import argparse
from pathlib import Path

import pandas as pd
import numpy as np


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=str)
    args = parser.parse_args()

    return args

    
def main():

    args = parse_args()

    summary_files = Path(args.folder).glob('*/quality_summary.tsv')

    summaries = {}

    for filename in summary_files:
        name = filename.parent.name.split('-')[-1]
        summary = pd.read_csv(
            filename, index_col=0, sep='\t',
            usecols=['contig_id', 'contig_length', 'viral_genes', 'checkv_quality',
                     'miuvig_quality', 'completeness', 'contamination'],
            low_memory=False
        )
        summaries[name] = summary

    summaries = pd.concat(summaries)
    assembly_quality = dict(
        checkv=summary.checkv.value_counts(),
        miuvig=summary.miuvig_quality.value_counts()
    )
        
    # distr_plot(summary.viral_genes)
    # distr_plot(summary.contamination)
    # distr_plot(summary.completeness)
    
    

if __name__ == '__main__':
    main()
