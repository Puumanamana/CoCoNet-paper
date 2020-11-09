#!/usr/bin/env python

import argparse
import pandas as pd
import re


def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    args = parser.parse_args()

    return args

def main():
    args = parse_args()

    ctg_info = get_meta(args.fasta)
    rename_fasta(args.fasta, ctg_info)

def get_meta(fasta, output='metadata.csv'):

    pattern = re.compile(r'>(.*)_from_([0-9]+)_to_([0-9]+)_total_([0-9]+)')
    ctg_info = []
    for line in open(fasta):
        if line.startswith('>'):
            (name, start, end, size) = re.findall(pattern, line)[0]
            ctg_info.append([name, int(start), int(end), int(size)])

    ctg_info = pd.DataFrame(ctg_info, columns=['V_id', 'start', 'end', 'size'])

    suffixes = ctg_info.groupby('V_id').cumcount().astype(str)
    ctg_info['C_id'] = ctg_info.V_id + '|' + suffixes

    ctg_info.to_csv(output, index=False)

    return ctg_info
    

def rename_fasta(fasta, metadata, output='assembly.fasta'):
    i = 0
    writer = open(output, 'w')
    for line in open(fasta):
         if not line.startswith('>'):
             writer.write(line)
         else:
             writer.write('>'+metadata.C_id[i]+'\n')
             i += 1
             
    writer.close()
    
if __name__ == '__main__':
    main()
