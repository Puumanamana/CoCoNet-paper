#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--bins', type=str)
    parser.add_argument('--min-dtr-size', type=int, default=10)
    parser.add_argument('--max-dtr-size', type=int, default=300)
    parser.add_argument('--min-dtr-id', type=float, default=0.95)
    
    args = parser.parse_args()

    args.blast_kw = dict(min_dtr_size=args.min_dtr_size,
                         max_dtr_size=args.max_dtr_size,
                         min_dtr_id=args.min_dtr_id)
    return args

                    
def main():

    args = parse_args()

    contigs = {ctg.id: ctg for ctg in SeqIO.parse(args.fasta, 'fasta')}

    assignments = pd.read_csv(args.bins, header=None, names=['contig', 'bin_id'])
    assignments.bin_id = pd.factorize(assignments.bin_id)[0]
    assignments = assignments.set_index('contig').bin_id.reindex(contigs.keys())
    # Add singletons
    assignments[assignments.isnull()] = 1 + assignments.max() + range(assignments.isnull().sum())
    assignments = assignments.astype(int)
    
    prefix = Path(args.bins).stem
    order_bins(assignments, contigs, output=f'{prefix}-merged.fasta', **args.blast_kw)
    assignments.to_csv(f'{prefix}-complete.csv')

def order_bins(assignments, contigs, output=None, **blast_kw):

    handle = open(output, 'w')
    grouped = assignments.reset_index().groupby('bin_id').contig.agg(list)
    
    for (bin_id, ctg_names) in grouped.iteritems():

        contigs_in_bin = [contigs[name] for name in ctg_names]
        if len(ctg_names) > 1:
            ends_fa = cut_ends(contigs_in_bin, max_dtr_size=blast_kw['max_dtr_size'])
            blast_xml = self_blast(ends_fa)

            matches = filter_blast(blast_xml, **blast_kw)

            if matches is not None:
                print(f'Bin {bin_id}: DTR found for {matches[0]} and {matches[1]}')
                first = next(ctg for ctg in contigs_in_bin if ctg.id == matches[0])
                last = next(ctg for ctg in contigs_in_bin if ctg.id == matches[1])
                contigs_in_bin = [ctg for ctg in contigs_in_bin if ctg.id not in matches]
                contigs_in_bin = [first] + contigs_in_bin + [last]
                                   
        (name, seq) = concatenate_contigs(contigs_in_bin, bin_id=bin_id)
        handle.write(f'>{name}\n{seq}\n')

    handle.close()

def cut_ends(contigs, max_dtr_size=300, output='extr.fasta'):
    extremities = []

    for ctg in contigs:
        start = SeqRecord(
            id=f'{ctg.id}@start',
            description='',
            seq=ctg.seq[:max_dtr_size]
        )
        end = SeqRecord(
            id=f'{ctg.id}@end',
            description='',
            seq=ctg.seq[-max_dtr_size:]
        )
        extremities += [start, end]

    SeqIO.write(extremities, output, 'fasta')

    return output

def self_blast(fasta, output='blastn.xml'):
    make_db = NcbimakeblastdbCommandline(
        dbtype="nucl", input_file=fasta, out='db'
    )
    blastn_on_db = NcbiblastnCommandline(
        query=fasta, db='db', out=output, outfmt=5,
        task='blastn-short'
    )

    make_db()
    blastn_on_db()

    # Cleaning
    for f in Path('.').glob('db*'):
        f.unlink()

    Path(fasta).unlink()
    
    return output

def filter_blast(blast_xml, min_dtr_id=0.95, min_dtr_size=10, max_dtr_size=300):

    for record in NCBIXML.parse(open(blast_xml)):
        # 1 record is one blast result from the query sequences
        if not record.alignments:
            continue
        
        (query_name, query_pos) = record.query.split('@')

        for hit in record.alignments:
            # 1 alignment is one hit (with potentially multiple HSPs)
            (hit_name, hit_pos) = hit.hit_def.split('@')

            # We need different contigs and a (start, end) pair
            if hit_name == query_name or hit_pos == query_pos:
                continue

            for hsp in hit.hsps:
                # 1) check if match is long enough and with enough identity
                if (hsp.identities / hsp.align_length >= min_dtr_id and
                    hsp.align_length >= min_dtr_size):
                    
                    # 2) check if it's at the extremities and return in correct order
                    if (query_pos == 'end' and
                        hsp.sbjct_start == 1 and hsp.query_end == max_dtr_size):
                        return [hit_name, query_name]
                    
                    elif (query_pos == 'start' and
                          hsp.query_start == 1 and hsp.sbjct_end == max_dtr_size):
                        return (query_name, hit_name)

def concatenate_contigs(contigs, bin_id=0):
    return (
        f'bin_{bin_id} size={len(contigs)}',
        ''.join(str(ctg.seq) for ctg in contigs)
    )
            
if __name__ == '__main__':
    main()
