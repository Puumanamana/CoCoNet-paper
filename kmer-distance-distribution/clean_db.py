from random import choice
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

AMBIGUOUS_CODES = {
    'R': 'AG',
    'Y': 'CT',
    'S': 'CG',
    'W': 'AT',
    'K': 'GT',
    'M': 'AT',
    'B': 'CGT',
    'V': 'ACG',
    'D': 'AGT',
    'H': 'ACT',
    'N': 'ACGT',
}

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--output', type=str)
    parser.add_argument('--min-length', type=int, default=3000)
    args = parser.parse_args()

    return args


def clean(fasta, output, min_length=None):

    n_sequences = sum(1 for line in open(fasta) if line.startswith('>'))

    writer = open(output, 'w')
    with open(fasta, 'r') as handle:
        for i, (title, seq) in enumerate(SimpleFastaParser(handle)):
            clean_seq = ''.join([nucl if nucl in 'ACGT' else choice(AMBIGUOUS_CODES[nucl])
                                 for nucl in seq])
            if len(clean_seq) > min_length:
                writer.write(f'>{title}\n{clean_seq}\n')

            print(f'{i:,}/{n_sequences:,}', end='\r')

    writer.close()

def main():
    args = parse_args()
    clean(args.fasta, args.output, min_length=args.min_length)


if __name__ == '__main__':
    main()
