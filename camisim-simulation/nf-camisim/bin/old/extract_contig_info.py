from Bio import SeqIO
import numpy as np
import pandas as pd
import sys

assembly_file = sys.argv[1]

ids = pd.DataFrame([ np.array(seq.id.split('_'))[[0,2,4,6]]
                     for seq in SeqIO.parse(assembly_file,'fasta') ],
                   columns=['V_id','start','end','length'])
ids.start = ids.start.astype(int)-1
ids.end = ids.end.astype(int)-1

suffixes = ids.groupby(ids.V_id).cumcount().astype(str)

ids['C_id'] = ids.V_id + '|' + suffixes
ids.to_csv('metadata.csv',index=False)

assembly_formatted = []
for i,seq in enumerate(SeqIO.parse(assembly_file,'fasta')):
    seq.id = ids.C_id[i]
    assembly_formatted.append(seq)
SeqIO.write(assembly_formatted,'assembly.fasta','fasta')

