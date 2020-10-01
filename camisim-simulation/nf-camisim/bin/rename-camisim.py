#!/usr/bin/env python2

from pathlib2 import Path

for folder in Path('.').glob('*/*_sample_*'):
    sample = folder.name.split('_')[-1]

    for bam in folder.glob('bam/*.bam'):
        if '-' in bam.stem:
            continue
        bam.rename('%s/%s-sample_%s.bam' % (bam.parent, bam.stem, sample))
        Path(bam.parent, bam.stem + '.bam.bai').rename('%s/%s-sample_%s.bai' % (bam.parent, bam.stem, sample))
