#!/usr/bin/env bash

set -e

for f in `ls input_data | shuf `; do
    if [ -f results/test-${f}.csv ]; then
        echo "$f already processed"
        continue
    elif [[ $f =~ ^[0-9] ]] ; then
        h5=input_data/$f/coverage_contigs.h5
        [ -f "${h5}.gz" ] \
            && unpigz -p 20 "${h5}.gz"
        make sim FOLDER=$f
        [ -f "${h5}" ] && pigz -p 20 $h5
    elif [ $f != Station_Aloha-metaspades ]; then
        make SA FOLDER=$f
    fi
done

set +e
