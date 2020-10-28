#!/usr/bin/env bash

set -e

[ -z "$1" ] && exit 1 || data="$1"

for f in $(ls -d $data/* | shuf); do
    if [[ $(basename f) =~ ^[0-9] ]] ; then
        h5=$f/coverage_contigs.h5
        [ -f "${h5}.gz" ] && is_gz=1 || is_gz=0
        [ is_gz -eq 1 ] && unpigz -p 20 "${h5}.gz"
        make sim FOLDER=$f
        [ is_gz -eq 1 ] && pigz -p 20 $h5
    elif [ $f != Station_Aloha-metaspades ]; then
        make SA FOLDER=$f
    fi
done

set +e
