#!/usr/bin/env bash

set -e

[ -z "$1" ] && exit 1 || data="$1"

for f in $(ls -d $data/camisim* | shuf); do
    if [[ $(basename $f) =~ [0-9]$ ]] ; then
        make sim FOLDER=$f
    elif [ $f != station-ALOHA ]; then
        make SA FOLDER=$f
    fi
done

set +e
