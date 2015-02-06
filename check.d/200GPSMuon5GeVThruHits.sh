#! /bin/sh

# Enable error trapping.
set -e
set -x

INPUT=detsim-100GPSMuon5GeVThru.root
if [ ! -f $INPUT ]; then
    echo "Missing file: " ${INPUT}
    echo MISSING INPUT
    exit
fi

OUTPUT=elecsim-200GPSMuon5GeVThruHits.root
if [ -f $OUTPUT ]; then
    rm $OUTPUT
fi

ELECSIM.exe -n 1 -o $OUTPUT $INPUT

