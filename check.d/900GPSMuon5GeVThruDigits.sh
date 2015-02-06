#!/bin/bash

INPUT=elecsim-200GPSMuon5GeVThruDigits.root

# Check that the pmt digit container was created
if ! (dump-event.exe -O dump ${INPUT} | grep 'TDigitContainer.*pmt'); then
    echo FAILED Missing pmt digit container in ${INPUT}
    exit 1
fi

# Check the number of PMTs created
if [ $(dump-event.exe -O dump ${INPUT} | grep 'MC-P-' | wc -l) -lt 5 ]; then
    echo FAILED Missing pmt digits in ${INPUT}
    exit 1
fi

# Check that the drift digit container was created
if ! (dump-event.exe -O dump ${INPUT} | grep 'TDigitContainer.*drift'); then
    echo FAILED Missing drift digit container in ${INPUT}
    exit 1
fi

# Check the number of X wires created
if [ $(dump-event.exe -O dump ${INPUT} | grep 'MC-W-X-' | wc -l) -lt 10 ]; then
    echo FAILED Missing X wire digits in ${INPUT}
    exit 1
fi

# Check the number of U wires created
if [ $(dump-event.exe -O dump ${INPUT} | grep 'MC-W-U-' | wc -l) -lt 10 ]; then
    echo FAILED Missing U wire digits in ${INPUT}
    exit 1
fi

# Check the number of V wires created
if [ $(dump-event.exe -O dump ${INPUT} | grep 'MC-W-V-' | wc -l) -lt 10 ]; then
    echo FAILED Missing V wire digits in ${INPUT}
    exit 1
fi

echo SUCCESS
