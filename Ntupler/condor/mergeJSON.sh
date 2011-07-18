#!/bin/bash

infile=$1

outfile=`sed -n '1 s/\.root/\.json/gp' <$infile`
firstfile=`sed -n '2 s/\.root/\.json/gp' <$infile`

cp $firstfile $outfile

for json in `sed -n '3,$ s/\.root/\.json/gp' <$infile`; do
    ./combine_JSON.py -r or -o $outfile -a $outfile -b $json
done
