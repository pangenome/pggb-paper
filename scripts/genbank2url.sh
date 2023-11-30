#!/bin/bash

GENBANK_IDS=$1
OUTPUT_FILE=$2

cut -f 1 $GENBANK_IDS | while read -r acc ; do
  esearch -db assembly -query $acc </dev/null \
    | esummary \
    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
    | while read -r url ; do
        fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
        echo $acc
        echo "$url/$fname" >> $OUTPUT_FILE
      done;
done
