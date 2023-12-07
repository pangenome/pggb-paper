#!/bin/bash

# Check if a FASTA file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <prefix>"
    exit 1
fi

FASTA="$1"
ACC="$2"

awk -v ACC=$ACC '
/^>/ {
    # Extract the contig name (first field), chromosome number, and description
    split($1, arr, ".");
    contig = arr[1];
    contig = substr(contig, 2)

    if (match($0, /chromosome ([XY0-9]+)/, chr)) {
        chromosome = "chr" chr[1];
    } else if (match($0, /chromosome: ([XY0-9]+)/, chr)) {
        chromosome = "chr" chr[1];
    } else {
        chromosome = "unknown";
    }

    description = (index($0, "unlocalized") > 0) ? "unlocalized" : "primary";

    # Print the new formatted header
    print ">" ACC "#1#" contig "-" chromosome "-" description;
}
!/^>/ {
    # Print sequence lines as they are
    print $0;
}
' "$FASTA"
