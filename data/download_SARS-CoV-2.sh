#!/bin/bash

batch_size=10  # number of lines to process at a time
line_count=0  # counter for current line number
batch=()  # array to store current batch of lines
batch_num=0

while read -r line; do
    if [[ $line_count -eq $batch_size ]]; then
        # process batch
        #echo "Processing batch: ${batch[*]}"
        
        # loop over elements of batch
        for element in "${batch[@]}"; do
            # do something with element
            #echo "Processing element: $element"
	        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${element}&rettype=fasta" >> ${batch_num}.fa
        done
        
        # clear batch
        batch=()
        line_count=0
        ((batch_num++))
    fi
    # add line to batch
    batch+=("$line")
    ((line_count++))

    # wait for 10 seconds
    sleep 10
done < accession_numbers.SARS-CoV-2.txt

# process last batch (if any)
if [[ $line_count -gt 0 ]]; then
    #echo "Processing batch: ${batch[*]}"
    
    # loop over elements of last batch
    for element in "${batch[@]}"; do
        # do something with element
        #echo "Processing element
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${element}&rettype=fasta" >> ${batch_num}.fa
    done
fi

#sed '/^$/d' x.fa > sars_cov_2_sequences.fa
#rm x.fa
