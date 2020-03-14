#!/bin/bash

./simLibrary -x ${coverage} bcr-seq.transcripts..fasta | ./simNGS -o fastq -p single -n ${read_length} ${runfile} > simulated_reads_${read_length}.fq
