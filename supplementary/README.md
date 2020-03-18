## How to run repertoire assembly tools

### MIXCR

We run MIXCR in RNA-Seq mode. Details are here http://mixcr.readthedocs.io/en/latest/rnaseq.html

We use the following commands to runMIXCR:

```
mixcr align -p rna-seq -r log_align.txt $1 $2 alignments.vdjca
mixcr assemblePartial-p alignments.vdjca alignments_rescued_1.vdjca 
mixcr assemblePartial alignments_rescued_1.vdjca alignments_rescued_2.vdjca
mixcr assemble -r log_assemble.txt alignments_rescued_2.vdjca mixcr.clns
mixcr exportClones mixcr.clns mixcr.txt
```

where $1,$2 are fastq files with paired-end RNA-Seq reads. Please note, RNA-Seq reads were simulated as a mixture of transcriptomic reads and reads derived from BCR and TCR transcripts.

### IMSEQ

We run IMSEQ with default parameters for IGH, IGK, IGK, TCRA, and TCRB chains using the following commands:

```
imseq -ref Homo.Sapiens.TRA.fa -o output-file_TCRA.tsv 
imseq -ref Homo.Sapiens.TRB.fa -o output-file_TCRB.tsv 
imseq -ref Homo.Sapiens.IGH.fa -o output-file_IGH.tsv 
imseq -ref Homo.Sapiens.IGK.fa -o output-file_IGK.tsv 
imseq -ref Homo.Sapiens.IGL.fa -o output-file_IGL.tsv 

```

IMSEQ cannot be applied directly to RNA-Seq reads because it was originally designed for targeted sequencing of B or T cell receptor loci. Thus, to independently assess and compare accuracy with ImReP, we only ran IMSEQ with the receptor-derived reads.  When we run IMSEQ on simulated data, we provide IMSEQ with simulated reads derived from BCR or TCR transcripts as a fastq file. When we run IMSEQ on real RNA-Seq reads,  we provide fastq file with unmapped reads plus reads mapped to BR and TCR loci.

## Extract CDR3 assembled by each of the tools

We have extracted full-length CDR3 sequences based on the definition of CDR3. CDR3  is defined as the sequence of amino acids between the cysteine ( C ) on the right of the junction and phenylalanine ( F ) (for all TCR chains and immunoglobulin light chains) or tryptophan ( W ) (for IGH) on the left of the junction.


### ImReP


```
awk '{if ($2=="IGH") print}' SRR5248347.cdr3 | awk '{print $1}' | sort | uniq >SRR5248347_IGH.cdr3
awk '{if ($2=="IGK") print}' SRR5248347.cdr3 | awk '{print $1}' | sort | uniq >SRR5248347_IGK.cdr3
awk '{if ($2=="IGL") print}' SRR5248348.cdr3 | awk '{print $1}' | sort | uniq >SRR5248348_IGL.cdr3
```


### MIXCR

To extract CDR3s from each of the chains we use the following commands:

```
grep TRA mixcr.txt | grep -v "TRB" | grep -v "TRG" | grep -v "TRD" | cut -f 33 | sort | uniq  | grep -v "*" >mixcr_TCRA.txt
grep IGH mixcr.txt | grep -v "IGL" | grep -v "IGK" | cut -f 33 | grep -v "*" |  grep "^C" | grep "W$" | sort | uniq >mixcr_IGH.txt 
grep IGK mixcr.txt | grep -v "IGL" | grep -v "IGH" | cut -f 33 | grep -v "*" |  grep "^C" | grep "F$" | sort | uniq >mixcr_IGK.txt 
grep IGL mixcr.txt | grep -v "IGK" | grep -v "IGH" | cut -f 33 | grep -v "*" |  grep "^C" | grep "F$" | sort | uniq >mixcr_IGL.txt 
```


### IMSEQ

To extract CDR3s from each of the chains we use the following commands:

```
grep IGH output-file.tsv | awk -F "\t" '{print $11}' | grep -v "cdr | grep -v "*" | sort | uniq >imseq_igh.txt
grep TRA output-file.tsv | awk -F "\t" '{print $11}' | grep -v "cdr | grep -v "*" | sort | uniq >imseq_tcra.txt
```


## Simulate RNA-Seq data as mixture of transcriptomic and receptor-derived reads

* Script to generate IGH and TCRA transcripts is available : [simulateTranscriptsClonotypes.py](https://github.com/smangul1/ImRep-Gtex/blob/master/simulateTranscriptsClonotypes.py)

* Script to simulate reads from IGH and TCRA transcripts is available : [simulateReads.sh](https://github.com/smangul1/ImRep-Gtex/blob/master/simulateReads.sh)

More details on how simulated data is generated are available in the manuscript

## Compare CDR3s from RNA-Seq and TCRB-Seq data

We have downloaded TCRB-Seq data from here https://clients.adaptivebiotech.com/pub/Liu-2016-NatGenetics. Data was prepared by Li, Bo, et al. "Landscape of tumor-infiltrating T cell repertoire of human cancers." Nature genetics 48.7 (2016): 725-732.

It contains 3 TCRB-Seq samples from 3 individuals from TCGA study. From TCRB-Seq data we have selected full-length CDR3s and excluded partial CDR3s using this command:

```
awk '{print $2}' TCGA-CZ-5463.tsv | grep -v Out | sort | uniq | grep -v "amino" | grep "^C" | grep "F$"
```

* For sample TCGA-CZ-4862 we have excluded 2083 partial CDR3s (e.g. CARSLI) and saved 22391 full-length CDR3s
* For sample TCGA-CZ-4862 we have excluded 3 partial CDR3s (e.g. CARSLI) and saved 748 full-length CDR3s
* For sample TCGA-CZ-4862 we have excluded 35 partial CDR3s (e.g. CARSLI) and saved 5959 full-length CDR3s

Full-length CDR3 assembled by repertoire assembly tools from RNA-Seq reads were extracted as described in Section "Extract CDR3 assembled by each of the tools"

One should note, the number of complete CDR3s fully matching CDR3s obtained by TCRB-Seq in our study  are not fully comparable with the results reported in Li, Bo, et al. , where CDR3 sequences are considered to match CDR3s from TCRB-Seq if at least 6 amino acids are matched.

## Compare CDR3s from RNA-Seq and BCR-Seq data

We perfrom a direct comparison between RNA-Seq and targeted BCR-Seq, by using RNA-Seq and BCR-Seq samples from the same individuals (n=18). Data was prepared by : Lombardo, Katharine A., et al. "High-throughput sequencing of the B-cell receptor in African Burkitt lymphoma reveals clues to pathogenesis." Blood Advances 1.9 (2017): 535-544.

CDR3 clonotypes from BCR-Seq were assembled by immunoSEQ Analyzer (https://www.adaptivebiotech.com/). CDR3 clonotypes from RNA-Seq were assembled by Imrep and MixCR. 

We have downloaded BCR-Seq CDR3s from here (https://clients.adaptivebiotech.com/pub/lombardo-2017-bloodadvances, Analysis 3).
RNA-Seq data was downoaded from SRA archive, https://www.ncbi.nlm.nih.gov/bioproject/PRJNA374464. 


The metadata, including SRA and bio names was downloaded here
- https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP099346

SRA files were converged to fastq file (_1.fastq.gz and _2.fastq.gz) using the following command:

```
while read line ; do echo "$PWD/../../sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump --gzip --split-files -O ${line} ${line}.sra">run_${line}.sh;done<samples.txt 
```


We extract IGH-CDR3s from immunoSEQ Analyzer files using the folowwing command:

```
grep IGH_ 009-0103.raw >009-0103_IGH.txt
And across all samples:
while read line; do grep IGH_ ${line}.raw >${line}_IGH.txt;done<samples_18.txt
while read line;do awk '{if ($2=="IGH") print}' ${line}.cdr3>${line}_IGH.cdr3;done<samples.txt
```



We calculated immune diversity of the individuals based on BCR-Seq and RNA-Seq data


### Imrep

To rename the files from SRA to real names
```
while read line; do new=$(echo $line| awk '{print $2}');old=$(echo $line| awk '{print $1}'); mv ${old}.cdr3 ${new}.cdr3;done<../sample_SRA_bio_new.txt 
ls *cdr3 | awk -F ".cdr3" '{print $1}'>samples.txt
for f in *_IGH.cdr3;do n=$(awk '{s+=$3} END {print s}' $f); echo $f,$n;done

```
