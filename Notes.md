## Test Data

**Transcript reference file**  
gencode.v28.transcripts.fa.gz ~65MB  
Downloaded from GENCODE release 28 (GRCh38.p12)  
```bash
wget "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"
```
* Note: refused connection error being encountered, currently unresolved  

**Genome reference file**  
GRCh38.primary_assembly.genome.fa.gz ~840MB  
Downloaded from GENCODE release 28 (GRCh38.p12)  
```bash
wget "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz"
```


**Paired read files**  
* fastq read files were taken from the lab directory /homes/wheelerlab2/Data/gEUVADIS_RNASeq  

## Kallisto
current version 0.44.0  
  
**references**  
https://pachterlab.github.io/kallisto/  
  
**Kallisto installation**  
requires conda installed and the bioconda channel be opened

```bash
conda config --add channels conda-forge
conda config --add channels bioconda

conda install kallisto
```
* Note: permission denied error encountered, resolved by running these commands with **sudo**  
* Note: conda is not currently in the PATH. Typing out the full path to conda resolves this, but is not optimal.

**running Kallisto**  
* Initial testing was run on the using the transcript fastq file taken from the GENCODE project, see **Test Data** for download instructions
* Paired reads ERR188030_1.fastq.gz (1151 MB) and ERR188030_2.fastq.gz (1136 MB) were used to test this data

1. Creating an index file   
Kallisto requires the proccessing of a transcriptome file. This is a one time process for a given reference file and should only take up to 10 minutes.
```bash
kallisto index -i name_of_index_file.idx gencode.v28.transcripts.fa.gz
```
* Runtime of this step on GENCODE data 3m38s
* Name of index file should be established here and used from this point on for .idx arguments
* Note: Kallisto and transcript file are not currently in the PATH. Typing out the full PATH to these items resolves this, but is not optimal.

2. Quantifying transcript abundances
 ```bash
 kallisto quant -i name_of_index_file.idx -o output_directory -b 100 ERR188030_1.fastq.gz ERR188030_2.fastq.gz
 ```
* Runtime with bootstrap = 100 (-b 100) 42m47s
* Runtime without bootstrap (-b 0 is the default) 2m37s
* Manual notes it is possible to increase the runtime by up to 15% but this has not been explored
* Note: read files may not be in a different directory than the wd. Typing out the full PATH to these items resolves this, but is not optimal.
* Note: Uncertain if kallisto has the ability to name the output files for this step may interfere with downstream processing - may be adequate to name a unique directory for each output but not optimal  
  
should be possible to assign genomic coordinates to the transcripts - possible solution to reads per gene question?


## STAR
version 2.6.0c   

**Running Star** 

1. Creating an index file  
first create a directory to place the index file
```bash
mkdir STAR_out
```
Next, STAR will not accept compressed genome files as input, unzip the reference genome with zcat or gunzip -c and append it to a new file.   
```bash
gunzip -c GRCh38.primary_assembly.genome.fa.gz > GRCh38.primary_assembly.genome.fa
```
Next create the index file   
```bash 
STAR --runMode genomeGenerate --genomeDir STAR_out --genomeFastaFiles Data/Gencode/GRCh38.primary_assembly.fa
```
* Note: This process takes a long time recommended that it's run with nohup or similar. To time with nohup run as follows
```bash
nohup bash -c "time STAR --runMode genomeGenerate --genomeDir STAR_out --genomeFastaFiles Data/Gencode/GRCh38.primary_assembly.fa"
```

* Note: STAR does perform WASP filtering - appends a tag to the end of alignments indicating whether it passed or failed.
  * Need to check if this requires WASP installation or if this feature is built in   
remains to be tested - couldn't download GENCODE test files due to connection refusal error
