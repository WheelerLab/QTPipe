
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

###running Kallisto  
* Initial testing was run on the using the transcript fastq file taken from the GENCODE project, see **Test Data** for download instructions
* Paired reads ERR188030_1.fastq.gz (1151 MB) and ERR188030_2.fastq.gz (1136 MB) were used to test this data

####1. Creating an index file  
Kallisto requires the proccessing of a transcriptome file. This is a one time process for a given reference file and should only take up to 10 minutes.
```bash
kallisto index -i name_of_index_file.idx gencode.v28.transcripts.fa.gz
```
* Runtime of this step on GENCODE data 3m38s
* Name of index file should be established here and used from this point on for .idx arguments
* Note: Kallisto and transcript file are not currently in the PATH. Typing out the full PATH to these items resolves this, but is not optimal.

**Creating an index file using a genome reference file**
Should be the same as creating an index file, just supply a genome file in place of the transcript file
* Important: Significant server delays were experienced during this, process do not attempt to run multiple jobs unless you want the server to be unuseable for a while
* May be dependent on the size of your reference file as this process uses significant memory

####2. Quantifying transcript abundances
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

## Useful commands

bash commands
```bash
wget "some_url_here" #retrieve data from a url

nohup your_commands_here --whatever --arguments #continues to run the process even if you leave the server

nohup taskset -c some_number your_command_here #runs your command on a specific core Note wheeler lab has 12 (0-11) cores

some_command & #runs command in the background

bash -c "Command_as_string" #interprets the string input as bash command with spaces delimiting positional arguments"

time Some_command #times how long it takes to execute Note does not work with nohup unless used in conjunction with bash -c

top #displays what processes are running typing 1 will display core use

`ps aux --sort=-%mem | head` #To see which processes are using the most memory:

#To free up memory, run these commands:
sync
su
echo 3 > /proc/sys/vm/drop_caches

//useful combination
nohup taskset -c # bash -c "time Some_Command --argument &" //runs your job in the background on a specified core, times it , and doesn't terminate when you leave the shell
```
