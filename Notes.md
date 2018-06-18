
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

**1. Creating an index file**  
Index generation is  a single step that should be a one time process for a given reference file. Should always be performed prior to the rest of the transcriptome quantification. This is a one time process for a given reference file and should only take up to 10 minutes.
```bash
kallisto index -i name_of_index_file.idx gencode.v28.transcripts.fa.gz
```
* Runtime of this step on GENCODE data 3m38s
* Name of index file should be established here and used from this point on for .idx arguments
* Note: Kallisto and transcript file are not currently in the PATH. Typing out the full PATH to these items resolves this, but is not optimal.

**2a. Quantifying transcript abundances**
 ```bash
 kallisto quant -i name_of_index_file.idx -o output_directory -b 100 ERR188030_1.fastq.gz ERR188030_2.fastq.gz
 ```
* Runtime with bootstrap = 100 (-b 100) 42m47s
* Runtime without bootstrap (-b 0 is the default) 2m37s
* Manual notes it is possible to increase the runtime by up to 15% but this has not been explored
* Note: read files may be in a different directory than the wd. Typing out the full PATH to these items resolves this, but is not optimal.
* Note: Uncertain if kallisto has the ability to name the output files for this step may interfere with downstream processing - may be adequate to name a unique directory for each output but not optimal  
  
**2b. Looping kallisto quantification**

```bash
./loop_kallisto /path/to/fastq/directory/
```
* Index generation should always be performed prior to this
* Code should be updated to ask if this index has been generated yet
* Code calls the reads_per_gene_kallisto.py after each loop
* Currently this code is oriented to only perform on the first 10 sample pairs in a library but can be easily updated
* currently does generate bamfiles - update as user option?
* Uses a default bootstrap of 100 - can be easily changed within code - update as user option?

Speeds at different levels of bootstrap using 10 sample pairs
* @bootstrap = 0, run time = 47m30s
* @bootstrap = 10, run time = 99m11s
* @bootstrap = 100, run time = 531m41s

**3. Quantifying gene abundance**
```bash
python3 reads_per_gene_kallisto.py /path/to/abundance.tsv /desired/output/path/out_name.tsv
```
* Does not currently have default arguments - can either be added or simply fed withing loop_kallisto
* Number of bootstraps does not affect TPM only accuracy estimates

## STAR & RSEM
version 2.6.0c   

**Running Star** 

**1. Creating an index file**  
first create a directory to place the index file
```bash
mkdir StarTest
```
Next, STAR will not accept compressed genome files as input, unzip the reference genome with zcat or gunzip -c and append it to a new file.   
```bash
gunzip -c GRCh38.primary_assembly.genome.fa.gz > GRCh38.primary_assembly.genome.fa
```
Next create the index file   
```bash 
STAR --runMode genomeGenerate --genomeDir StarTest --genomeFastaFiles Data/Gencode/GRCh38.primary_assembly.fa
```
* Note: This process takes a long time recommended that it's run with nohup or similar. To time with nohup run as follows
```bash
nohup bash -c "time STAR --runMode genomeGenerate --genomeDir StarTest --genomeFastaFiles Data/Gencode/GRCh38.primary_assembly.fa"
```
* Approximate time 112m

**2.Map the gzipped FASTQ files outputting unsorted and coordinate-sorted BAMs**  
\* This process was taken from Alternate Protocol 7 of Dobin & Gingeras (2016)  
[Mapping RNA-Seq Reads with STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/)
* Requires a gtf file
  * Should be possible to directly pass STAR the transcript file in place of the genome file and gtf file but this has not been tested

```bash
STAR --genomeDir ~/StarTest/ --sjdbGRFfile ~/DATA/GenCode/gencode.v28.annotations.gtf --readFilesIn /home/wheelerlab2/Data/gEUVADIS_RNASeq/ERR188030_1.fa.gz /home/wheelerlab2/Data/gEUVADIS_RNASeq/ERR188030_2.fa.gz --readFilesCommand zcat --quantMode TranscriptomeSAM
```
 * Time for this process on one pair of fastq files: 26m7.41s  
 
 **3. Prepare the RSEM reference files**
 
 ```bash
 ~/RSEM/rsem-prepare-reference --gtf ~/DATA/GenCode/gencode.v28.annotations.gtf ~/DATA/GenCode/GRCh38.primary_assembly.genome.fa ~/RSEM/ref
 ```
 
* Note: STAR does perform WASP filtering - appends a tag to the end of alignments indicating whether it passed or failed.
  * Need to check if this requires WASP installation or if this feature is built in   

**4. Run RSEM quantification on the STAR transcriptomic BAM file**
```bash
~/RSEM/rsem-calculate-expression --bam --no-bam-output -p 12 --paired-end --forward-prob 0 ~/transcriptStar_out/Aligned.toTranscriptome.out.bam ~/RSEM/ref ~/transcriptStar_out/Quant >& ~/transcriptStar_out/rsem.log
```
* Approximate time 106m53.44s

## Useful commands

bash commands
```bash
/usr/local/bin/anaconda3/bin/python3 #path to python3 on wheeler lab

chmod a+x script_name #mak ea custom script executable

wget "some_url_here" #retrieve data from a url

nohup your_commands_here --whatever --arguments #continues to run the process even if you leave the server

nohup taskset -c some_number your_command_here #runs your command on a specific core Note wheeler lab has 12 (0-11) cores

some_command & #runs command in the background

bash -c "Command_as_string" #interprets the string input as bash command with spaces delimiting positional arguments"

time Some_command #times how long it takes to execute Note does not work with nohup unless used in conjunction with bash -c

top #displays what processes are running typing 1 will display core use

kill -9 PID #kill a process using its pid number
kill -9 $(pgrep -f command_name) #useful for commands that execute loops - exits the whole loop not just one iteration

`ps aux --sort=-%mem | head` #To see which processes are using the most memory:

#To free up memory, run these commands:
sync
su
echo 3 > /proc/sys/vm/drop_caches

//useful combination
nohup taskset -c # bash -c "time Some_Command --argument &" //runs your job in the background on a specified core, times it , and doesn't terminate when you leave the shell
```
## Bash Scripting Tips & Tricks     

**setting default arguments for options**

```bash
LINES=50       # Default number of lines saved.

if [ -n "$1" ]
# Test whether command-line argument is present (non-empty).
then
  lines=$1
else  
  lines=$LINES # Default, if not specified on command-line.
fi  
```

**Double quote your variables**
```bash
"${var}"
```
Prevents misinterpretation of special characters and values containing spaces

**$ argument interpolation**     

```bash
echo "Last program's return value: $?"
echo "Script's PID: $$"
echo "Number of arguments passed to script: $#"
echo "All arguments passed to script: $@"
echo "Script's arguments separated into different variables: $1 $2..."
```

**check the number of arguments is correct**    
```bash
E_WRONG_ARGS=85
script_parameters="-a -h -m -z"
#                  -a = all, -h = help, etc.

if [ $# -ne $Number_of_expected_args ]
then
  echo "Usage: `basename $0` $script_parameters"
  # `basename $0` is the script's filename.
  exit $E_WRONG_ARGS
fi
```
**reading user input**
```bash
read variable_name #Read one line from the standard input, (or from a file) and assign the word(s) to variable name(s).
```
