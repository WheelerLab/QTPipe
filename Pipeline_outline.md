# About    
this document is meant to give a detailed description of this pipeline. It will inform of three things:  
1. The intended workflow for RNA-seq data as well as the required input files
2. An outline of how each of the individual scripts of this pipeline works and descriptions of their full suite of options
3. Examples of how to run this pipeline, using both the individual parts as well as the composite superscript    


This outline is intended to be a living document, with examples and descriptions intended to be updated as the pipeline is updated. Where appropriate this outline will reference the manuals and studies that aided its development. Importantly there are many defaults and options set within this pipeline that are intended to make the lives of this lab in particular easier. Their purpose ranges from changing the number of samples for easy scalable testing to implementing defaults for quick command line runs. Where possible this document notes what is intended for in house usage and why. Individuals considering utilizing or referencing this pipeline may wish to make alterations accordingly.    

## Basic Flowgram

## Star_loop  

The star_loop script is intended to loop STAR alignment steps to both genomic and transcriptomic contexts. The basic required files are an uncompressed annotation file (.gtf), an uncompressed genome file (.fa), and any number of paired end compressed read files (.fa.gz). Because this script is an implementation of STAR, many of the options included are based on the standard STAR options with reasonable defaults provided for this lab.

The star loop now has a built in method to rename samples so that they match with SNP sample names. Currently this is hardcoded in for in lab usage. More flexible applications need to be explored still. Similarly the script now filters the input samples to ensure they are within the genotyped pool.

**Options**  
-a or --annotation
>This option specifies the full path to the annotation file. This is a **_required_** option with no defaults currently. In future implementations it is reasonable to expect that the annotation file may not differ between between runs of this pipeline. As such it is reasonable to expect that a default will be provided for in lab usage.    
      
-g or --genomefastafiles
>This option specifies the full file path to the genome files to be used in indexing. This option is **_sometimes required_, only when the --runindexing option is also specified.** Meaning if indexing is being run then this file must be provided. Similar to the annotation file it is reasonable to expect that this file may not differ between between runs of this pipeline. As such it is reasonable to expect that a default will be provided for in lab usage. Right now this option only accepts one file for genomic fasta files. Future implementations may accept multiple if this proves necessary in lab.

-i or --inputdirectory 
>This option specifies the directory containing the fastq files to be processed. This is a **_required_** option with no defaults currently. Rather than individual fasta files, it is assumed that this pipeline will be run on entire directories containing an arbitrarily large number of samples. Similar to the annotation file it is reasonable to expect that this file may not differ between between runs of this pipeline. As such it is reasonable to expect that a default will be provided for in lab usage.    

-id or --indexingdirectory 
>This option has two uses. When supplied without the --runindexing option, this option should specify the directory containing the correct indexing files to be used by STAR and indexing will not be run. The second purpose of option is to specify the desired indexing output directory when the --runindexing option is provided. In both cases this option is **_optional when only 1 genome has been/is being indexed._** This option is provided with a default of `~/starIndex/` since index files are expected to be reused between runs.  If only 1 genome is being studied then this argument is not necessary. Additionally a warning is provided for redundant names when inexing is being run. This is to help mitigate unecessary indexing runs, which are time and space inefficient.

-o or --outputdirectory)
>This option specifies the output path of the STAR Alignments. Becuase STAR alignments are output into individual subdirectories for each read file, this option is mostly organizational. This argument is **_optional_** as a default output of ~/starAlignments/ is provided.

-ri or --runindexing
>This option declares whether the user wishes to run indexing or not. This flag does not have any positional arguments .As this is a time intensive step, by default indexing will not be run.  As such this argument is explicitly **_required only when indexing has not been run._** Additionally, when supplied it is **_required that the --genomefastafiles flag is supplied as well._**

-s or --samplenumber)
>This option specifies the number of samples within your fastq directory you wish to run. This option is **_optional_** with a default value of 10. This option is largely provided for in lab testing of scripts. It is meant to make testing scalable for comparison of different sample sizes. The default is set at 10 since 1)This sample size allows for reasonable estimation of larger sample sizes and 2) It is low enough that if one forgets to set the sample size it will execute in a reasonable amount of time. Future implementations should expect an argument to run on the entire directory.

**features not yet implemented**
argument to run on an entire directory  
option for different types of star output (SAM vs BAM vs sorted vs unsorted)
defaults for several options
option to not run alignment (if already run) in super script
