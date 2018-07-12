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
>This option has two uses. When supplied without the --runindexing option, this option should specify the directory containing the correct indexing files to be used by STAR and indexing will not be run. The second purpose of option is to specify the desired indexing output directory when the --runindexing option is provided. In both cases this option is **_optional when only 1 genome has been/is being indexed._** This option is provided with a default of `~/starIndex/` since index files are expected to be reused between runs.  If only 1 genome is being studied then this argument is not necessary. Additionally a warning is provided for redundant names when inexing is being run. This is to help mitigate unecessary indexing runs, which are time and space inefficient. However it is specifically because of this prompt that usage of this option is  **_strongly recommended when running in the background_** as if the command prompt comes up while in th

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
overwrite/do not replace option  

## leaf_loop  

-i or --inputdirectory
>Directory containing the star alignment files. Assumes that the alignment files are organized as /$InputDirectory/${SAMPLE}\_star/Alignment.bam etc. for every individual sample. For each run of the pipeline this input directory should always be the same as the the output directory for STAR. It is currently planned (but not yet implemented) for this option to default as the STAR output. As such this option is **_required if running independent of the pipeline_** otherwise will always use the star output of the given run. This option is coded in to allow for increased independence of this script. This option allows for ease of comparison of different alignment types (say those made with different aligners or different reference files)

-l or --lengthmaximum  
>Defines the maximum length of the introns to be used in clustering. This argument is **_optional_** with a default value set at 100kb. 100kb is the default defined within the leafcutter script.

-m or --minimumclustering 
>Defines the minimum number of reads that must map to a cluster to support it. This argument is **_optional_** with a default set at 30. 30 is the default that is defined within the leafcutter script.

-p or --prefix   
>Defines the prefix or tag you wish to use for clustering script. The prefix to be used to define your perindividual counts file. This argument is **_optional but recommended_** when running multiple iterations of leafcutter. The default output for this option is simply defined as leafcutter ie outputfile will look like leafcutter_perind.counts.gz

-pc or --principalcomponents
>Defines the number of Principle components that leafcutter generates to use as covariates when running QTL mapping.  **_optional_** This argument is included because it was provided for by the original leafcutter script. The utility of this option is questionable and needs to be further researched for use with matrixEQTL.

-r or --ratiominimum
>Defines the minumum proportion of reads that must map to a cluster to support it. This argument is **_optional_** with a default set at 0.001. 0.001 is the default that is defined within leafcutter.

-s or --samplelist
>Defines the list of alignment samples you wish to use. this option is **_required if running independent of the pipeline_** otherwise will always use all those samples aligned by star in a given run. This option is intended to increase the flexibility of this script and allow for ease of comparison. This option makes it easy to run leafcutter on different subsets of alignments. Sample lists assume the structure of the fastq list generated by the star_loop script.


-w or --writeto 
>Defines the output directory for leafcutter to write to. This argument is **_optional but recommended._** Current Default is set at ~/leaf_out/ This option is largely organizational to allow for ease of use.

**bugs and tbi**  
argument to run on an entire directory  
means to easily generate sample list  
overwrite/do not replace option  



# Salmon_loop  

-a or --annotation 
>Defines the full path of the annotation file for the genome being used. This is a **_required argument_**.

-i | --inputdirectory)
>Defines the full path to the Directory containing all the star aligned.to.transcriptome files to be used. Assumes that the alignment files are organized as /$InputDirectory/${SAMPLE}\_star/Aligned.toTranscriptome.bam etc. for every individual sample. For each run of the pipeline this input directory should always be the same as the the output directory for STAR. It is currently planned (but not yet implemented) for this option to default as the STAR output. As such this option is **_required if running independent of the pipeline_** otherwise will always use the star output of the given run. This option is coded in to allow for increased independence of this script. This option allows for ease of comparison of different alignment types (say those made with different aligners or different reference files)

-l | --librarytype)
>The lib type defines a feature granted by Salmon. Salmon requires that you input what type of alignment you have, whether paired end, single end etc. Salmon provides an A option by which it can automatically infer what type of reads are used. This is used as a default here. As such this argument is **_optional_** but recommended if the lib type is known. For full list of library types please see the [salmon documentation](http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype).

-o | --outputdirectory) #directory where you'd like to send all your quantification files
>Output directory where you would like to send your quantification files. Outputs will be structured as $OutputDirectory/${Sample}\_salmon etc. with each ${Sample}\_salmon being a directory that contains the quant.sf files for each sample. These files will be independently compiled into a single expression file in order for matrixEQTL processing. **_optional_**

-s | --samplelist) 
>The list of samples you would like to process with salmon. This option is provided in order to let one easily subset the alignments you wish to quantify. To run on all samples, simply feed it a complete list of samples. Sample lists assume the structure of the fastq list generated by the star_loop script. This is a **_required_** argument.

-t | --transcript)
>Defines the full path to the transcript file of the genome. This is a **_required_** option.

