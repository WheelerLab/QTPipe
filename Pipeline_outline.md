# About    
this document is meant to give a detailed description of this pipeline. It will inform of three things:  
1. The intended workflow for RNA-seq data as well as the required input files
2. An outline of how each of the individual scripts of this pipeline works and descriptions of their full suite of options
3. Examples of how to run this pipeline, using both the individual parts as well as the composite superscript    


This outline is intended to be a living document, with examples and descriptions intended to be updated as the pipeline is updated. Where appropriate this outline will reference the manuals and studies that aided its development. Importantly there are many defaults and options set within this pipeline that are intended to make the lives of this lab in particular easier. Their purpose ranges from changing the number of samples for easy scalable testing to implementing defaults for quick command line runs. Where possible this document notes what is intended for in house usage and why. Individuals considering utilizing or referencing this pipeline may wish to make alterations accordingly.    

## Basic Flowgram

1. Subset your samples (filter)
2. Generate FastQC report
3. Align samples to both genomic and transcriptomic coordinates
4. eQTL & sQTL analysis  
      4a. Salmon quantification (parse & filter)  
      4b. Leafcutter spllicing analysis & file processing (parse & filter)  
      4c. SNP file processing (parse &filter)  
      4d. MatrixeQTL x2  

## Pipeline options  
In general The options of this pipeline fall into two categories, shared and specific. As they may sound, shared options are those that connect the different scripts of this pipeline. These options entail both raw data and its flow from one ource to another. This includes reference files, sample data, and the inputs and outputs from one script to another. On the pipeline level, these are all meant to be set with adequate defaults so that the only information the user typically must provide is the raw input data as well as the necessary reference files. These options include all required options that the user must supply the pipeline. On the level of the individual scripts, these are options that do not have a set default (although they are set in the pipeline superscript).   The majority of options fall into this category. 

In contrast are specific options which correspond to features found only in a particular script. These options include such features as deciding the minimum clustering size for leafcutter analysis or specifying what library type Salmon should use for its analysis. These arguments tend to be optional, with most if not all having defined within their own scripts a default value. They can still easily be accessed on the super script level.

# Individual Scripts & Their Options

## list_samples

The list_samples script, so aptly named, lists the subset of samples you wish to process in the ofrmat required by the rest of the pipeline. While creating this pipeline it was found that the sample RNA data was tagged differently from how they were labelled within the SNP files. This posed a problem for downstream processing as matrixEQTL requires that sample names be consistent between our SNP and our gene files. As a solution a simple "translation" step was added within the star_loop script. This would simply rename a given sample to its corresponding genotyped sample. this requires a two column file, with the translated names in the first column and the corresponding untranslated names in the second cloumn. This script makes it simple to generate this list that will be used by star as well as all other codes downstream of it. Currently the file that it uses to translate samples is hardcoded. expect this as an options update in the future. Additionally this script also filters out those samples that are not already genotyped and thus cannot be used for downstream analysis. This filtering step is not currently optional and anyone wishing to process said samples will need to do so independently of genotyped samples.

The renaming step creates a problem in portabiliy. The translation step had to be incorporated in all code downstream of this application. Anyone wishing to run this pipeline on samples that do not need this name change (or whose name change follows a different format than the one specified) will run into errors while running this pipeline. A more portable code is in the works, but for now do not hesitate to contact the developer who will gladly help you with this issue. 

**options**
-i or --inputdirectory
> Specifies the directory containing the fastq files to be processed. This is a **_required option._** Currently this option defaults to the developers directory for quick and easy in lab usage.

-o or --outputdirectory) 
>specifies the outtput file you'd like to send file. Can include the full path of the file. Otherwise defualts to generating a file at ~/sample_list.txt. **_optional argument_**

-s or --samplenumber
>Specifies how many samples you wish to process. Must be less than or equal to the number of samples/sample pairs. By default this option selects 89 samples as that is the total number of test samples available. **_recommended argument._** An option to simply run the pipeline on the entire directory of fastq samples is not currently in placed, but this can be achieved by inputing the total number of samples.

## fastqc_loop

This script iterates FastQC over the entirety of the sample list including all paired end samples that share a sample tag. It then outputs two summary files, one containing the pass/fail stats for all samples and another containing the pass/fail stats for only those fastq files with at least 1 failing statistic.

-f or --fastqdir
>directory containing the fastq files to be processed. This is a **_required option_**

-o or --output
>directory where you'd like to send all your quantification files

-s or --samplelist
>a list of all the samples you're going to process. This sample list should follow the format as laid out in the list_samples code


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

## qqnormparser.py 

## Salmon_loop  

-a or --annotation 
>Defines the full path of the annotation file for the genome being used. This is a **_required argument_**.

-i | --inputdirectory)
>Defines the full path to the Directory containing all the star aligned.to.transcriptome files to be used. Assumes that the alignment files are organized as /$InputDirectory/${SAMPLE}\_star/Aligned.toTranscriptome.bam etc. for every individual sample. For each run of the pipeline this input directory should always be the same as the the output directory for STAR. It is currently planned (but not yet implemented) for this option to default as the STAR output. As such this option is **_required if running independent of the pipeline_** otherwise will always use the star output of the given run. This option is coded in to allow for increased independence of this script. This option allows for ease of comparison of different alignment types (say those made with different aligners or different reference files)

-l | --librarytype)
>The lib type defines a feature granted by Salmon. Salmon requires that you input what type of alignment you have, whether paired end, single end etc. Salmon provides an A option by which it can automatically infer what type of reads are used. This is used as a default here. As such this argument is **_optional_** but recommended if the lib type is known. For full list of library types please see the [salmon documentation](http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype).

-o | --outputdirectory) 
>Output directory where you would like to send your quantification files. Outputs will be structured as $OutputDirectory/${Sample}\_salmon etc. with each ${Sample}\_salmon being a directory that contains the quant.sf files for each sample. These files will be independently compiled into a single expression file in order for matrixEQTL processing. **warning** it is very important that this directory not have any unwanted quantification files already within it, particularly when running within the pipeline as this will create problems downstream. This is thus a **_recommended_** option

-s | --samplelist) 
>The list of samples you would like to process with salmon. This option is provided in order to let one easily subset the alignments you wish to quantify. To run on all samples, simply feed it a complete list of samples. Sample lists assume the structure of the fastq list generated by the star_loop script. This is a **_required_** argument.

-t | --transcript)
>Defines the full path to the transcript file of the genome. This is a **_required_** option.

## Salmon_Parser.R

When Salmon generates its quantification outputs, it does so by generating a new directory for each sample. In order to process these files with MatrixEQTL it is necessary to combine each output into one single file. Salmon parser runs on an entire directory containing the subdirectories output by salmon. It is important to note that salmon_parser will parse _all_ quantification files found within this directory and not just the ones within the original sample list. It is thus important to keep all your quantification subsets in different directories if you don't want them included in your quantification. At the same time it is important to filter out any genes with significantly low variance to avoid large volume of non significant genes-SNP pairs from being identified by MatrixEQTL. The parser requires the input of an annotation (.gtf) file as well as the quantification directory containing all of the sample quantifications to be processed. Will output an expression file and a location file for each chromosome

-a or --annotation file
>Defines the full path of the annotation file for the genome being used. This is a **_required argument_**. Salmon_parser uses the annotation to assign gene locations and chr num to the quantified genes as salmon fails to write this out.

-m or --meanthreshold
>The minimum mean expression that the gene must have across all samples. By default this threshold is 0.1. Any gene with mean expression below this threshold is removed from both the expression and the location file. This argument is **_optional_** however it is recommended that threholds above 0 be used if provided.

-o or --outputdir
>The path to the output directory. This is **_required_** as no default exists.

-q or --quantdir
>The path to the input directory containing the salmon quantifications you would like to process. Note that Salmon_parser will run on *all* quantification files within this directory and its subdirectories and not just those most recently generated by salmon. This is a **_required_** option.

-s or --scaledthreshold
>The threshold for normalized variance across gene samples for filtering. As variance is both positive and negative, this filters out those genes whose normalized variation lies between +threshold to -threshold inclusively. Variation is normalized using R's scale(0 function. This option has a defualt of 0. In testing it has been found that this option is very sensitive in filtering out genes and that it is difficult to judge an appropriate threshold for filtering. As such it is actually **_not recommended_** for use and will likely be removed in later iterations.

-v or --variance threshold
>The threshold for variance across gene samples for filtering. As variance is both positive and negative, this filters out those genes whose variation lies between +threshold to -threshold inclusively. This option is set with a default of 0.1. This default was found to be adequate to remove overrepresented insignificant gene-snp pairs from matrix-eqtl analysis that would skew our results. This argument is **_optional_**

## SNP Parser
