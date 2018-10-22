#origianlly written by Brianne Coffey Github@breecoffey
#modified by Ryan Schubert

import argparse
import gzip
import re
from os.path import expanduser
current = expanduser("~")

parser = argparse.ArgumentParser(description='Input & Output Files') #create the argument parser
parser.add_argument('--VCF', help='The VCF file to open') #variable for VCF file
parser.add_argument('--maf', help='Minor allele frequency threshold', type = float)
parser.add_argument('--r2', help='R2 threshold', type = float)
parser.add_argument('--outputdir', "-o", default=current, help='The VCF file to open') #variable for VCF file
parser.add_argument('--tag', help='tag') #variable for VCF file
parser.add_argument('--sanger',action='store_true',help='files were sanger imputed. Default interprets as Michigan')
args = parser.parse_args() #parse the arguments

#geno = "SNPGenotypes_" + args.VCF
#snpsloc = "SNPLoc_" + args.VCF

filtered = [] #array to store filtered snp info
loc_output = [] #array for snp chrom location output file

if args.sanger == True:
    acc_metric="INFO="
    metric_length=5
    DS_index=2
else:
    acc_metric="R2="
    metric_length=3
    DS_index=1
with gzip.open(args.VCF, "rt") as file: #open the vcf file
    for line in file:
        if line.startswith('##'): #strip the meta-information
            continue
        if line.startswith('#'): #store the header
            header = line
            continue
        words = line.strip('\n').split('\t') #split row by tab 
        chrom_num = words[0] #grab chromosome number
        chrom_loc = words[1] #grab chromosome location
        if words[2] == '.':
            snp_id = chrom_num + ":" + chrom_loc #assign SNP ID if empty
        else:
            snp_id = words[2] #grab snp id
        acc_index = -1
        af_index = -1
        values = words[7].split(';') #grab the 7th column with GT info for each sample
#        af_index = 2 #set initial allele freq to index 2	
#        if(str(values[2][0:4]) == 'MAF='): #check first 3 indicies for AF info, because some snps are missing first 2 values
#            af_index = 2
        for i in values:
#            print(i)
            if( 'MAF=' in str(i)):
                af_index = values.index(i)
#                print("TRUE")
            if( acc_metric in str(i)):
                acc_index = values.index(i)
 	#if(str(values[1][0:4]) == 'MAF='):
        #    af_index = 1
        #if(str(values[0][0:4]) == 'MAF='):
        #    af_index = 0
        #if(len(values) == 3 and str(values[2][0:metric_length]) == acc_metric): #check first 3 indicies for AF info, because some snps are missing first 2 values
        #    acc_index = 2
        #if(str(values[1][0:metric_length]) == acc_metric):
        #    acc_index = 1
        #if(str(values[0][0:metric_length]) == acc_metric):
        #    acc_index = 0
        R2_float = [float(i) for i in values[acc_index][metric_length:].split(',')]
        maf_float = [float(i) for i in values[af_index][4:].split(',')]
        if (af_index != -1 and acc_index != -1) and (max(maf_float) >= args.maf and max(R2_float) >= args.r2): #check if AF is above 0.01 frequency
            unfiltered_dose = words[9:] #grab unfiltered genotype info
            dosages = [] #create array for all genotype values to be stored
            #unfiltered_geno = unfiltered_geno[2:]

            for i in unfiltered_dose:
                j = i.split(':') #split info by colon
                dose = j[DS_index]
                dosages.append(dose)
		
            filt = str(snp_id) #convert snpid to string
            if dosages == []: #if no values, go to next snp
                continue
            for k in dosages:
                filt += '\t' + str(k)
            filtered.append(filt) #append filtered DS values by tab delimited
            loc_output.append(str(snp_id) + '\t' + "chr" + str(chrom_num) + '\t' + str(chrom_loc))
	#write to output file

dose_out = args.outputdir + "/Chr" + args.tag + "MAF" + str(args.maf) + "R2" + str(args.r2) + "Dosages.txt.gz"
snpsloc = args.outputdir + "/Chr" + args.tag + "MAF" + str(args.maf) + "R2" + str(args.r2) + "SNPloc.txt.gz"

with gzip.open(dose_out, 'wb') as output_file:
    samples = header.split('\t') #sample IDs
    output_file.write(('id').encode('utf-8')) 
    samples = samples[9:] #remove first columns from header
    for name in samples:
        output_file.write(('\t'+str(name)).encode('utf-8'))

    for i in filtered:
        output_file.write((str(i) + '\n').encode('utf-8')) #write out filtered genotype info
with gzip.open(snpsloc, 'wb') as output_file: 
    output_file.write(("snp\tchr\tpos" + '\n').encode('utf-8')) #concat info for chrom position 
    for snp in loc_output:
        output_file.write((snp + '\n').encode('utf-8'))
