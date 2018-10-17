#origianlly written by Brianne Coffey Github@breecoffey
#modified by Ryan Schubert

import argparse
import gzip
import re
from os.path import expanduser
current = expanduser("~")

parser = argparse.ArgumentParser(description='Input & Output Files') #create the argument parser
parser.add_argument('--VCF', help='The VCF file to open') #variable for VCF file
parser.add_argument('--maf', help='Minor allele frequency threshold', const=1, type = float)
parser.add_argument('--r2', help='R2 threshold', const=1, type = float)
parser.add_argument('--outputdir', "-o", default=current, help='The VCF file to open') #variable for VCF file
args = parser.parse_args() #parse the arguments

#geno = "SNPGenotypes_" + args.VCF
#snpsloc = "SNPLoc_" + args.VCF

filtered = [] #array to store filtered snp info
loc_output = [] #array for snp chrom location output file

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
        snp_id = words[2] #grab snp id
	#loc_output.append(str(snp_id) + '\t' + "chr"+str(chrom_num) + '\t' + str(chrom_loc)) #append together for the location output file
        R2_index = -1
	af_index = -1
        values = words[7].split(';') #grab the 7th column with GT info for each sample
        af_index = 2 #set initial allele freq to index 2	
        if(str(values[2][0:4]) == 'MAF='): #check first 3 indicies for AF info, because some snps are missing first 2 values
            af_index = 2
        if(str(values[1][0:4]) == 'MAF='):
            af_index = 1
        if(str(values[0][0:4]) == 'MAF='):
            af_index = 0
	if(str(values[2][0:3]) == 'R2='): #check first 3 indicies for AF info, because some snps are missing first 2 values
            R2_index = 2
        if(str(values[1][0:3]) == 'R2='):
            R2_index = 1
        if(str(values[0][0:3]) == 'R2='):
            R2_index = 0
	R2_float = [float(i) for i in values[R2_index][3:].split(',')]
        maf_float = [float(i) for i in values[af_index][3:].split(',')]
        if (af_index != -1 and R2_index != -1) and (max(maf_float) >= args.maf and max(R2_float) >= args.R2): #check if AF is above 0.01 frequency
            unfiltered_geno = words[9:] #grab unfiltered genotype info
            genotypes = [] #create array for all genotype values to be stored
            #unfiltered_geno = unfiltered_geno[2:]

            for i in unfiltered_geno:
                j = i.split(':') #split info by colon
                zygosity = j[0].count("1")
		#append 2, 1, or 0 based on corresponding genotype
                if zygosity == 2:
                    genotypes.append('2')
                elif zygosity == 1:
                    genotypes.append('1')
                elif j[0] == '0/0' or j[0] == '0|0':
                    genotypes.append('0')
                elif j[0] == './.' or j[0] == '.|.':
                    genotypes.append('-1')
                else:
                    genotypes.append('-1')
                    continue
		
            filt = str(snp_id) #convert snpid to string
            if genotypes == []: #if no values, go to next snp
                continue
            for k in genotypes:
                filt += '\t' + str(k)
            filtered.append(filt) #append filtered GT values by tab delimited
            loc_output.append(str(snp_id) + '\t' + "chr"+str(chrom_num) + '\t' + str(chrom_loc))
	#write to output file

geno = args.outputdir + "SNPGenotypes.txt.gz"
snpsloc = args.outputdir + "SNPLoc.txt.gz"

with gzip.open(geno, 'wb') as output_file:
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
