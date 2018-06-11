import argparse as ap
parser = ap.ArgumentParser()

parser.add_argument('transcript_Abundance', help = 'file path of abundance.tsv file to be parsed')
parser.add_argument('Gene_abundance_out', help = 'desired file path for gene abundance output file') 
args = parser.parse_args()

genes = {} #empty dictionary to append genes to
with open(args.transcript_Abundance, 'r') as abundance:
	next(abundance)
	for line in abundance:
		fields = line.strip('\n').split('\t') #split the file by tabs for easy access to each column value. The last two fields contain values of interest.
		IDs = fields[0].split('|') #Now we check the list of IDs from the first column. The first ID in the list is the transcript ID. The second item in the list is the gene ID it maps to.
		if IDs[1] in genes: #check if the gene ID is already in the dictionary
			#if yes we are going to update three values
			genes[IDs[1]][0] += '|' + IDs[0] #update the list of transcripts that map to this gene
			genes[IDs[1]][1] += float(fields[-2]) #update the est_counts of the gene
			genes[IDs[1]][2] += float(fields[-1]) #update the tpm of the gene
		else:
			genes[IDs[1]] = [IDs[0],float(fields[-2]),float(fields[-1])] #otheriwse create a new dict entry with those three values as a list

with open(args.Gene_abundance_out, 'w') as output:
	output.write('Gene_ID\tTransctipt_ID(s)\test_counts\ttpm\n') #write out header
	for key in genes.keys():#now write out the four fields
		output.write(key + '\t')
		output.write(genes[key][0] + '\t')
		output.write(str(genes[key][1]) + '\t')
		output.write(str(genes[key][2]) + '\n')	

