from itertools import islice
import gzip as gz
import argparse as ap
parser = ap.ArgumentParser()

parser.add_argument('target_counts', help = 'file path of targe counts file generated in between leafcutter scripts')
parser.add_argument('associated_names', help = 'the file E-GEUV-1.sdrf.txt containing all the name associations') 
parser.add_argument('outfile', help = 'file path of the renamed headers file')
args = parser.parse_args()

names_dict = {}
with open(args.associated_names, 'r', encoding = 'utf-8') as names:
	next(names)
	for line in names:
		fields = line.strip('\n').split('\t')
		names_dict[fields[27] + ".staraligned"] = fields[0]


with gz.open(args.target_counts, 'rt', encoding='utf-8') as target:
	header = target.readline().strip('\n').split()
	new_header = ['chrom']
	for colname in header:
		if colname in names_dict:
			print(names_dict[colname])
			new_header.append(names_dict[colname])
	with gz.open(args.outfile, 'wt', encoding='utf-8') as out:
		out.write(" ".join(new_header) + "\n") 
		#next(target)
		for line in target:
			out.write(line)
      

  
