import argparse as ap
parser = ap.ArgumentParser()

parser.add_argument('target_counts', help = 'file path of targe counts file generated in between leafcutter scripts')
parser.add_argument('associated_names', help = 'the file E-GEUV-1.sdrf.txt containing all the name associations') 
parser.add_argument('outfile', help = file path of the renamed headers file)
args = parser.parse_args()

names_dict = {}
with open(args.associated_names, 'r') as names:
  for line in names:
    fields = line.strip('\n').split('\t')
    names_dict[fields[27] + ".staraligned "] = fields[0]
    
with open(args.target_counts, 'r') as target:
	header = target.readline().strip('\n').split('\t')
  for colname in header:
    if colname in names_dict:
      colname = dict[colname]
  with open(args.outfile, 'w') as out:
    out.write("\t".join(header) + "\n") 
    #next(target)
    for line in target:
      out.write(line)
      

  
