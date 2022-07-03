#! /usr/bin/env python3

#Example input: ./extract_weights.py -s <summary_dir> -o <output_dir>

import pandas as pd
import argparse as ap
import os
import time

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description="creation of SNP.txt for each trait")
parser.add_argument("-s", "--summary_dir", help="directroy with summary files", required=True)
parser.add_argument("-o", "--out_dir", help="output directory where SNP.txt files will be stored", required=True)
args = parser.parse_args()


def create_snp_files(summary):
	#Load into pandas dataframe and sort
	summary_df = pd.read_csv(summary, delimiter = '\t')
	summary_df['effect_weight'] = summary_df['effect_weight'].apply(lambda x: abs(x)) #calculating absolute value for the weight 
	summary_df.sort_values(['effect_weight'], ascending=False, inplace=True) #ordering by abs(effect weight)
	#print(summary_df)

	#Appending to rsList
	rsList = []
	for index in summary_df.index[0:100]:
		rsList.append((summary_df['rsID'][index], summary_df['A2'][index])) #Tuple to be used in the creation of snp_list.txt

	#print(rsList)

	#Make directory if does not exists
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)
	else:
		pass

	os.chdir(args.out_dir) #move to snp_list_dir

	#Creation of snp_list filename
	count = 0
	new_filename = []
	filename = summary
	for i in range(len(filename) - 1, -1, -1):
		if filename[i] == '.':
			count = 1
			new_filename.append(filename[i])
		elif filename[i] == '/':
			break
		elif count == 1:
			new_filename.append(filename[i])
		else:
			pass

	new_filename.reverse()
	new_filename = new_filename[:-1] #removal of dot on the end
	new_filename = ''.join(new_filename)
	new_filename += '_snp_list.txt'

	#Writing to output
	output = open(new_filename, 'w')
	for item in rsList:
		item = '\t'.join(item)
		output.write(item)
		output.write('\n')
	output.close()

#Generating list of all summary files 
directory = args.summary_dir
f = []
for filename in os.listdir(directory):
    file = os.path.join(directory, filename)
    f.append(file)

#Creating SNP files from summary files
for file in f:
	print(f"Running file {file}")
	create_snp_files(file)


