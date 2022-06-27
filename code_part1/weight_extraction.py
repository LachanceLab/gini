#! /usr/bin/env python3

#Example input: ./weight_extraction.py -s /Volumes/T7/prive_sf_211_PLR_9ancestries/153-1KG_PLR.txt -w /Users/adrianharris/Desktop/

import pandas as pd
import argparse as ap
import os
import time

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description="creation of SNP.txt for each trait")
parser.add_argument("-s", "--summary", help="path to summary file", required=True)
parser.add_argument("-w", "--working", help="working directory where SNP.txt files will be stored", required=True)
args = parser.parse_args()

#Load into pandas dataframe and sort
summary_df = pd.read_csv(args.summary, delimiter = '\t')
summary_df['effect_weight'] = summary_df['effect_weight'].apply(lambda x: abs(x)) #calculating absolute value for the weight 
summary_df.sort_values(['effect_weight'], ascending=False, inplace=True) #ordering by abs(effect weight)
#print(summary_df)

#Appending to rsList
rsList = []
for index in summary_df.index[0:100]:
	rsList.append((summary_df['rsID'][index], summary_df['A2'][index])) #Tuple to be used in the creation of snp_list.txt

#print(rsList)

os.chdir(args.working)

#Make directory if does not exists
if not os.path.exists('snp_list_dir'):
	os.makedirs('snp_list_dir')
else:
	pass

os.chdir('snp_list_dir') #move to snp_list_dir

#Creation of snp_list filename
count = 0
new_filename = []
filename = args.summary
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


