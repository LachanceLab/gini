#! /usr/bin/env python3

#Example input: ./calc_effect_dose.py -s <summary_files_directory> -g <directory_with_genotype_matrices> -o <output_directory>

import pandas as pd
import argparse as ap
import os
import time

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Calculate effect dose for each UKBB individual using the matrix')
parser.add_argument("-s", "--summary_dir", help="directory of summary files", required=True)
parser.add_argument("-g", "--genotype_dir", help="directory of genotype matrices for all traits", required=True)
parser.add_argument("-o", "--out_dir", help="directory of genotype matrices for all traits", required=True)
args = parser.parse_args()

def calc_effect(summary):
	#Load into pandas dataframe and sort
	summary_df = pd.read_csv(summary, delimiter = '\t')
	summary_df['abs_effect_weight'] = summary_df['effect_weight'].apply(lambda x: abs(x)) #addition of new_column to sort by while retaining info of 'effect_weight' col
	summary_df.sort_values(['abs_effect_weight'], ascending=False, inplace=True)
	#print(summary_df)

	#Creation of dictionary for each chromosome
	chrDictionary = {}
	count = 1
	while count < 23: #change back to 23 once done
		chrDictionary[count] = {} #creation of nested-dictionary within dictionary
		count += 1

	for index in summary_df.index[0:100]:
		chromosome = summary_df['chrom'][index] #chromosome number
		rsID_allele = str(summary_df['rsID'][index]) + '_' + str(summary_df['A2'][index])
		chrDictionary[chromosome][rsID_allele] = summary_df['effect_weight'][index]#dictionary with rs_A/T/G/C and effect_size (key:value)

	#Parsing out the directory name in preparation for directory change soon
	dirname = summary
	dirname=dirname.split('/')
	dirname=dirname[-1]
	dirname=''.join(dirname)
	trait=dirname[:-4] #creation of variable trait used for final output csv
	dirname=args.genotype_dir + trait

	os.chdir(dirname) #change to trait-specific directory

	#Calculating effect dosages per individual
	final_df = pd.DataFrame() #empty dataframe
	for chrom in chrDictionary.keys(): 
		inputFile = 'genotype_matrix_chr' + str(chrom) + '.raw'
		if os.path.exists(inputFile):
			matrix_df = pd.read_csv(inputFile, delimiter = ' ')
			matrix_df = matrix_df.fillna(0) #replace NA with 0
			for key in chrDictionary[chrom].keys(): #iterate through rs_A/T/G/C
				matrix_df[key] = matrix_df[key].apply(lambda x: x*(chrDictionary[chrom][key])) #multiply current dataframe by effect sizes
			if final_df.empty:
				matrix_df = matrix_df.drop(['PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis=1)
				final_df = matrix_df #overwrite the empty one
			else:
				matrix_df = matrix_df.drop(['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis=1)
				final_df = pd.concat([final_df, matrix_df], axis=1)	
		else:
			pass

	#print(final_df)

	#Make directory if does not exists
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)
	else:
		pass

	os.chdir(args.out_dir)

	ouptut_name = trait + '.csv'
	final_df.to_csv(ouptut_name)

#List of files to run
directory = args.summary_dir
f = []
for filename in os.listdir(directory):
    file = os.path.join(directory, filename)
    f.append(file)

for file in f:
	print(f"Running file {file}")
	calc_effect(file)
