#! /usr/bin/env python3

#Example input: ./downsample.py -p <file_with_all_pop_iids_path> -f <up-to-date_csv_of_ukbb_iids> -o <ouptut_filename>

import pandas as pd
import argparse as ap
import os
import time

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description="randomly select 1224 individuals for each ancestry")
parser.add_argument("-p", "--pop", help="path to population iids file", required=True)
parser.add_argument("-f", "--filtered", help="most up-to-date list of individuals (iids) who still consent in UKBB", required=True)
parser.add_argument("-o", "--output", help="output list (csv) of 1224 sampled individuals for each population", required=True)
args = parser.parse_args()

pop_df = pd.read_csv(args.pop, delimiter = '\t') 
#print(pop_df)

filtered_df = pd.read_csv(args.filtered) #most up-to-date list of individuals who still consent in UKBB 

#merge filtered_iids.csv with pops
merged_df = pd.merge(pop_df, filtered_df, on='IID') #merging on the filtered_df to remove those that have opted out of UKBB

#merged_df.to_csv('filtered_pops.csv')

#Subset the merged_df for each ancestry and randomly sample 1224 individuals
uk_df = merged_df.loc[(merged_df['ancestry'] == 'United Kingdom')&(merged_df['white_brit'] == 1)] 
final_df = uk_df.sample(n = 1224)

italy_df = merged_df.loc[merged_df['ancestry'] == 'Italy']
italy_df = italy_df.sample(n = 1224)
final_df = final_df.append(italy_df, ignore_index=True)

india_df = merged_df.loc[merged_df['ancestry'] == 'India']
india_df = india_df.sample(n = 1224)
final_df = final_df.append(india_df, ignore_index=True)

poland_df = merged_df.loc[merged_df['ancestry'] == 'Poland']
poland_df = poland_df.sample(n = 1224)
final_df = final_df.append(poland_df, ignore_index=True)

nigeria_df = merged_df.loc[merged_df['ancestry'] == 'Nigeria']
nigeria_df = nigeria_df.sample(n = 1224)
final_df = final_df.append(nigeria_df, ignore_index=True)

caribbean_df = merged_df.loc[merged_df['ancestry'] == 'Caribbean']
caribbean_df = caribbean_df.sample(n = 1224)
final_df = final_df.append(caribbean_df, ignore_index=True)

ashkenazi_df = merged_df.loc[merged_df['ancestry'] == 'Ashkenazi']
ashkenazi_df = ashkenazi_df.sample(n = 1224)
final_df = final_df.append(ashkenazi_df, ignore_index=True)

china_df = merged_df.loc[merged_df['ancestry'] == 'China']
china_df = china_df.sample(n = 1224)
final_df = final_df.append(china_df, ignore_index=True)

iran_df = merged_df.loc[merged_df['ancestry'] == 'Iran']
final_df = final_df.append(iran_df, ignore_index=True)

#print(final_df)

csv_ouptut = args.output + '.csv'
final_df.to_csv(csv_ouptut, index=None)

final_df = final_df[['IID']]
final_df['FID'] = final_df.loc[:, 'IID']

txt_output = args.output + '.txt'
final_df.to_csv(txt_output, header=None, index=None, sep='\t')

#Formation of two files: (1) IID, ancestry, white_brit, unamed (2) IID and FID

