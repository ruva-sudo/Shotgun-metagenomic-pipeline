#!/usr/bin/env python

#This script is for filtering bins based on percentage completeness and genome size
#Note that it may not be as computationally efficient for samples that generate large number of bins i.e >= 500

import os
import re
import shutil
import pandas as pd

path = os.listdir()
for files in path:
    match = re.search("\.summary$", files)
    if match:
        print(files)
        df = pd.read_csv(files, sep="\t")

def filter_bins(dataframe):
    dataframe.rename(columns = {'Bin name':'Bin_name', 'Genome size':'Genome_size', 'GC content':'GC_content'}, inplace = True)
    dataframe['Completeness'] = dataframe['Completeness'].str.replace(r'%', '').astype(float)
    final_df = dataframe[(dataframe.Completeness >= 30) & (dataframe.Genome_size >= 100000)]
    x = final_df['Bin_name']
    new_dir = 'passed_bins'
    os.mkdir(new_dir)
    for names in x:
        new_path = new_dir
        shutil.move(names, new_path)
filter_bins(df)

directory = 'failed_bins'
os.mkdir(directory)
failed_path = directory
filtered_path = os.listdir()
for fastas in filtered_path:
    failed_match = re.search("\.fasta$", fastas)
    if failed_match:
        shutil.move(fastas, failed_path)

