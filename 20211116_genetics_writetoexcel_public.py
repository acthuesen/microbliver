# this script writes genetic data (SNP and riskscores) to a single excelfile, which also includes a metadata sheet
# dependencies
import pandas as pd
import os 
import glob
import re
import openpyxl as op
import shutil

# data files must adhere to the naming convention [datatype]_[cohort] - ensure that this happens in split by cohort step (set up this way if you use ACs split by cohort scripts)
# only csv files are included in the filelist (again, set up this way if you use ACs split by cohort scripts)

# in the following, replace PATH with path
# Input files
filedir = 'PATH\\Genetics\\SNPs_PRSs_for_datahub\\preprocessing\\split\\'

# this is the list of input files (all csv)
filenames = glob.glob(os.path.join(filedir, '*.csv'))

# this is the metadata
metadata = 'PATH\\Genetics\\SNPs_PRSs_for_datahub\\preprocessing\\metadata\\metadata.xlsx'

# Output dir
outputdir = 'PATH\\Genetics\\SNPs_PRSs_for_datahub\\'

# first, initialize excel files
for filename in filenames:
    cohort = re.split('[_.\\\\]',filename)[-2] # find cohort from filename because it's part of the file extension
    extension = "geneticdata_"+cohort+".xlsx" # create the fileextension that the file will be saved in
    shutil.copy(metadata, os.path.join(outputdir,extension)) # copy to a fresh excel
    # rest is cleanup
    del cohort
    del extension

# then write data to excelfiles
for filename in filenames:
    df = pd.read_csv(os.path.join(filedir,filename)) # read in csv file that holds the data
    df_type = re.split('[_.\\\\]',filename)[-3] # define the type of data in file by input filename (for sheetname)
    cohort = re.split('[_.\\\\]',filename)[-2] # find cohort from input filename (for output filename)
    extension = "geneticdata_"+cohort+".xlsx" # create the fileextension of the output filename
    with pd.ExcelWriter(os.path.join(outputdir,extension), engine='openpyxl', mode='a') as writer: # call the excelwriter
        df.to_excel(writer, index=False, sheet_name=df_type) # write data
    # rest is cleanup
    del df
    del df_type
    del cohort
    del extension

