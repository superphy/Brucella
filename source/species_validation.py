import os 
import pandas as pd
from Bio import SeqIO
import re 

j=0

rmash = pd.read_csv("rmash_files/rmash_filtered_output.csv")
samples = list(rmash['Sample'])
fna_sp_col = []
rec_id_col = []

#gathering metadata information and fna Strain info
for sample in samples:
	location = "refseq/bacteria/"+sample[0:15]+"/"+sample+".fna"
	for seq_record in SeqIO.parse(location, "fasta"):
		description = seq_record.description
		match= re.search('Brucella\s\w+',description)
		if match!= None:
			fna_strain = match.group(0)
		if match== None:
			fna_strain = "no strain found"	
	rec_id_col.append(sample[4:13])
	fna_sp_col.append(fna_strain)	


# adding a Strain and Record ID column to the dataframe 
rmash.insert(loc = 3, column='Strain', value = fna_sp_col)
rmash.insert(loc = 2, column= 'Record ID', value = rec_id_col)

# comapring the rmash results and the fna Strain
while j<len(rmash):
	sample = rmash['Sample'].iat[j]
	# if the fna Strain is sp, it takes on the Strain found by rmash
	if rmash['Strain'].iat[j] == 'Brucella sp':
		if rmash['Max Distance'].iat[j] < 0.001:
			rmash['Strain'].iat[j] = rmash['Rmash Species'].iat[j]
		else:
			rmash['Rmash Species'].iat[j] = rmash['Strain'].iat[j]
	# if the fna and rmash Strain info does not match the file is deleted  
	if rmash['Strain'].iat[j] != rmash['Rmash Species'].iat[j]:
		sysin = 'rm -rf refseq/bacteria/'+sample[0:15]
		os.system(sysin)
	j+=1 

#deleting any entries where the fna and rmash Strain info does not match
rmash = rmash[rmash['Strain']==rmash['Rmash Species']]

# deleting the Rmash column and creating metadata.csv 
rmash.drop(['Rmash Species'], axis = 1, inplace = True)
rmash.to_csv('rmash_files/Validated_Species.csv', index = False)