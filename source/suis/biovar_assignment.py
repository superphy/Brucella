import pandas as pd 

metadata = pd.read_csv("Metadata.csv", index_col = 'Sample' )
bpm = pd.read_csv('suis_biovar_files/Blast_Primer_Matching.csv', converters={'File ID': lambda x: str(x)})
files = list(bpm['File ID'].unique())

sample_matching = {} # key = shortened sample name, value = full sample name
for sample in metadata.index.values:
	sample_matching[sample[4:13]]=sample

biovar_1 = [773, 774, 775, 424, 425, 426, 196, 197, 198]
biovar_2 = [773, 774, 775, 550, 551, 552, 277, 278, 279]
biovar_3 = [773, 774, 775, 298, 299, 300, 196, 197, 198]
biovar_4 = [773, 774, 775, 613, 614, 615, 196, 197, 198]
biovar_5 = [773, 613, 277, 196]

def biovar(file):
	file_df = bpm[bpm['File ID']==file].fillna(0)
	p1 = file_df['Primer One (BMEI1426-BMEI1427)'].to_list() # amplicons derived from primer 1
	p2 = file_df['Primer Two (BR1080f-BR1080r)'].to_list() # amplicons derived from primer 2
	p3 = file_df['Primer Three (BMEI1688-BMEI1687)'].to_list() # amplicons derived from primer 3
	p4 = file_df['Primer Four (BMEI0205f-BMEI0205r)'].to_list() # amplicons derived from primer 4

	amplicons = p1+p2+p3+p4
	amplicons = [int(amp) for amp in amplicons if amp != 0] # removing the nan values 
	amplicons = list(set(amplicons))

	if len(amplicons)<3: # every biovar class requires at least three primer matches, therefore if there are less then three amplicons
		return('No Matched Biovar')
	else:
		bv1 = all(amps in biovar_1 for amps in amplicons)
		bv2 = all(amps in biovar_2 for amps in amplicons)
		bv3 = all(amps in biovar_3 for amps in amplicons)
		bv4 = all(amps in biovar_4 for amps in amplicons)
		bv5 = all(amps in amplicons for amps in biovar_5)

		if bv1: # if all of the amplicons found in the given file are found in the biovar 1 criteria
			return('Biovar 1')
		if bv2: # if all of the amplicons returned for one sample are found in the biovar 2 criteria
			return('Biovar 2')
		if bv3: # if all of the amplicons returned for one sample are found in the biovar 3 criteria
			return('Biovar 3')
		if bv5: # if all of the amplicons returned for one sample are found in the biovar 5 criteria
			return('Biovar 5')
		if bv4: # if all of the amplicons returned for one sample are found in the biovar 4 criteria
			return('Biovar 4')
		else: 
			return('No Matched Biovar')

df = pd.DataFrame(index = ['Biovar Classification'])
for file in files:
	df[sample_matching[file]] = [biovar(file)]

df = df.T
df.to_csv('suis_biovar_files/Suis_Biovar.csv')