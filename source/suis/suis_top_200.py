import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

ranks_df = pd.read_csv('../../suis_biovar_files/Suis_Ranks.csv', index_col = 0)
ranks_df = ranks_df.filter(axis = 1, like='RANK')
ranks_index = ranks_df.index.values
refrence_genomes = pd.read_csv('../../suis_biovar_files/Suis_Refrence_Genomes.csv')
k_loc = pd.read_csv('../../suis_biovar_files/Suis_Kmer_Locations.csv', index_col = 0)
suis_biovar = pd.read_csv('../../suis_biovar_files/Suis_Biovar.csv')
suis_biovar = suis_biovar[suis_biovar['Biovar Classification'] != 'No Matched Biovar']
biovars = list(suis_biovar['Biovar Classification'].unique())

chroms = {'C1':'C2', 'C2':'C1'}
in_list = [(x,y) for y in biovars for x in ['self', 'alt']]

def rolling(rank_df): # takes the rolling average of the rank values sorted by location
	roll_mean = rank_df.rolling(window=200).mean()
	return(list(roll_mean))

def biovar_df_f(biovar, soa): # returns a biovar specific df with the rank of each kmer and their location on either chromosome
	biovar_df = pd.DataFrame(index = ranks_index)
	if soa == 'self':
		biovar_df['Rank'] = ranks_df[biovar+' RANK']
		biovar_df['C1 Loc'] = k_loc['C1 '+biovar+' self']
		biovar_df['C2 Loc'] = k_loc['C2 '+biovar+' self']
		biovar_df = biovar_df.dropna(axis = 0, thresh=2) # drop any rows where the kmer isnt found on any chromosome
		return(biovar_df)
	if soa == 'alt':
		biovar_df['Rank'] = ranks_df[biovar+' RANK']
		biovar_df['C1 Loc'] = k_loc['C1 '+biovar+' alt']
		biovar_df['C2 Loc'] = k_loc['C2 '+biovar+' alt']	
		biovar_df = biovar_df.dropna(axis = 0, thresh=2) # drop any rows where the kmer isnt found on any chromosome
		return(biovar_df)

def loc_dict_f(biovar_df): # generates an empty dictionary to set rank values for locations without a match to zero
	loc_dict = {}
	locs = list(biovar_df['C1 Loc'].dropna())+list(biovar_df['C2 Loc'].dropna())
	for i in range(0, int(max(locs))): # zero to the largest number found as a location on either C1 or C2 
		loc_dict[i] = (0, np.inf)
	return(loc_dict)

def top_row(biovar, chrom, soa):
	biovar_df= biovar_df_f(biovar, soa)
	loc_dict = loc_dict_f(biovar_df)
	 
	biovar_df = biovar_df.drop(columns=chroms[chrom]+' Loc').dropna(axis =0, how='any')
	for index, row in biovar_df.iterrows():
		loc_dict[row[chrom+' Loc']]=(row['Rank'], row.name) # key = location, value = (rank, kmer)
	df = pd.DataFrame()
	
	kmers = [value[1] for value in loc_dict.values()] # pulling the kmers from the dict
	ranks = [value[0] for value in loc_dict.values()] # pulling the rank values from the dict
	
	df['Location'] = loc_dict.keys() # adding location column to df 
	df['Kmer'] = kmers # adding kmer column to df 
	df['Ranks'] = ranks # adding ranks column to df 

	roll_mean= rolling(df['Ranks'])
	df['Rolling Avg'] = roll_mean # adding rolling average column to df 
	
	if soa =='self':
		df = df.sort_values(by=['Rolling Avg'], ascending=False).head(1) # returns the row in the df with the largest rolling average value
		return(df, loc_dict)
	if soa == 'alt':
		#MAKE SURE THIS DOES WHAT U THINK IT DOES!
		df = df.sort_values(by=['Rolling Avg'], ascending=True).head(1) # returns the row in the df with the smallest rolling average value
		return(df, loc_dict)

df, loc_dict = top_row('Biovar 3', 'C2', 'self')
print(df)

def intermediate(in_list):
	soa, biovar = in_list
	c1_top_row, c1_loc_dict = top_row(biovar, 'C1', soa)
	c2_top_row, c2_loc_dict = top_row(biovar, 'C2', soa)

	ref_genome = refrence_genomes[refrence_genomes['biovar']==biovar]
	ref_genome = list(ref_genome['Sample'])
	ref_genome = ref_genome[0]

	max_ra = max(abs(float(c1_top_row['Rolling Avg'])), abs(float(c2_top_row['Rolling Avg']))) # is the rolling average on C1 or C2 larger? 
	if max_ra == abs(float(c1_top_row['Rolling Avg'])): # if C1 is larger
		# the list returned to be the col in the df [chrom. start loc, end loc, rank, rolling avg, ref genome]
		out_col = ['C1',int(c1_top_row['Location'])-200, int(c1_top_row['Location']), float(c1_top_row['Ranks']), float(c1_top_row['Rolling Avg']), ref_genome]
		# lists used to plot the rank values from the start loc to the end loc of the top 200 bp
		locations = [i for i in range(int(c1_top_row['Location'])-200, int(c1_top_row['Location']))]
		ranks = [c1_loc_dict[loc][0] for loc in locations]
		return(biovar, out_col, locations, ranks, 'C1', ref_genome, soa)
	
	if max_ra == abs(float(c2_top_row['Rolling Avg'])): #If C2 is larger
		# the list returned to be the col in the df [chrom. start loc, end loc, rank, rolling avg, ref genome]
		out_col = ['C2', int(c2_top_row['Location'])-200, int(c2_top_row['Location']), float(c2_top_row['Ranks']), float(c2_top_row['Rolling Avg']), ref_genome]
		# lists used to plot the rank values from the start loc to the end loc of the top 200 bp
		locations =[i for i in range(int(c2_top_row['Location'])-200, int(c2_top_row['Location']))]
		ranks = [c2_loc_dict[loc][0] for loc in locations]
		return(biovar, out_col, locations, ranks, 'C2', ref_genome, soa)

'''
final_df = pd.DataFrame(index = ['Chromosome','Start Location', 'End Location', 'Rank', 'Rolling Avg', 'Refrence Genome'])

with ProcessPoolExecutor(cpu_count()-10) as ppe:
	for biovar, column, locations, ranks, chrom, ref_genome, soa in ppe.map(intermediate, in_list):
		final_df[biovar+" "+soa] = column # adding a column to the df for the given species

		# generating figure to show the rank values across the top 200 bp
		plt.figure(figsize = (17,7))
		plt.plot(locations, ranks, color=col_dict[biovar])
		plt.xlabel('Location')
		plt.ylabel('Rank')
		plt.title('Rank Distribution for Top 200 bp on '+chrom+' of '+biovar+' '+soa+'\n'+ref_genome)

		plt.savefig('Visualizations/Top_200_bp/'+biovar+'_'+soa+'.png')
		plt.clf()

final_df.to_csv('Top_200_bp.csv')
'''