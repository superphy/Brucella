import pandas as pd 
import re, sys
import scipy.stats as stats
import time 
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

if __name__ == '__main__':
	
	start_time = time.time()
	metadata = pd.read_csv("Metadata.csv", index_col = 'Sample')
	strain_occurances = pd.read_csv("Species_Occurrence.csv")
	kmer_count = pd.read_csv("Kmer_Counts.csv", index_col = 'Unnamed: 0')
	
	print('Dataframes Loaded')

	species_occ_dict = {} # key = species, value = # of samples of that species
	for i in range(0,len(strain_occurances)):
		species_occ_dict[strain_occurances.iat[i,0]] = strain_occurances["Species Occurrence"].iat[i]

	species_gen_dict = {} # key = genome file name, value = species
	for j in range(0,len(metadata)):
		species_gen_dict[metadata.index.values[j]] = metadata["Strain"].iat[j]

	rows = []
	rank_dict = {}
	p_dict = {}
	for k in range(0,len(kmer_count)): # for each kmer in the dataframe
		row = kmer_count.iloc[k]
		rows.append(row)
		rank_dict[row.name] = 0
		p_dict[row.name] = 0

	print('Dictionaries Created')

	cpu_count = cpu_count()
	rank_cpu = len(set(species_occ_dict))
	sp_match_cpu = cpu_count - rank_cpu - 5

	def sp_match(row): # for a given kmer
		row , species = row
		sp_present = 0 # counter for the amount of times a kmer is present in a sequence of the same species as the function input 
		nsp_present = 0 # counter for the amount of times a kmer is present in a sequence of a different species besides the function input
		kmer = row.name
		for l in range(0,len(row)):	# for every entry in the row
			if row[l] != 0: 
				sample = kmer_count.columns[l]
				sample = sample[0:len(sample)-4]
				if species_gen_dict[sample] == species:
					sp_present +=1
				else: 
					nsp_present+=1 
		return(kmer, sp_present, nsp_present)

	def rank(species):
		for key in rank_dict: 
			rank_dict[key] = 0
		for key in p_dict: 
			p_dict[key] = 0
		divisor_sp = species_occ_dict[species] # number of samples with a species matching the input species
		divisor_nsp = len(kmer_count.columns) - divisor_sp # number of samples with a species not matching the input species
		zip_list = zip(rows,[species for i in rows])
		count = 0 
		inter_count = 0
		print(species, inter_count)
		with ProcessPoolExecutor(sp_match_cpu) as ppe:
			for kmer, sp_present, nsp_present in ppe.map(sp_match, zip_list):
				rank = (sp_present/divisor_sp)-(nsp_present/divisor_nsp)
				odds_ratio, p_value = stats.fisher_exact([[sp_present, nsp_present],[divisor_sp, divisor_nsp]])
				p_dict[kmer]=p_value
				rank_dict[kmer]=rank
				count+=1
				if count == 3000:
					inter_count+=count
					count=0
					print(species, inter_count)
				if inter_count+count == len(rows):
					print(species, 'complete')
		return(rank_dict,p_dict, species)

	df = pd.DataFrame(index = kmer_count.index.values)

	with ProcessPoolExecutor(rank_cpu) as ppe:
		for rank_dict, p_dict, species in ppe.map(rank, set(species_occ_dict)):
			df[species+' RANK'] = rank_dict.values()
			df[species+' P_VAL'] = p_dict.values()
	
	df.to_csv('Ranks.csv')
	end_time = time.time()
	print('Time Elapsed: ', end_time-start_time)
	