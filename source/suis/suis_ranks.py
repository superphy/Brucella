import pandas as pd 
import re, sys
import scipy.stats as stats
import time 
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

if __name__ == '__main__':
	
	suis_biovar = pd.read_csv('Suis_Biovar.csv', converters={'Unnamed: 0': lambda x: str(x)})
	suis_biovar = suis_biovar.set_index('Unnamed: 0')
	suis_biovar = suis_biovar[suis_biovar['Biovar Classification'] != 'No Matched Biovar']
	kmer_count = pd.read_csv("suis_kmer_counts.csv", index_col = 'Unnamed: 0')

	biovar_occ_dict = {} # key = biovar, value = # of samples of that biovar
	for biovar in list(suis_biovar['Biovar Classification'].unique()):
		biovar_occ_dict[biovar]=len(suis_biovar[suis_biovar['Biovar Classification']== biovar])

	biovar_gen_dict = {} # key = genome file name, value = biovar
	for i, row in enumerate(suis_biovar.index.values):
		biovar_gen_dict[row]=suis_biovar['Biovar Classification'].iat[i]

	rows = []
	rank_dict = {}
	p_dict = {}
	for k in range(0,len(kmer_count)): # for each kmer in the dataframe
		row = kmer_count.iloc[k]
		rows.append(row)
		rank_dict[row.name] = 0
		p_dict[row.name] = 0

	cpu_count = cpu_count()
	rank_cpu = len(set(biovar_occ_dict))
	sp_match_cpu = cpu_count - rank_cpu - 5

	def sp_match(row): # for a given kmer
		row , biovar = row
		sp_present = 0 # counter for the amount of times a kmer is present in a sequence of the same biovar as the function input 
		nsp_present = 0 # counter for the amount of times a kmer is present in a sequence of a different biovar besides the function input
		kmer = row.name
		for l in range(0,len(row)):	# for every entry in the row
			if row[l] != 0: 
				sample = kmer_count.columns[l]
				sample = sample[0:len(sample)-4]
				if biovar_gen_dict[sample] == biovar:
					sp_present +=1
				else: 
					nsp_present+=1 
		return(kmer, sp_present, nsp_present)

	def rank(biovar):
		for key in rank_dict: 
			rank_dict[key] = 0
		for key in p_dict: 
			p_dict[key] = 0
		divisor_sp = biovar_occ_dict[biovar] # number of samples with a biovar matching the input biovar
		divisor_nsp = len(kmer_count.columns) - divisor_sp # number of samples with a biovar not matching the input biovar
		zip_list = zip(rows,[biovar for i in rows])
		with ProcessPoolExecutor(sp_match_cpu) as ppe:
			for kmer, sp_present, nsp_present in ppe.map(sp_match, zip_list):
				rank = (sp_present/divisor_sp)-(nsp_present/divisor_nsp)
				odds_ratio, p_value = stats.fisher_exact([[sp_present, nsp_present],[divisor_sp, divisor_nsp]])
				p_dict[kmer]=p_value
				rank_dict[kmer]=rank
		return(rank_dict,p_dict, biovar)

	df = pd.DataFrame(index = kmer_count.index.values)

	with ProcessPoolExecutor(rank_cpu) as ppe:
		for rank_dict, p_dict, biovar in ppe.map(rank, set(biovar_occ_dict)):
			df[biovar+' RANK'] = rank_dict.values()
			df[biovar+' P_VAL'] = p_dict.values()
	
	df.to_csv('Suis_Ranks.csv')
