import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

suis_ranks = pd.read_csv('suis_biovar_files/Suis_Ranks.csv')
suis_bv = pd.read_csv('suis_biovar_files/Suis_Biovar.csv')
suis_bv = suis_bv[suis_bv['Biovar Classification'] != 'No Matched Biovar']
biovar = list(suis_bv['Biovar Classification'].unique())


bins = range(-10,12,1)
bins = [x/10 for x in bins]
plt.xlabel('Rank Value')
plt.ylabel('Number of Kmers')

def hist(species):
	if 'Biovar' in species:
		suis_sp_rank = suis_ranks[species+' RANK'].to_list()
		histo = plt.hist(suis_sp_rank, bins=bins, rwidth = 0.5, color = '#000000')
		plt.title('Suis '+species)
		plt.savefig('Visualizations/Rank_Histograms/Suis '+species+'.png')

with ProcessPoolExecutor(9) as ppe:
	for plt in ppe.map(hist, biovar):
		x=0
