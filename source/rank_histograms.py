import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

ranks = pd.read_csv('Ranks.csv', index_col = 0)
strain_occurances = pd.read_csv("Species_Occurrence.csv", index_col = 0)
species = list(strain_occurances.index.values )
col_dict = {'Brucella sp':'#ebc04b', 'Brucella suis':'#f8a4d0','Brucella abortus': '#72a24a','Brucella canis': '#b8642f','Brucella ceti':'#59cf9e', 'Brucella melitensis':'#e8382e', 'Brucella neotomae':'#6c86e2', 'Brucella ovis':'#ba6eb4', 'Brucella pinnipedialis':'#f4792a'}

bins = range(-10,12,1)
bins = [x/10 for x in bins]
plt.xlabel('Rank Value')
plt.ylabel('Number of Kmers')

def hist(species):
		sp_rank = ranks[species+' RANK'].to_list()
		histo = plt.hist(sp_rank, bins=bins, color = col_dict[species], rwidth = 0.5)
		plt.title(species)
		plt.savefig('Visualizations/Rank_Histograms/'+species+'.png')

with ProcessPoolExecutor(9) as ppe:
	for plt in ppe.map(hist, species):
		x=0
