import pandas as pd 
import numpy as np
import ahocorasick
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.SeqIO.FastaIO import SimpleFastaParser
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

metadata = pd.read_csv('Metadata.csv')
ranks = pd.read_csv('Ranks.csv', index_col = 0)
strain_occurances = pd.read_csv("Strain_Occurances.csv", index_col = 0)
species = strain_occurances.index.values 
col_dict = {'Brucella sp':'#ebc04b', 'Brucella suis':'#f8a4d0','Brucella abortus': '#72a24a','Brucella canis': '#b8642f','Brucella ceti':'#59cf9e', 'Brucella melitensis':'#e8382e', 'Brucella neotomae':'#6c86e2', 'Brucella ovis':'#ba6eb4', 'Brucella pinnipedialis':'#f4792a' }
meta_ref = metadata[metadata['Number of Contigs'] == 2] # only complete sequences (a.k.a 2 contigs)
meta_ref = meta_ref.drop_duplicates(subset = 'Species', keep='first') # one of each species
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

c1_dict, c2_dict = {}, {}
for kmer in ranks.index: 
	c1_dict[kmer] = np.inf
	c2_dict[kmer] = np.inf

#Choosing and locating a refrence genome 
def refrence_location(species):
	ref_row = meta_ref[meta_ref['Species']==species]
	ref_file = ref_row['Sample'].to_list()
	ref_file_location = 'Approved_Sequences/'+ref_file[0]+'.fna'
	return(ref_file_location)

#Creates 2 haystacks from the refrence genome, one for each chromosome
def haystacks(species):
	location = refrence_location(species)
	sequence = []
	with open(location) as handle:
		for title, seq in SimpleFastaParser(handle):
			sequence.append(seq)
	chromosome1 = sequence[0]
	chromosome2 = sequence[1]
	return(chromosome1, chromosome2)

#Generates the reverse complement of a given kmer
def reverse_complement(kmer):
	return ''.join([complement[base] for base in kmer[::-1]])

#Creates the needles (all kmers)
def needles(species):
	kmers = ranks.index.values
	rank = ranks[species+' RANK']
	aho = ahocorasick.Automaton()
	for i in range(0,len(kmers)):
		kmer = kmers[i]
		rev_c = reverse_complement(kmer)
		aho.add_word(kmer,(kmer,rank[kmer]))
		aho.add_word(rev_c,(kmer,rank[kmer]))
	aho.make_automaton()
	return(aho)

#Runs the ahocorasick algorithm to find the location of the kmers (needles) in the refrence genome (haystack) 
def find_needles(species):
	aho = needles(species)
	chromosome1, chromosome2 = haystacks(species)
	c1_index, c1_rank, c2_index, c2_rank = [], [], [], []
	for end_index, (kmer, rank) in aho.iter(chromosome1):
		c1_index.append(end_index)
		c1_dict[kmer] = end_index
		c1_rank.append(rank)
	for end_index, (kmer, rank) in aho.iter(chromosome2):
		c2_index.append(end_index)
		c2_dict[kmer] = end_index
		c2_rank.append(rank)
	return(species, c1_index, c1_rank, c2_index, c2_rank, c1_dict, c2_dict)

df = pd.DataFrame(index = ranks.index)

with ProcessPoolExecutor(len(species)) as ppe:
	for species, c1_index, c1_rank, c2_index, c2_rank, c1_dict, c2_dict in ppe.map(find_needles, species):
		#For any kmer that is not found in the refrence genome it is assigned a value of zero
		c1_index_padded = range(0, max(c1_index))
		c1_rank_padded = np.zeros(max(c1_index))

		for i in range(0,len(c1_rank)):
			index = c1_index[i]
			c1_rank_padded[index-1]=c1_rank[i]

		c2_index_padded = range(0, max(c2_index))
		c2_rank_padded = np.zeros(max(c2_index))

		for i in range(0,len(c2_rank)):
			index = c2_index[i]
			c2_rank_padded[index-1]=c2_rank[i]
		
		#Creating the dataframe of kmer locations
		df['C1 '+species] = c1_dict.values()
		df['C2 '+species] = c2_dict.values()

		#Generating the moving average of the rank values
		df_c1_rank = pd.DataFrame(c1_rank_padded)
		roll_c1 = df_c1_rank.rolling(window=1000).mean()
		df_c2_rank = pd.DataFrame(c2_rank_padded)
		roll_c2 = df_c2_rank.rolling(window=1000).mean()

		#Generating plots
		fig, (ax1, ax2) = plt.subplots(1,2, sharey = True, figsize=(17,7))
		fig.suptitle(species)
		
		ax1.plot(c1_index_padded, c1_rank_padded, color= '#dedcdf', linewidth=0.25, ls = ':')
		ax1.plot(c1_index_padded, roll_c1, color=col_dict[species], linewidth=1)
		ax1.set_ylim([-1,1])
		ax1.set_title('Chromosome 1')
		
		ax2.plot(c2_index_padded, c2_rank_padded, color='#dedcdf', linewidth=0.25, ls = ':')
		ax2.plot(c2_index_padded, roll_c2, color=col_dict[species], linewidth=1)
		ax2.set_ylim([-1,1])
		ax2.set_title('Chromosome 2')
		
		plt.savefig('Kmer_Location_Visualizations/'+species+'.png')

df.to_csv('Kmer_Locations.csv')
