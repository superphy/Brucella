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
ranks = ranks.filter(axis = 1, like='RANK')
strain_occurances = pd.read_csv("Species_Occurrence.csv", index_col = 0)
strain = strain_occurances.index.values 
col_dict = {'Brucella sp':'#ebc04b', 'Brucella suis':'#f8a4d0','Brucella abortus': '#72a24a','Brucella canis': '#b8642f','Brucella ceti':'#59cf9e', 'Brucella melitensis':'#e8382e', 'Brucella neotomae':'#6c86e2', 'Brucella ovis':'#ba6eb4', 'Brucella pinnipedialis':'#f4792a' }
meta_ref = metadata[metadata['Number of Contigs'] == 2] # only complete sequences (a.k.a 2 contigs)
meta_ref = meta_ref.drop_duplicates(subset = 'Strain', keep='first') # one of each strain
meta_ref.to_csv('Refrence_Genomes.csv')
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
alt_ref_strain = {}

c1_dict, c2_dict = {}, {}
for kmer in ranks.index: 
	c1_dict[kmer] = np.inf
	c2_dict[kmer] = np.inf

#Choosing and locating a refrence genome 
def refrence_location(strain):
	ref_row = meta_ref[meta_ref['Strain']==strain]
	ref_file = list(ref_row['Sample'])
	ref_file_location = 'Approved_Sequences/'+ref_file[0]+'.fna'
	return(ref_file_location, ref_file[0])

def alt_ref(strain):
	sort_rank = ranks.sort_values(by = [strain+" RANK"]).head(1)
	kmer = list(sort_rank.index.values)
	sort_rank = sort_rank.T.sort_values(by = kmer, ascending=False).head(1)
	alt_st = sort_rank.index.values[0][0:len(sort_rank.index.values[0])-5]
	alt_ref_strain[strain]=alt_st
	alt_ref_file_location, ref_file = refrence_location(alt_st)
	return(alt_ref_file_location, ref_file)

#Creates 2 haystacks from the refrence genome, one for each chromosome
def haystacks(strain, soa):
	if soa == 'self':
		location, ref_file= refrence_location(strain)
	if soa == 'alt':
		location, ref_file = alt_ref(strain)
	sequence = []
	with open(location) as handle:
		for title, seq in SimpleFastaParser(handle):
			sequence.append(seq)
	chromosome1 = sequence[0]
	chromosome2 = sequence[1]
	return(chromosome1, chromosome2, ref_file)

#Generates the reverse complement of a given kmer
def reverse_complement(kmer):
	return ''.join([complement[base] for base in kmer[::-1]])

#Creates the needles (all kmers)
def needles(strain):
	kmers = ranks.index.values
	rank = ranks[strain+' RANK']
	aho = ahocorasick.Automaton()
	for i in range(0,len(kmers)):
		kmer = kmers[i]
		rev_c = reverse_complement(kmer)
		aho.add_word(kmer,(kmer,rank[kmer]))
		aho.add_word(rev_c,(kmer,rank[kmer]))
	aho.make_automaton()
	return(aho)

#Runs the ahocorasick algorithm to find the location of the kmers (needles) in the refrence genome (haystack) 
def find_needles(tup):
	soa, strain = tup
	aho = needles(strain)
	chromosome1, chromosome2, ref_file = haystacks(strain, soa)
	c1_index, c1_rank, c2_index, c2_rank = [], [], [], []
	for end_index, (kmer, rank) in aho.iter(chromosome1):
		c1_index.append(end_index)
		c1_dict[kmer] = end_index
		c1_rank.append(rank)
	for end_index, (kmer, rank) in aho.iter(chromosome2):
		c2_index.append(end_index)
		c2_dict[kmer] = end_index
		c2_rank.append(rank)
	return(soa, strain, ref_file, c1_index, c1_rank, c2_index, c2_rank, c1_dict, c2_dict)

df = pd.DataFrame(index = ranks.index)
in_list = [(x,y) for y in strain for x in ['self', 'alt']]

with ProcessPoolExecutor(len(strain)) as ppe:
	for soa, strain, ref_file, c1_index, c1_rank, c2_index, c2_rank, c1_dict, c2_dict in ppe.map(find_needles, in_list):
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
		df['C1 '+strain+' '+soa] = c1_dict.values()
		df['C2 '+strain+' '+soa] = c2_dict.values()

		#Generating the moving average of the rank values
		df_c1_rank = pd.DataFrame(c1_rank_padded)
		roll_c1 = df_c1_rank.rolling(window=1000).mean()
		df_c2_rank = pd.DataFrame(c2_rank_padded)
		roll_c2 = df_c2_rank.rolling(window=1000).mean()

		#Generating plots
		fig, (ax1, ax2) = plt.subplots(1,2, sharey = True, figsize=(17,7))
		fig.suptitle(strain+" "+soa+'\n'+ref_file)
			
		ax1.plot(c1_index_padded, c1_rank_padded, color= '#dedcdf', linewidth=0.25, ls = ':')
		ax1.plot(c1_index_padded, roll_c1, color=col_dict[strain], linewidth=1)
		ax1.set_ylim([-1,1])
		ax1.set_title('Chromosome 1')
			
		ax2.plot(c2_index_padded, c2_rank_padded, color='#dedcdf', linewidth=0.25, ls = ':')
		ax2.plot(c2_index_padded, roll_c2, color=col_dict[strain], linewidth=1)
		ax2.set_ylim([-1,1])
		ax2.set_title('Chromosome 2')
		
		name =strain+'_'+soa+'.png'.replace(' ', '_')
		plt.savefig('Visualizations/Kmer_Locations/'+name)
				
df.to_csv('Kmer_Locations.csv')
