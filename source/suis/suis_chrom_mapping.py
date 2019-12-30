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
ranks = pd.read_csv('suis_biovar_files/Suis_Ranks.csv', index_col = 0)
suis_biovar = pd.read_csv("suis_biovar_files/Suis_Biovar.csv", index_col = 0)
suis_biovar = suis_biovar[suis_biovar['Biovar Classification'] != 'No Matched Biovar']
suis_files = suis_biovar.index.values

suis_metadata = metadata.loc[metadata['Sample'].isin(suis_files)]
suis_metadata = suis_metadata[suis_metadata['Number of Contigs'] == 2] # only complete sequences (a.k.a 2 contigs)

suis_files_contig_filtered = list(suis_metadata['Sample'])
for i, row in suis_biovar.iterrows():
	if row.name not in suis_files_contig_filtered:
		meta_ref = suis_biovar.drop(index = row.name)
meta_ref = meta_ref.drop_duplicates(subset = 'Biovar Classification', keep='first') # one of each species
meta_ref.to_csv('suis_biovar_files/Suis_Refrence_Genomes.csv')

biovars = list(meta_ref['Biovar Classification'])
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
alt_ref_strain = {}
#print(meta_ref)
c1_dict, c2_dict = {}, {}
for kmer in ranks.index: 
	c1_dict[kmer] = np.inf
	c2_dict[kmer] = np.inf

#Choosing and locating a refrence genome 
def refrence_location(biovar):
	ref_row = meta_ref[meta_ref['Biovar Classification']==biovar]
	ref_file = ref_row.index.values
	ref_file_location = 'Approved_Sequences/'+ref_file[0]+'.fna'
	return(ref_file_location, ref_file[0])

def alt_ref(biovar):
	sort_rank = ranks.sort_values(by = [biovar+" RANK"]).head(1)
	kmer = list(sort_rank.index.values)
	sort_rank = sort_rank.T.sort_values(by = kmer, ascending=False).head(1)
	alt_st = sort_rank.index.values[0][0:len(sort_rank.index.values[0])-6]
	alt_ref_file_location, ref_file = refrence_location(alt_st)
	return(alt_ref_file_location, ref_file)

#Creates 2 haystacks from the refrence genome, one for each chromosome
def haystacks(biovar, soa):
	if soa == 'self':
		location, ref_file= refrence_location(biovar)
	if soa == 'alt':
		location, ref_file = alt_ref(biovar)
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
def needles(biovar):
	kmers = ranks.index.values
	rank = ranks[biovar+' RANK']
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
	soa, biovar = tup
	aho = needles(biovar)
	chromosome1, chromosome2, ref_file = haystacks(biovar, soa)
	c1_index, c1_rank, c2_index, c2_rank = [], [], [], []
	for end_index, (kmer, rank) in aho.iter(chromosome1):
		c1_index.append(end_index)
		c1_dict[kmer] = end_index
		c1_rank.append(rank)
	for end_index, (kmer, rank) in aho.iter(chromosome2):
		c2_index.append(end_index)
		c2_dict[kmer] = end_index
		c2_rank.append(rank)
	return(soa, biovar, ref_file,  c1_index, c1_rank, c2_index, c2_rank, c1_dict, c2_dict)

col_list = [x+y+z for y in biovars for x in ['C1 ', 'C2 '] for z in [' self', ' alt']]
df = pd.DataFrame(index = ranks.index, columns = col_list)
in_list = [(x,y) for y in biovars for x in ['self', 'alt']]


with ProcessPoolExecutor(len(biovars)) as ppe:
	for soa, biovar, ref_file, c1_index, c1_rank, c2_index, c2_rank, c1_dict, c2_dict in ppe.map(find_needles, in_list):
		#For any kmer that is not found in the refrence genome it is assigned a value of zero
		
		fig, (ax1, ax2) = plt.subplots(1,2, sharey = True, figsize=(17,7))
		fig.suptitle(biovar)
		ax1.set_ylim([-1,1])
		ax1.set_title('Chromosome 1')
		ax2.set_ylim([-1,1])
		ax2.set_title('Chromosome 2')

		if len(c1_index) > 0:
			c1_index_padded = range(0, max(c1_index))
			c1_rank_padded = np.zeros(max(c1_index))

			for i in range(0,len(c1_rank)):
				index = c1_index[i]
				c1_rank_padded[index-1]=c1_rank[i]
			
			df['C1 '+biovar+' '+soa] = c1_dict.values()
			df_c1_rank = pd.DataFrame(c1_rank_padded)
			roll_c1 = df_c1_rank.rolling(window=1000).mean()

			ax1.plot(c1_index_padded, c1_rank_padded, color= '#dedcdf', linewidth=0.25, ls = ':')
			ax1.plot(c1_index_padded, roll_c1, linewidth=1)

		if len(c2_index) > 0:	
			c2_index_padded = range(0, max(c2_index))
			c2_rank_padded = np.zeros(max(c2_index))

			for i in range(0,len(c2_rank)):
				index = c2_index[i]
				c2_rank_padded[index-1]=c2_rank[i]
		
			df['C2 '+biovar+' '+soa] = c2_dict.values()
			df_c2_rank = pd.DataFrame(c2_rank_padded)
			roll_c2 = df_c2_rank.rolling(window=1000).mean()
		
			ax2.plot(c2_index_padded, c2_rank_padded, color='#dedcdf', linewidth=0.25, ls = ':')
			ax2.plot(c2_index_padded, roll_c2, linewidth=1)

		name =biovar+'_'+soa+'.png'.replace(' ', '_')
		plt.savefig('Visualizations/Kmer_Locations/'+name)

df.to_csv('suis_biovar_files/Suis_Kmer_Locations.csv')
