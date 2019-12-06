import os 
import pandas as pd
import re
sp_dict = {}
out_list = []
species = []
j=0
n=5

rmash = pd.read_csv("rmash_output.csv")
strains= ['Brucella abortus','Brucella canis','Brucella ceti', 'Brucella inopinata', 'Brucella melitensis', 'Brucella microti', 'Brucella neotomae', 'Brucella ovis', 'Brucella pinnipedialis', 'Brucella sp', 'Brucella suis', 'Ochrobactrum anthropi', 'Ochrobactrum intermedium']

while j< len(rmash):
	sample = rmash['sample'].iat[j]
	group = rmash[j:j+n]
	k, l, most, distance =0, 0, 0, 0 #initializing the counters to zero
	
	for key in strains:
		sp_dict[key] = 0
	
	#incrementing the correct counter for every occurance of a species name within a group 
	while k< len(group):
		tax_sp = group['taxonomic_species'].iat[k]
		match= re.search('Brucella sp.',tax_sp)
		if match != None: 
			tax_sp_ver = 'Brucella sp'
		if match == None:
			tax_sp_ver = tax_sp
		sp_dict[tax_sp_ver]+=1
		k+=1

	#retreving the strain name with the highest counter value
	for key in sp_dict:
		if sp_dict[key] > most:
			most = sp_dict[key]
			strain = key

	#retreving the distance from the sample to the farthest match with of the correct species
	while l < len(group):
		if group['taxonomic_species'].iat[l] == strain:
			dist = group['distance'].iat[l]
			if dist > distance:
				distance = dist
		l+=1
	out_in = [sample, strain, distance]
	out_list.append(out_in)
	j+=n

df = pd.DataFrame.from_records(out_list, columns = ['Sample','Rmash Species', 'Max Distance'])
df.to_csv('rmash_files/rmash_filtered_output.csv', index = False)

