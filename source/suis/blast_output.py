import skbio.io
import pandas as pd
import numpy as np 
pd.options.mode.use_inf_as_na = True

blast_output_path = 'suis_biovar_files/blast/blast_search_output.tsv'
with open(blast_output_path) as fh:
	blast_df = skbio.io.read(fh, format="blast+6",into=pd.DataFrame,default_columns=True) # turn the blast output into a pandas df 

blast_df = blast_df[blast_df['pident']==100]
sseqids = list(blast_df['sseqid'].unique())
blast_df.to_csv('suis_biovar_files/Blast_Output.csv')

def ampliconf(start_a, end_a, start_b, end_b):
	primer_min = min(start_a, start_b) # determine if the start val of sequence A or sequence B is earlier in the contig
	if primer_min == start_a: # if A is first
		amp_start = min(start_a, end_a) # start is the smallest of the start and end values of seq A
		amp_end = max(start_b, end_b) # end is largest of of the start and end values of seq B
	if primer_min == start_b: # if B is first
		amp_start = min(start_b, end_b) # start is the smallest of the start and end values of seq B
		amp_end = max(start_a, end_a) # end is largest of of the start and end values of seq A
	return(amp_end-amp_start) # amplicon size = amp_end - amp_start

def primer_one(df, sseqid):
	BMEI1426 = df[df['qstart']==1.0] # df with hits that match the BMEI1426 sequence
	BMEI1427 = df[df['qstart']==51.0] # df with hits that match the BMEI1427 sequence
	if len(BMEI1426)>=1 and len(BMEI1427)>=1: # as long as the df is not empty
		start_26 = int(BMEI1426['sstart']) # start value of the BMEI1426 sequence
		end_26 = int(BMEI1426['send']) # end value of the BMEI1427 sequence
		start_27 = int(BMEI1427['sstart']) # start value of the BMEI1426 sequence
		end_27 = int(BMEI1427['send']) # end value of the BMEI1427 sequence
		amplicon = ampliconf(start_26, end_26, start_27, end_27)
		return(amplicon)
	else:
		return(np.inf) # if either of the df are empty, there is no match for this primer so an empty cell is returned

def primer_two(df, sseqid):
	#multicases!!
	BR1080f = df[df['qstart']==1.0] # df with hits that match the BR1080f sequence
	BR1080r = df[df['qstart']==51.0] # df with hits that match the BR1080r sequence

	if len(BR1080f)==1 and len(BR1080r)==1: # if there is only one set of primers found, behave like the other primer functions
		start_f = int(BR1080f['sstart'])
		end_f = int(BR1080f['send'])
		start_r = int(BR1080r['sstart'])
		end_r = int(BR1080r['send'])
		amplicon = ampliconf(start_f, end_f, start_r, end_r)
		return(amplicon)
	
	if len(BR1080f)==2 and len(BR1080r)==2: # if there are two sets of primers:
		start_f_1 = int(BR1080f.iat[0,8])
		end_f_1 = int(BR1080f.iat[0,9])
		start_r_1 = int(BR1080r.iat[0,8])
		end_r_1 = int(BR1080r.iat[0,9])
		amplicon_1 = ampliconf(start_f_1, end_f_1, start_r_1, end_r_1)
		
		start_f_2 = int(BR1080f.iat[1,8])
		end_f_2 = int(BR1080f.iat[1,9])
		start_r_2 = int(BR1080r.iat[1,8])
		end_r_2 = int(BR1080r.iat[1,9])
		amplicon_2 = ampliconf(start_f_2, end_f_2, start_r_2, end_r_2)
		if amplicon_1 == amplicon_2:
			return(amplicon_1)
		else:	
			return([amplicon_1, amplicon_2]) # return a list of the amplicons of each primer		
	else:
		return(np.inf)

def primer_three(df, sseqid):
	BMEI1688 = df[df['qstart']==1.0] # df with hits that match the BMEI1688 sequence
	BMEI1687 = df[df['qstart']==51.0] # df with hits that match the BMEI1687 sequence
	if len(BMEI1688)>=1 and len(BMEI1687)>=1: # as long as the df is not empty
		start_88 = int(BMEI1688['sstart']) # start value of the BMEI1688 sequence
		end_88 = int(BMEI1688['send']) # end value of the BMEI1688 sequence
		start_87 = int(BMEI1687['sstart']) # start value of the BMEI1687 sequence
		end_87 = int(BMEI1687['send']) # end value of the BMEI1687 sequence
		amplicon = ampliconf(start_88, end_88, start_87, end_87) 	
		return(amplicon)
	else:
		return(np.inf) # if either of the df are empty, there is no match for this primer so an empty cell is returned

def primer_four(df, sseqid):
	BMEI0205f = df[df['qstart']==1.0] # df with hits that match the BMEI0205f sequence
	BMEI0205r = df[df['qstart']==52.0] # df with hits that match the BMEI0205r sequence
	if len(BMEI0205f)>=1 and len(BMEI0205r)>=1: # as long as the df is not empty
		start_f = int(BMEI0205f['sstart']) # start value of the BMEI0205f sequence
		end_f = int(BMEI0205f['send']) # end value of the BMEI0205f sequence
		start_r = int(BMEI0205r['sstart']) # start value of the BMEI0205r sequence
		end_r = int(BMEI0205r['send']) # send value of the BMEI0205r sequence
		amplicon = ampliconf(start_f, end_f, start_r, end_r)
		return(amplicon)
	else:
		return(np.inf) # if either of the df are empty, there is no match for this primer so an empty cell is returned

def primer_search(sseqid):
	unique_df = blast_df[blast_df['sseqid']==sseqid] # df containing only hits from a single contig
	primer_one_df = unique_df[unique_df['qseqid']=='BMEI1426-BMEI1427'] # subset of unique_df with only primer 1 data
	primer_two_df = unique_df[unique_df['qseqid']=='BR1080f-BR1080r'] # subset of unique_df with only primer 2 data
	primer_three_df = unique_df[unique_df['qseqid']=='BMEI1688-BMEI1687'] # subset of unique_df with only primer 3 data
	primer_four_df = unique_df[unique_df['qseqid']=='BMEI0205f-BMEI0205r'] # subset of unique_df with only primer 4 data

	if len(primer_one_df) >1:
		p1_amp = primer_one(primer_one_df, sseqid)
	else:
		p1_amp = np.inf

	if len(primer_two_df)>1:
		p2_amp = primer_two(primer_two_df, sseqid)
	else:
		p2_amp = np.inf
	
	if len(primer_three_df) >1:
		p3_amp = primer_three(primer_three_df, sseqid)
	else:
		p3_amp = np.inf

	if len(primer_four_df)>1:
		p4_amp = primer_four(primer_four_df, sseqid)
	else:
		p4_amp = np.inf	

	return(sseqid, [p1_amp, p2_amp, p3_amp, p4_amp])

df = pd.DataFrame(index = ['Primer One (BMEI1426-BMEI1427)', 'Primer Two (BR1080f-BR1080r)', 'Primer Three (BMEI1688-BMEI1687)', 'Primer Four (BMEI0205f-BMEI0205r)'])
for seqid in sseqids:
	sseqid, amplicons = primer_search(seqid)
	df[sseqid] = amplicons

df = df.dropna(axis='columns', how='all') # drops any columns that did not return any primer matches
df = df.T
indexes = df.index.values # all of the sseqids that have at leat one primer match
files = [str(index[0:9]) for index in indexes] # the file id's, parsed from the sseqid name
contigs = [index[10:len(index)] for index in indexes] # the contigs, parsed from the sseqid name
df.insert(loc =0, column = 'File ID',value = files) 
df.insert(loc =1, column = 'Contig' , value = contigs)
df.to_csv('suis_biovar_files/Blast_Primer_Matching.csv', index = False)
