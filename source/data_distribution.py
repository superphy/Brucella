import pandas as pd 
import matplotlib.pyplot as plt

metadata = pd.read_csv('../Metadata.csv')
sp_occ = pd.read_csv('../Species_Occurrence.csv', index_col = 0)
col_dict = {'Brucella sp':'#ebc04b', 'Brucella suis':'#f8a4d0','Brucella abortus': '#72a24a','Brucella canis': '#b8642f','Brucella ceti':'#59cf9e', 'Brucella melitensis':'#e8382e', 'Brucella neotomae':'#6c86e2', 'Brucella ovis':'#ba6eb4', 'Brucella pinnipedialis':'#f4792a' }


species = sp_occ.index.values 
occurances = list(sp_occ['Species Occurrence'])

bar_plot = plt.bar(species, occurances)

for i, index in enumerate(species):
	bar_plot[i].set_color(col_dict[index])
	#bar_plot.text(occurances[i]+3, i+3, str(occurances[i]))
plt.show()