import matplotlib.pyplot as plt
import sys
import pickle as pickle
import glob
import seaborn as sns

"""This script compares the water contacts of silk fibroin turns generated via
the bin_waters.py script. Essentially a Radial Distribution Function"""

pwd = '{}/'.format(sys.path[0])

###Load in pickled files, labels, color palette, etc.
pickles = [name.replace(pwd, '') for name in sorted(glob.glob(pwd + '*.pickle'))] #sorted() defaults to alphabetical order
colors = ['red', 'cornflowerblue', 'orange']
labels = ['anti 1', 'para', 'anti-2']

fig = plt.figure(figsize=(11,7))
i = 0
dfs = []

###Retrieve the pickles (open the jar):
for file in pickles: #anti -> para -> stable
    color = colors[i]

    with open(file, 'rb') as f:
        dist_array = pickle.load(f)
    
    #plot histogram of water contacts
    #plt.hist(dist_array, bins=60, range=(0, 15), density=True, color=color, alpha=0.2, 
    #edgecolor=color, linewidth=0.5, label=str(file)) #each bin=0.25Å
    sns.kdeplot(dist_array, color=color, common_norm=False, label=labels[i]) #KDE
    print(file)
    i += 1


###Plot Decoration
plt.legend()
plt.title('KDE of Water Contacts', fontsize=20)
plt.xlabel('Distance (Å)', fontsize=15)
plt.ylabel('normalized count', fontsize=15)
plt.xlim(0, 15)

plt.show()
