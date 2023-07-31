from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import argparse
import sys
from re import search
import pickle as pickle
pwd = '{}/'.format(sys.path[0])

###Parse args
parser = argparse.ArgumentParser(description="""calculating dihedral angles of turn residues. 
Will make subplots to compare if supplied extra arguments for -index""")

parser.add_argument('-gro', help='pdb file for MDAnalysis', type=str, default='prod.gro')
parser.add_argument('-xtc', help="""xtc file for MDAnalysis. 
Truncated trajectory should come LAST and is used to asses equilibrium values""", type=str, default='prod.xtc')

parser.add_argument('-index', help='atom index to calculate dihedral angles', type=str, default='1-3130 3131-6260', nargs='+')
args = parser.parse_args()

indexes = np.asarray(args.index)
gro = args.gro
xtc = args.xtc
fig, ax = plt.subplots(3, len(indexes), figsize=(18,15))


###This is a pretty bespoke script, but allows the Ramachandran plots of multiple different
###Experimental setups to be compared side by side
for i in range(len(indexes)):
    atom_index = indexes[i]
    u = mda.Universe(gro, xtc)
    
    ###Turns facing water is not universal. However, turn type is.

    ###Hot turn analysis 
    hot = [30,31, 62,63, 94,95, 126,127, 158,159] # universal
    selection_hot = 'bynum '+atom_index+' and('
    for turn in hot:
        selection_hot += "resid "+str(turn)
        if turn != 159:
            selection_hot +=" or "
    selection_hot += ")"
    print(selection_hot)

    r = u.select_atoms(selection_hot)
    R = Ramachandran(r).run()
    
    if i == 0:
        title_col0 = 'Cold Turns (water facing) {}'.format(atom_index)
    if i == 1:
        title_col0 = 'Cold Turns (protein facing) {}'.format(atom_index)
    
    ###Generating the plt.subplots() for main figure
    ax[i, 0].hist2d(R.angles[:,:,0].flatten(), R.angles[:,:,1].flatten(), bins=[75,75], cmap='Purples', alpha=0.8)
    ax[i, 0].set_title(title_col0, size="small", fontweight='bold')
    ax[i, 0].set_xlabel(r"$\phi$", weight="bold", size="small")
    ax[i, 0].set_ylabel(r"$\psi$", weight="bold", size="small")


    ###Cold turn analysis
    selection_cold = 'bynum '+atom_index+' and('
    cold = [14,15, 46,47, 78,79, 110,111, 142,143] #universal 
    for turn in cold:
        selection_cold += "resid "+str(turn)
        if turn != 143:
            selection_cold +=" or "
    selection_cold += ")"
    print(selection_cold)    
    
    r_cold = u.select_atoms(selection_cold)
    R_cold = Ramachandran(r_cold).run()

    if i == 0:
        title_col1 = 'Hot Turns (protein facing) {}'.format(atom_index)
    if i == 1:
        title_col1 = 'Hot Turns (water facing) {}'.format(atom_index)

    #generating the plt.subplots() for the main figure
    #histcold = ax[i, 1].hist2d(R_cold.angles[:,:,0].flatten(), R_cold.angles[:,:,1].flatten(), bins=[75,75], cmap='Purples')[0]
    ax[i, 1].hist2d(R_cold.angles[:,:,0].flatten(), R_cold.angles[:,:,1].flatten(), bins=[75,75], cmap='Purples')
    ax[i, 1].set_title(title_col1, size="small", fontweight='bold')
    ax[i, 1].set_xlabel(r"$\phi$", weight="bold", size="small")
    ax[i, 1].set_ylabel(r"$\psi$", weight="bold", size="small")


    ###Creating arrays to analyze water and protein interactions
    if search('1-3130', atom_index):
        coldWater_interaction = R_cold.angles
        hotProtein_interaction = R.angles
    if search('3131-6260', atom_index):
        coldProtein_interaction = R_cold.angles
        hotWater_interaction = R.angles
    
    ###Doing math with 3D arrays to study water-protein interactions

    #diff = (np.abs(histhot-histcold))
    #plt.hist2d(np.abs(histhot-histcold))
    #plt.hist2d(diff[:, 0], diff[:, 1].flatten(), bins=[75,75], cmap='viridis')
    #print(np.shape(diff))

    ###Export saved dihedral analysis with pickle. 
    with open('{}.pickle'.format(atom_index), 'wb') as f:
        pickle.dump([R.angles, R_cold.angles], f)

    #get the objects back:
    #with open('{}.pickle'.format(atom_index), 'rb') as f:
        #R.angles, R_cold.angles = pickle.load(f)


for a in ax.flat:
    a.tick_params(axis='both',labelsize=6)
plt.show()
plt.tight_layout()
#plt.savefig('dihedrals.png')