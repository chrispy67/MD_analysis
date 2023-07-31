import csv
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
import csv
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='calculating root-mean square fluctuations of residues')
parser.add_argument('-pdb', help='pdb file with correct index', type=str, default='prod.pdb')
parser.add_argument('-xtc', help='xtc file from trajectory', type=str, default='prod.xtc')
parser.add_argument('-index', help='atom index to calculate rmf (cannot have index that reaches across the pbc)', type=str, default='1-1565', nargs='+')
args = parser.parse_args()

indexes = np.asarray(args.index)
pdb = args.pdb
xtc = args.xtc


for i in range(len(indexes)):    
    i = indexes[i]
    print(i)
    u = mda.Universe(pdb, xtc)
    avg = align.AverageStructure(u, u, select="protein and name CA", ref_frame=0).run()
    ref = avg.universe
    aligner = align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    c_alphas = u.select_atoms('protein and name CA and bynum {}'.format(i))
    R = rms.RMSF(c_alphas).run() 
    plt.plot(c_alphas.resids, R.rmsf, label='{} (300ns)'.format(i))

    
    plt.xlabel('Residue Number')
    plt.ylabel('RMSF ($\AA$)')

    zip(c_alphas.resids, R.rmsf)

    with open('rmf{}.dat'.format(str(i)), 'a') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(c_alphas.resids, R.rmsf))

    plt.legend()


plt.title('Root Mean Fluctuations of 11c silk fibroin')
turns = [[14,15], [30,31], [46,47], [62,63], [78, 79], [94, 95], [110,111], [126, 127], [142,143], [158, 159], [174, 175]]

hot = [[30,31], [62,63], [94,95], [126,127], [158,159]]
cold = [[14,15], [46,47], [78,79], [110,111], [142,143]]

for i in range(len(hot)):
    arr1 = hot[i]
    plt.axvspan(arr1[0], arr1[-1], zorder=0, alpha=0.2, color='red', label='high rmf')
    arr2 = cold[i]
    plt.axvspan(arr2[0], arr2[-1], zorder=0, alpha=0.2, color='blue', label='low rmf')
plt.show()
#zip(c_alphas.resids, R.rmsf)


