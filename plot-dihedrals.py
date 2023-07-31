import matplotlib.pyplot as plt
import sys
import pickle as pickle

pwd = '{}/'.format(sys.path[0])

files = ['anti/1-3130.pickle', 'anti-stable/1-3130.pickle', 'para/1-3130.pickle']
labels = ['anti-parallel 1', 'anti-parallel 2', 'parallel']

for file, title in zip (files, labels):
    #retrieve pickled objects
    with open(file, 'rb') as f:
        hot, cold = pickle.load(f)
#'hot' turns
    plt.figure(figsize=(10, 8))
    plt.hist2d(hot[:,:,0].flatten(), hot[:,:,1].flatten(),
    bins=[75,75], cmap='cividis', alpha=0.6, label='hot')

#plot settings
    plt.title('{} Hot Turns'.format(title), fontsize=20)
    plt.xlabel(r"$\phi$", weight="bold", fontsize=12)
    plt.ylabel(r"$\psi$", weight="bold", fontsize=12)
    plt.show()
#'cold' turns
    plt.figure(figsize=(10, 8)) 
    plt.hist2d(cold[0:,:,0].flatten(), cold[:,:,1].flatten(),
    bins=[75,75], cmap='cividis', alpha=0.6, label='cold')

#plot settings
    plt.title('{} Cold Turns'.format(title), fontsize=20)
    plt.xlabel(r"$\phi$", weight="bold", fontsize=12)
    plt.ylabel(r"$\psi$", weight="bold", fontsize=12)
    plt.show()
