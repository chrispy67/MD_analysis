import matplotlib.pyplot as plt
import glob
import pandas as pd
import sys
print(sys.path[0])

pwd = '{}/'.format(sys.path[0])
data_files = [name.replace(pwd, '') for name in sorted(glob.glob(pwd + '*.dat'))] #loads all data files
#sorted() sorts everything by alphabetical order!


plt.figure(figsize=(9, 6))
turns = [14, 30, 46, 62, 78, 94, 110, 126, 142, 158, 174]

colors = ['red', 'cornflowerblue', 'orangered']

#turn labels
hot = [[30,31], [62,63], [94,95], [126,127], [158,159]]
cold = [[14,15], [46,47], [78,79], [110,111], [142,143]]

for i in range (len(data_files)):
    file = data_files[i]
    print(file)
    df = pd.read_csv(file, delimiter="\t", header=0, names=['resid', 'rmf'])
    plt.plot(df.resid, df.rmf, label=str(file), alpha=0.8, color=colors[i], linewidth=1.75)
    plt.xlabel('Residue Index', fontsize=18)
    plt.ylabel('RMF (Å)', fontsize=18)


###Plots vertical lines that indicate key turn residues
for i in range(len(hot)):
    arr1 = hot[i]
    plt.axvspan(arr1[0], arr1[-1], zorder=0, alpha=0.2, color='orange')
    arr2 = cold[i]
    plt.axvspan(arr2[0], arr2[-1], zorder=0, alpha=0.2, color='orange')


plt.xticks(turns, fontsize=11)
plt.yticks(fontsize=11)

plt.legend(bbox_to_anchor=(1, 1), loc='upper right')
#plt.tight_layout()
plt.title('RMF Analysis of Bound Fibers (300ns)', fontsize=20)


plt.show()


##plotting derivatives

        #blue_x, blue_y = blue.array[0], blue.array[i]
        #blue_dy = np.diff(blue_y / np.diff(blue_x))
        #blue_dx = (np.array(blue_x)[:-1] + np.array(blue_x)[1:]) / 2
        #plt.plot(blue_dx, blue_dy)