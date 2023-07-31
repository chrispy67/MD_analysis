import matplotlib.pyplot as plt
from sys import argv
from math import pi
import numpy as np
from os.path import exists
import shlex

## -----------------------------------------------------------------------------
## Sript to plot histogram of CVs (written by Marlo Zorman (2022))
## -----------------------------------------------------------------------------
n = 73
filenames = ["../colvar/prod/colvar_multi_"+str(i) for i in range(n)]

#for placing vertical lines to quantify drift of extra windows
og_sep = np.linspace(0, 2.25, num=64)
window_num = np.asarray([64, 65, 66, 67, 68, 69, 70, 71, 72])
spliced_window = np.asarray([3, 4, 10, 13, 0, 1, 6, 7, 8])
delta = og_sep[1]-og_sep[0]


# fig = plt.figure()
for i, filename in enumerate(filenames):
    if exists(filename):
        file = open(filename, "r")
        data = []
        for line in file:
            if "FIELDS" not in line:
                data.append(float(line.split(" ")[2]))
        print(filename, len(data))
        file.close()
        plt.hist(data, bins=50, alpha=0.75)
        if i >= 64: #plot extra added windows thicker and easier to see. only works because extra windows => 64 and are in order.
            plt.hist(data, bins=50, alpha=0.9, edgecolor='black', linewidth=0.15, label=str(i))
            plt.legend()
            print(filename, len(data))

###Plotting center of new windows to achieve convergence in reweighting.
for i in range(len(spliced_window)): #plotting where new windows are going to be centered
    #(window[i]*regular interval) + half-delta from add_window_driver + min separation
    if spliced_window[i] in [0, 6]:
        x_val = spliced_window[i]*delta + (delta/2) + 4.95
    if spliced_window[i] == 1:
        x_val = spliced_window[i]* delta +(delta*0.75) + 4.95
    if spliced_window[i] == 3:
        x_val = spliced_window[i]*delta +(delta/2) + 4.95
    if spliced_window[i] in [4, 7, 8]:
        x_val = spliced_window[i]*delta + (delta/3) + 4.95
    if spliced_window[i] in [10, 13]:
        x_val = spliced_window[i]*delta + (delta*(7/8)) + 4.95
    
    
    plt.axvline(x=x_val, color='black', alpha=0.8, linewidth=0.5, linestyle='--')

###Plotting center of original windows with standard separation
for i in np.linspace(4.95, 7.2, num=64):
    plt.axvspan(i, i, zorder=0, alpha=0.4, color='orange', label='turns')

    # plt.plot(x_axis, y_axis)
#plt.xlim([0,1.9])
plt.title('Umbrella Sampling Histogram (Parallel)', fontsize=22)
plt.xlabel("Center of Mass Separation (nm)", fontsize=16)
plt.ylabel("Count", fontsize=16)

plt.xticks(fontsize=14) #blowing up font size because this is a small figure


plt.show()
