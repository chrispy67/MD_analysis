import matplotlib.pyplot as plt
import sys
import gromacs.formats as gmx
import numpy as np  


columns = ['time','AA-interactions', 'protAA-hopg', 'protAA-surfG__G' 'protAA-SOL',
'surfAA-hopg', 'surfAA-protG__G', 'surfAA-SOL', ]
titles = ['N/A', 'turn-turn interactions', 'AA-hopg', 'AA-G__G', 'AA-SOL', 
'AA-hopg', 'AA-G__G', 'AA-SOL']


###Loading in data
blue = gmx.XVG(filename='energy-blue.xvg', names=columns)
red = gmx.XVG(filename='energy-red.xvg', names=columns)
para = gmx.XVG(filename='energy-para.xvg', names=columns)
#not sure what this gromacs library does, but it does offer a great way to load in 
#.xvg files that are otherwise difficult to graph outside of gnuplot. 

#print(blue.names[1])
#print(blue.array[1])access one indexed array inside the gmx.XVG() function.


for i in range(len(columns)):
    if i != 0:
        plt.figure(figsize=(10,8))
        plt.plot(blue.array[0], blue.array[i], label='anti parallel_1 {}'.format(columns[i]), color='red', alpha=0.8, linewidth=0.75)
        plt.plot(red.array[0], red.array[i], label='anti parallel_2 {}'.format(columns[i]), color='orange', alpha=0.8, linewidth=0.75)
        plt.plot(para.array[0], para.array[i], label='parallel {}'.format(columns[i]), color='cornflowerblue', alpha=0.8, linewidth=0.75)
        plt.title(titles[i], fontsize=22)

        #approx average values
        blue_avg = np.average(blue.array[i][1000:-1])
        red_avg = np.average(red.array[i][1000:-1])
        para_avg = np.average(para.array[i][1000:-1])

        #plotting hlines to approximate equil values
#        if i == 1:
#            plt.hlines(blue_avg, 0, 300000, color='maroon', linestyles='dashed', lw=0.5)
#            plt.text(1, blue_avg+0.2, '{:.1f}'.format(blue_avg), color='black') #labelling right on the line by hand, not too bad in this case

#            plt.hlines(red_avg, 0, 300000, color='darkorange', linestyles='dashed', lw=0.5)
#            plt.text(1, red_avg+0.2, '{:.1f}'.format(red_avg), color='black')

#            plt.hlines(para_avg, 0, 300000, color='darkblue', linestyles='dashed', lw=0.5)
#            plt.text(1, para_avg+0.2, '{:.1f}'.format(para_avg), color='black')

        
        plt.ylabel('SR LJ (kJ/mol)', fontsize=16)
        plt.xlabel('time (ps)', fontsize=16)
        plt.xlim(0, 300000)
        plt.legend()
        plt.show()
    else:
        continue

plt.show()




