import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import sys
import glob

pwd = '{}/'.format(sys.path[0])

#####CREATE AN ARRAY OF DATAFRAMES READ FROM A DATAFILE!!! Yes it really is this easy
metrics = ['d.x', 'protein_lat.x', 'protein_lat.y', 'surface_lat.x', 'surface_lat.y']


###COLVAR file -> dataframe with indexes 
#In order for pandas to recognize headers as keys, delimiter must be uniform and perfect

dfs = [] #anti-stable-para 
data_files = [name.replace(pwd, '') for name in sorted(glob.glob(pwd + 'COLVAR-*'))]
for file in data_files:
    print(file)
    dfs.append(pd.read_csv(file, header=0, delimiter=' ', index_col=False))


for metric in metrics:
    steps = np.linspace(0, 15000, 15000)
    
    if metric != 'd.x':
        corrected_dfs = []
        #generating correction factor because position of frozen atom on surface is arbitrary
        for df in dfs:
            initial_value = df[metric][0]
            new_array = []
            #within each entry of array, I want to subtract the original difference
            for i in range(len(df)):
                corrected_values = df[metric][i] - initial_value
                new_array.append(float(corrected_values))
            corrected_dfs.append(new_array) #create new array of arrays the same way I did above
        #plot params
        plt.figure(figsize=(9,6))
        plt.scatter(steps, corrected_dfs[0], label='anti-parallel_1', s=0.5, color='red')
        plt.scatter(steps, corrected_dfs[1], label='para', s=0.5, color='cornflowerblue')
        plt.scatter(steps, corrected_dfs[2], label='ant_parallel_2', s=0.5, color='orange')

        #testing out derivatives
        dydx_0 = np.diff(corrected_dfs[0]) / np.diff(steps)
        dydx_1 = np.diff(corrected_dfs[1]) / np.diff(steps)
        dydx_2 = np.diff(corrected_dfs[2]) /  np.diff(steps)
        x = np.linspace(0, len(dydx_0), len(dydx_0))
        #plt.plot(x[0:-1:100], dydx_0[0:-1:100])
        plt.xlabel('steps (ts = 2 fs)', fontsize=15)
        plt.ylabel('ΔΔ COM (nm)', fontsize=15)
        plt.title(metric, fontsize=20)
        
        lgnd = plt.legend(loc='lower right',)
        plt.show()
