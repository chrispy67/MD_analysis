from shutil import copyfile
import numpy as np

"""Written by Marlo Zorman (2022)
----------------------------------
This script fixes the filenames of colvar files required by mbar.py reweighting scheme.
Filenames include the center of the restraint for each umbrella window
"""

cvs = np.linspace(4.95, 7.2, num=64)

for i in range(len(cvs)):
	try:
		com = cvs[i]
		copyfile("colvar_multi_"+str(i), "rcwin_"+str(com)+"_0_us.dat")
		#print(cut_cvs[i])
	except:
		print("cant find file: "+str(i))
