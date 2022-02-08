#!/usr/bin/env python3
from pathlib import Path as path
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
# import plotly.express as px

def modifLine(path, riga, modifica):
	"""
    modifica una linea
    """
	infile = open(path, "r")
	contents = infile.readlines()
	infile.close()
	contents[riga] = modifica+"\n"
	outfile = open(path, 'w')
	for item in contents:
		outfile.write("%s" % item)
	outfile.close()
	pass

workingDir = str(path(os.path.dirname(os.path.abspath(__file__)))) #abspath gestisce il link simbolico

print("\n"+(150*"*")+"\nWORKING DIR : {0}\n".format(workingDir))
caseDir = workingDir
transportProperties = caseDir + "/constant/electrokineticProperties"

z1_num = 4
z1_vals = np.linspace(0.1, 0.5, num=z1_num)

fig = plt.figure()
ax = fig.add_subplot(111)

for i,val in enumerate(z1_vals):
	modifLine(transportProperties, 23, "    z1      z1      [0 -2  0 0 0 0 0]    {};".format(val))
	os.system("pnpSingleSpeciesPicard1")
	os.system("postProcess -func sampleDict -latestTime")
	with open('./postProcessing/sampleDict/2/line_C_V.xy','r') as file:
		lines = file.readlines()
		N = len(lines)
		x = np.zeros(N, dtype=np.float64)
		C = np.zeros(N, dtype=np.float64)
		V = np.zeros(N, dtype=np.float64)
		for k,line in enumerate(lines):
			tmp = list(map(np.float64, line.split()))
			x[k] = tmp[0]
			C[k] = tmp[1]
			V[k] = tmp[2]
		ax.plot(x, C)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$C$')
plt.grid()
plt.show()
