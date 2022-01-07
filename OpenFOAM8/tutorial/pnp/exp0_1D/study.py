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

NL = 5
NZ = 10

LCvalues = np.linspace(1, 6,num=NL)
Zvalues = np.linspace(1.29, 1.31,num=NZ)
MAX = np.zeros((NL,NZ))

for i,Lval in enumerate(LCvalues):

	# print(caseDir, val)
	modifLine(transportProperties, 19, "LC           LC       [ 0 0 -1 0 0 0 0 ] {};".format(Lval))

	for j,Zval in enumerate(Zvalues):

		modifLine(transportProperties, 24, "\tZ        Z [0 0 -1 0 0 0 0]     {};".format(Zval))

		os.system("pnpSingleSpeciesPicard1")

		os.system("postProcess -func sampleDict -latestTime")

		with open('./postProcessing/sampleDict/1/line_C_V.xy','r') as f:
			lines = f.readlines()

		N = len(lines)
		x = np.zeros(N, dtype=np.float64)
		C = np.zeros(N, dtype=np.float64)
		V = np.zeros(N, dtype=np.float64)

		for k,line in enumerate(lines):
			tmp = list(map(np.float64, line.split()))
			x[k] = tmp[0]
			C[k] = tmp[1]
			V[k] = tmp[2]
		
		MAX[i,j] = np.max(C)

LL, ZZ = np.meshgrid(LCvalues, Zvalues)

print(LL.shape, ZZ.shape, MAX.shape)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(LL, ZZ, MAX.T)
ax.set_xlabel('LC')
ax.set_ylabel('Z')
ax.set_zlabel('max(C)')

plt.show()
