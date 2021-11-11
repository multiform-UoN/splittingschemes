#!/usr/bin/env python3
from pathlib import Path as path
import os
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

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


parameters = np.linspace(0, 1,num=5)
cases = {}

for val in parameters:

	caseDir = workingDir

	transportProperties = caseDir + "/constant/electrokineticProperties"

	print(caseDir, val)

	# modifLine(transportProperties, 24, "\tZ        Z [0 0 -1 0 0 0 0]     {};".format(val))
	modifLine(transportProperties, 19, "LC           LC       [ 0 0 -1 0 0 0 0 ] {};".format(val))

	# os.system("pnpSplittingSchemeFoam > log.pnpSplittingSchemeFoam_Z{}".format(val))
	os.system("pnpSplittingSchemeFoam")

	os.system("postProcess -func sampleDict -latestTime")

	with open('./postProcessing/sampleDict/1/line_C_V.xy','r') as f:
		lines = f.readlines()

	N = len(lines)
	x = np.zeros(N, dtype=np.float64)
	C = np.zeros(N, dtype=np.float64)
	V = np.zeros(N, dtype=np.float64)

	for i,l in enumerate(lines):
		tmp = list(map(np.float64, l.split()))
		x[i] = tmp[0]
		C[i] = tmp[1]
		V[i] = tmp[2]
	
	cases.update({"L={}".format(val):[x, C, V]})

for key in cases.keys():
	x = cases[key][0]
	C = cases[key][1]
	V = cases[key][2]
	plt.plot(x, C, label=f"concentration {key}")
	# plt.plot(x, V, label="potential")
	plt.xlabel(r"arc length $[m]$")
plt.legend()
plt.grid()
plt.show()
