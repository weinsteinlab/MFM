import numpy as np
from gridData import OpenDX
from gridData import Grid
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib
matplotlib.use('Agg')

membrane_position = 40
#valency = -4
#area_per_lipid = 70
#everage_density = 0.1
valency = -1
area_per_lipid = 65
everage_density = 0.2

posilist = 'posi_list'
posilist = posilist.split(' ')

dx1 = OpenDX.field(0)
dx2 = OpenDX.field(0)

plotcount = len(posilist)
plt.figure(figsize=(10, 7.5))

for i in range(plotcount):
	dx1.read("./%s/cmbChargeDensity.dx" %posilist[i])
	histogram1, edges = dx1.histogramdd()
	xcor = (edges[0][1:]+edges[0][:-1])/2
	ycor = (edges[1][1:]+edges[1][:-1])/2
	zcor = (edges[2][1:]+edges[2][:-1])/2
	membrane_coord = np.where(zcor==membrane_position)[0][0]
	density_map = histogram1[:,:,membrane_coord]*area_per_lipid/valency
	Cpotential = np.log(density_map*(1-everage_density)/((1-density_map)*everage_density))
	dx2.read("./%s/Potential.dx" %posilist[i])
	histogram2, edges = dx2.histogramdd()
	Epotential = histogram2[:,:,membrane_coord]*valency
	ECpotential = Epotential+Cpotential
	reference = np.mean((ECpotential[0,0],ECpotential[-1,0],ECpotential[-1,-1],ECpotential[0,-1]))
	dEC=ECpotential-reference
	dEC_line = dEC[:,128]
	plt.plot(xcor,dEC_line,color=cm.jet(float(i)/plotcount,1))	
plt.legend(posilist)
plt.savefig("./view/ECpotential.png")
plt.close()
