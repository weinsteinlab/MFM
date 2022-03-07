import numpy as np
from gridData import OpenDX
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib
matplotlib.use('Agg')

membrane_position = 40
valency = -4
area_per_lipid = 70
average_density = 0.1
valency = -1
area_per_lipid = 65
#average_density = 0.2
boxsize=99
mapsize=257

posilist = 'posi_list'
posilist = posilist.split(' ')
path = 'foldername'
plotcount=len(posilist)

dx1 = OpenDX.field(0)
for i in range(plotcount):
	dx1.read("./%s/cmbChargeDensity.dx" %posilist[i])
	histogram1, edges = dx1.histogramdd()
	xcor = (edges[0][1:]+edges[0][:-1])/2
	ycor = (edges[1][1:]+edges[1][:-1])/2
	zcor = (edges[2][1:]+edges[2][:-1])/2
	membrane_coord = np.where(zcor==membrane_position)[0][0]
	density_map = histogram1[:,:,membrane_coord]*area_per_lipid/valency
	np.savetxt('view/density_%s_%s.dat'%(path,posilist[i]),density_map,fmt='%1.4f')

