import numpy as np
from gridData import OpenDX
from gridData import Grid
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
matplotlib.use('Agg')

membrane_position = 40
valency = -4
#valency = -3
area_per_lipid = 70
everage_density = 0.1
#everage_density = 0.025*70/65
#valency = -1
#area_per_lipid = 65
#everage_density = 0.2

posilist = 'posi_list'
posilist = posilist.split(' ')

dx1 = OpenDX.field(0)
dx2 = OpenDX.field(0)

plotcount = len(posilist)
v_count = int(np.ceil(np.sqrt(plotcount)))
h_count = int(np.ceil(plotcount/v_count))
fig, ax = plt.subplots(h_count, v_count, figsize=(10, 7.5))

for i in range(plotcount):
	dx1.read("./%s/cmbChargeDensity.dx" %posilist[i])
	histogram1, edges = dx1.histogramdd()
	xcor = (edges[0][1:]+edges[0][:-1])/2
	ycor = (edges[1][1:]+edges[1][:-1])/2
	zcor = (edges[2][1:]+edges[2][:-1])/2
	membrane_coord = np.where(zcor==membrane_position)[0][0]
	membxy = [xcor[0],xcor[-1],ycor[0],ycor[-1]]
	density_map = histogram1[:,:,membrane_coord]*area_per_lipid/valency
	Cpotential = np.log(density_map*(1-everage_density)/((1-density_map)*everage_density))
	dx2.read("./%s/Potential.dx" %posilist[i])
	histogram2, edges = dx2.histogramdd()
	Epotential = histogram2[:,:,membrane_coord]*valency
	ECpotential = Epotential+Cpotential
	reference = np.mean((ECpotential[0,0],ECpotential[-1,0],ECpotential[-1,-1],ECpotential[0,-1]))
	dEC=ECpotential-reference
	if i==0:
		linthresh = min(abs(np.min(dEC)),abs(np.max(dEC)))/30
		vmax = max(abs(np.min(dEC)),abs(np.max(dEC)))
		norm = colors.SymLogNorm(linthresh=linthresh,vmin=-vmax,vmax=vmax)
	pcm = ax.flatten()[i].imshow(dEC,norm=norm,cmap='RdBu_r',origin='lower',extent=membxy)

	extreme_point = np.argmax(abs(dEC))
	extreme_point_x = np.floor(extreme_point/len(ycor))
	extreme_point_y = extreme_point-extreme_point_x*len(ycor)
	ax.flatten()[i].plot(extreme_point_y+ycor[0],extreme_point_x+xcor[0],'k.')
	ax.flatten()[i].text(extreme_point_y+ycor[0],extreme_point_x+xcor[0],'%2.5f' %dEC[int(extreme_point_x),int(extreme_point_y)])

	ax.flatten()[i].text(0,len(xcor)/2,'step %s' %posilist[i], color='k',fontsize=12,bbox=dict(facecolor='w', edgecolor='k',boxstyle='round'), horizontalalignment='center')
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(pcm,cax=cbar_ax,extend='both')	

fig.savefig("./view/ECpotential2d.png")
