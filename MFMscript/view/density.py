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
#valency = -1
#area_per_lipid = 65
#average_density = 0.2
boxsize=99
mapsize=257

posilist = 'posi_list'
posilist = posilist.split(' ')
posilist.reverse()

plotcount = len(posilist)
v_count = int(np.ceil(np.sqrt(plotcount)))
h_count = int(np.ceil(plotcount/v_count))
fig, ax = plt.subplots(h_count, v_count,figsize=(10, 7.5))

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


dx1 = OpenDX.field(0)
for i in range(plotcount):
	dx1.read("./%s/cmbChargeDensity.dx" %posilist[i])
	histogram1, edges = dx1.histogramdd()
	xcor = (edges[0][1:]+edges[0][:-1])/2
	ycor = (edges[1][1:]+edges[1][:-1])/2
	zcor = (edges[2][1:]+edges[2][:-1])/2
	membxy = [xcor[0],xcor[-1],ycor[0],ycor[-1]]
	membrane_coord = np.where(zcor==membrane_position)[0][0]
	density_map = histogram1[:,:,membrane_coord]*area_per_lipid/valency
	singularity = np.where(density_map<0)
	if i==0:
		norm = MidpointNormalize(midpoint=average_density)
		norm.autoscale(density_map)
	pcm = ax.flatten()[-1-i].imshow(density_map,norm=norm,cmap='RdBu_r',origin='lower',extent=membxy)
	for i in range(len(singularity[0])):
		ax.flatten()[-1-i].plot(singularity[1][i],singularity[0][i],'r.',markersize=1)
	if 0==0:
		extreme_point = np.argmax(density_map)
		extreme_point_x = np.floor(extreme_point/len(ycor))
		extreme_point_y = extreme_point-extreme_point_x*len(ycor)
		#ax.flatten()[-1-i].plot(extreme_point_y,extreme_point_x,'k.')
		ax.flatten()[-1-i].text(extreme_point_y+ycor[0],extreme_point_x+xcor[0],'%0.6f' %density_map[int(extreme_point_x),int(extreme_point_y)],color='r')
		extreme_point = np.argmin(density_map)
		extreme_point_x = np.floor(extreme_point/len(ycor))
		extreme_point_y = extreme_point-extreme_point_x*len(ycor)
		#ax.flatten()[-1-i].plot(extreme_point_y,extreme_point_x,'k.')
		ax.flatten()[-1-i].text(extreme_point_y+ycor[0],extreme_point_x+xcor[0],'%0.6f' %density_map[int(extreme_point_x),int(extreme_point_y)],color='b')
		ax.flatten()[-1-i].plot((mapsize-boxsize)/2+np.array(range(boxsize))+ycor[0],np.ones(boxsize)*(mapsize-boxsize)/2+xcor[0],'k-',linewidth=1)
		ax.flatten()[-1-i].plot((mapsize-boxsize)/2+np.array(range(boxsize))+ycor[0],np.ones(boxsize)*(mapsize+boxsize)/2+1+xcor[0],'k-',linewidth=1)
		ax.flatten()[-1-i].plot(np.ones(boxsize)*(mapsize-boxsize)/2+ycor[0],(mapsize-boxsize)/2+np.array(range(boxsize))+xcor[0],'k-',linewidth=1)
		ax.flatten()[-1-i].plot(np.ones(boxsize)*(mapsize+boxsize)/2+1+ycor[0],(mapsize-boxsize)/2+np.array(range(boxsize))+xcor[0],'k-',linewidth=1)	
	ax.flatten()[-1-i].text(0,len(ycor)/2,'step %s' %posilist[i], color='k',fontsize=12,bbox=dict(facecolor='w', edgecolor='k',boxstyle='round'), horizontalalignment='center')
		
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(pcm,cax=cbar_ax)	
plt.savefig("./view/Density.png")
