import sys
import numpy as np
from gridData import OpenDX
from gridData import Grid
#import matplotlib.pyplot as plt


diel_bg = 80

dx1 = OpenDX.field(0)
dx1.read("../membranefolder/inputDialetricMapX.dx")
histogram1, edges = dx1.histogramdd()
xcor = (edges[0][1:]+edges[0][:-1])/2
ycor = (edges[1][1:]+edges[1][:-1])/2
zcor = (edges[2][1:]+edges[2][:-1])/2
origin = np.array([xcor[0],ycor[0],zcor[0]])
delta = np.diag((xcor[1]-xcor[0],ycor[1]-ycor[0],zcor[1]-zcor[0]))
dx2 = OpenDX.field(0)
dx2.read("../protein/dielx_proteinfilename.dx")
histogram2, edges = dx2.histogramdd()
diel = histogram1+histogram2-diel_bg
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, diel.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, diel.shape))
dx.add('data', OpenDX.array(3, diel))
dx.write("proteinfilename/cmbDialetricMapX.dx")

dx1 = OpenDX.field(0)
dx1.read("../membranefolder/inputDialetricMapY.dx")
histogram1, edges = dx1.histogramdd()
xcor = (edges[0][1:]+edges[0][:-1])/2
ycor = (edges[1][1:]+edges[1][:-1])/2
zcor = (edges[2][1:]+edges[2][:-1])/2
origin = np.array([xcor[0],ycor[0],zcor[0]])
delta = np.diag((xcor[1]-xcor[0],ycor[1]-ycor[0],zcor[1]-zcor[0]))
dx2 = OpenDX.field(0)
dx2.read("../protein/diely_proteinfilename.dx")
histogram2, edges = dx2.histogramdd()
diel = histogram1+histogram2-diel_bg
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, diel.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, diel.shape))
dx.add('data', OpenDX.array(3, diel))
dx.write("proteinfilename/cmbDialetricMapY.dx")

dx1 = OpenDX.field(0)
dx1.read("../membranefolder/inputDialetricMapZ.dx")
histogram1, edges = dx1.histogramdd()
xcor = (edges[0][1:]+edges[0][:-1])/2
ycor = (edges[1][1:]+edges[1][:-1])/2
zcor = (edges[2][1:]+edges[2][:-1])/2
origin = np.array([xcor[0],ycor[0],zcor[0]])
delta = np.diag((xcor[1]-xcor[0],ycor[1]-ycor[0],zcor[1]-zcor[0]))
dx2 = OpenDX.field(0)
dx2.read("../protein/dielz_proteinfilename.dx")
histogram2, edges = dx2.histogramdd()
diel = histogram1+histogram2-diel_bg
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, diel.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, diel.shape))
dx.add('data', OpenDX.array(3, diel))
dx.write("proteinfilename/cmbDialetricMapZ.dx")

kappa_bg = 1

dx1 = OpenDX.field(0)
dx1.read("../membranefolder/inputIonAccess.dx")
histogram1, edges = dx1.histogramdd()
xcor = (edges[0][1:]+edges[0][:-1])/2
ycor = (edges[1][1:]+edges[1][:-1])/2
zcor = (edges[2][1:]+edges[2][:-1])/2
origin = np.array([xcor[0],ycor[0],zcor[0]])
delta = np.diag((xcor[1]-xcor[0],ycor[1]-ycor[0],zcor[1]-zcor[0]))
dx2 = OpenDX.field(0)
dx2.read("../protein/kappa_proteinfilename.dx")
histogram2, edges = dx2.histogramdd()
kappa = histogram1+histogram2-kappa_bg
kappa = kappa.clip(min=0,max=1)
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, kappa.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, kappa.shape))
dx.add('data', OpenDX.array(3, kappa))
dx.write("proteinfilename/cmbIonAccess.dx")

charge_bg = 0

dx1 = OpenDX.field(0)
dx1.read("../membranefolder/inputChargeDensity.dx")
histogram1, edges = dx1.histogramdd()
xcor = (edges[0][1:]+edges[0][:-1])/2
ycor = (edges[1][1:]+edges[1][:-1])/2
zcor = (edges[2][1:]+edges[2][:-1])/2
origin = np.array([xcor[0],ycor[0],zcor[0]])
delta = np.diag((xcor[1]-xcor[0],ycor[1]-ycor[0],zcor[1]-zcor[0]))
dx2 = OpenDX.field(0)
dx2.read("../protein/charge_proteinfilename.dx")
histogram2, edges = dx2.histogramdd()
charge = histogram1+histogram2-charge_bg
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, charge.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, charge.shape))
dx.add('data', OpenDX.array(3, charge))
dx.write("proteinfilename/cmbChargeDensity.dx")
average_sigma = np.mean(histogram1[:,:,168])
print("average_sigma = %.15f" %average_sigma)
