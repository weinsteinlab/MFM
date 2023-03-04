import numpy as np
from gridData import OpenDX
from gridData import Grid
	
thickness = 10
charge_z = 40+128

origin = np.array([-128,-128,-128])
delta = np.diag((1,1,1))
ChargeDensity = np.zeros((257,257,257))
ChargeDensity[:,:,charge_z] = -0.005714286
#ChargeDensity[:,:,charge_z+thickness+1] = -0.005714286
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, ChargeDensity.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, ChargeDensity.shape))
dx.add('data', OpenDX.array(3, ChargeDensity))
dx.write("inputChargeDensity_memb.dx")


origin = np.array([-127.5,-128,-128])
delta = np.diag((1,1,1))
DialetricMapX = np.ones((257,257,257))*80 #78.54
DialetricMapX[:,:,charge_z+1:charge_z+1+thickness] = 2
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, DialetricMapX.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, DialetricMapX.shape))
dx.add('data', OpenDX.array(3, DialetricMapX))
dx.write("inputDialetricMapX_memb.dx")


origin = np.array([-128,-127.5,-128])
delta = np.diag((1,1,1))
DialetricMapY = np.ones((257,257,257))*80 #78.54
DialetricMapY[:,:,charge_z+1:charge_z+1+thickness] = 2
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, DialetricMapY.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, DialetricMapY.shape))
dx.add('data', OpenDX.array(3, DialetricMapY))
dx.write("inputDialetricMapY_memb.dx")


origin = np.array([-128,-128,-127.5])
delta = np.diag((1,1,1))
DialetricMapZ = np.ones((257,257,257))*80 #78.54
DialetricMapZ[:,:,charge_z:charge_z+1+thickness] = 2
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, DialetricMapZ.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, DialetricMapZ.shape))
dx.add('data', OpenDX.array(3, DialetricMapZ))
dx.write("inputDialetricMapZ_memb.dx")


origin = np.array([-128,-128,-128])
delta = np.diag((1,1,1))
IonAccess = np.ones((257,257,257))
IonAccess[:,:,charge_z+1:charge_z+1+thickness] = 0
dx = OpenDX.field('regular positions regular connections')
dx.add('positions', OpenDX.gridpositions(1, IonAccess.shape, origin, delta))
dx.add('connections', OpenDX.gridconnections(2, IonAccess.shape))
dx.add('data', OpenDX.array(3, IonAccess))
dx.write("inputIonAccess_memb.dx")


