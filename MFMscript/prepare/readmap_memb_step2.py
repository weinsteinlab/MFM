import numpy as np
from gridData import OpenDX
from gridData import Grid
#import matplotlib.pyplot as plt


dx = OpenDX.field(0)
dx.read("Potential_step1.dx")

histogram, edges = dx.histogramdd()

xcor = (edges[0][1:]+edges[0][:-1])/2
ycor = (edges[1][1:]+edges[1][:-1])/2
zcor = (edges[2][1:]+edges[2][:-1])/2


midline = histogram[128,128,:]
#boundary = histogram[0,128,:]

writedata = np.zeros((257,5))
writedata[:,2] = zcor
writedata[:,3] = midline
np.savetxt('boundary.txt', writedata, fmt='%.6f',delimiter='  ')

#writedata = np.zeros((257,5))
#writedata[:,2] = zcor
#writedata[:,3] = boundary
#np.savetxt('testboundary.txt', writedata, fmt='%.6f',delimiter='  ')



'''
dx = OpenDX.field(0)
dx.read("thick_b/Potential.dx")
histogram_thick_b, edges = dx.histogramdd()
dx = OpenDX.field(0)
dx.read("thin_b/Potential.dx")
histogram_thin_b, edges = dx.histogramdd()
dx = OpenDX.field(0)
dx.read("no_boundary/thick/Potential.dx")
histogram_thick, edges = dx.histogramdd()
dx = OpenDX.field(0)
dx.read("no_boundary/thin/Potential.dx")
histogram_thin, edges = dx.histogramdd()
'''
'''
midline_tkb = histogram_thick_b[128,128,:]
midline_tnb = histogram_thin_b[128,128,:]
midline_tk = histogram_thick[128,128,:]
midline_tn = histogram_thin[128,128,:]
plt.figure(figsize=(10,7))
plt.plot(zcor,midline_tn,'m',alpha=0.5)
plt.plot(zcor,midline_tnb,'r',alpha=0.5)
plt.plot(zcor,midline_tk,'b',alpha=0.5)
plt.plot(zcor,midline_tkb,'g',alpha=0.5)
plt.legend(('10A','10A-bound','40A','40A-bound'))
plt.grid()
plt.savefig('potential_at_mid.jpg')
'''
'''
sideline_tk = histogram_thick[0,128,:]
sideline_tkb = histogram_thick_b[0,128,:]
sideline_tn = histogram_thin[0,128,:]
sideline_tnb = histogram_thin_b[0,128,:]
plt.figure(figsize=(10,7))
plt.plot(zcor,sideline_tn,'m',alpha=0.5)
plt.plot(zcor,sideline_tnb,'r',alpha=0.5)
plt.plot(zcor,sideline_tk,'b',alpha=0.5)
plt.plot(zcor,sideline_tkb,'g',alpha=0.5)
plt.legend(('10A-nobound','10A-bound','40A-nobound','40A-bound'))
plt.grid()
plt.savefig('potential_at_side.jpg')
'''
