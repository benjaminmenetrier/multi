#!/usr/bin/env python3

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os
from scipy import interpolate

n1 = 10
x1 = []
y1 = []
for i in range (0,n1+1):
   x1.append(float(i)/float(n1))
   y1.append(0.1*math.sin(2.0*math.pi*float(i)/float(n1)) \
            -0.4*math.sin(4.0*math.pi*float(i)/float(n1)) \
            +0.6*math.sin(6.0*math.pi*float(i)/float(n1)) )
lin1 = interpolate.interp1d(x1, y1, kind='linear')
nn1 = interpolate.interp1d(x1, y1, kind='nearest')

n2 = 14
x2 = []
ysp2 = []
ylin2 = []
ynn2 = []
for i in range (0,n2+1):
   x2.append(float(i)/float(n2))
   ysp2.append(0.1*math.sin(2.0*math.pi*float(i)/float(n2)) \
              -0.4*math.sin(4.0*math.pi*float(i)/float(n2)) \
              +0.6*math.sin(6.0*math.pi*float(i)/float(n2)) )
   ylin2.append(lin1(float(i)/float(n2)))
   ynn2.append(nn1(float(i)/float(n2)))
lin2 = interpolate.interp1d(x2, ylin2, kind='linear')
nn2 = interpolate.interp1d(x2, ynn2, kind='nearest')
 

n3 = 18
x3 = []
ysp3 = []
ylin3 = []
ynn3 = []
ylindir = []
ynndir = []
for i in range (0,n3+1):
   x3.append(float(i)/float(n3))
   ysp3.append(0.1*math.sin(2.0*math.pi*float(i)/float(n3)) \
              -0.4*math.sin(4.0*math.pi*float(i)/float(n3)) \
              +0.6*math.sin(6.0*math.pi*float(i)/float(n3)) )
   ylin3.append(lin2(float(i)/float(n3)))
   ynn3.append(nn2(float(i)/float(n3)))
   ylindir.append(lin1(float(i)/float(n3)))
   ynndir.append(nn1(float(i)/float(n3)))

nref = 100
xref = []
yspref = []
ylin1ref = []
ylin2ref = []
ynn1ref = []
ynn2ref = []
for i in range (0,nref+1):
   xref.append(float(i)/float(nref))
   yspref.append(0.1*math.sin(2.0*math.pi*float(i)/float(nref)) \
                -0.4*math.sin(4.0*math.pi*float(i)/float(nref)) \
                +0.6*math.sin(6.0*math.pi*float(i)/float(nref)) )
   ylin1ref.append(lin1(float(i)/float(nref)))
   ylin2ref.append(lin2(float(i)/float(nref)))
   ynn1ref.append(nn1(float(i)/float(nref)))
   ynn2ref.append(nn2(float(i)/float(nref)))


fig, ax = plt.subplots(5,figsize=(4,12))
fig.subplots_adjust(hspace=0.5)
ax[0].set_title("Low-resolution function")
ax[0].plot(x1, y1, 'r.')
ax[0].set(xlim=(0, 1))
ax[1].set_title("Spectral upscaling 1 -> 2")
ax[1].plot(x1, y1, 'r.')
ax[1].plot(x2, ysp2, 'b.')
ax[1].plot(xref, yspref, 'k-',linewidth=0.3)
ax[1].set(xlim=(0, 1))
ax[2].set_title("Spectral upscaling 2 -> 3")
ax[2].plot(x2, ysp2, 'b.')
ax[2].plot(x3, ysp3, 'g.')
ax[2].plot(xref, yspref, 'k-',linewidth=0.3)
ax[2].set(xlim=(0, 1))
ax[3].set_title("Spectral upscaling 1 -> 3")
ax[3].plot(x1, y1, 'r.')
ax[3].plot(x3, ysp3, 'g*')
ax[3].plot(xref, yspref, 'k-',linewidth=0.3)
ax[3].set(xlim=(0, 1))
ax[4].set_title("Comparison at scale 3")
ax[4].plot(x3, ysp3, 'g.')
ax[4].plot(x3, ysp3, 'g*')
ax[4].set(xlim=(0, 1))
plt.savefig("interpolation_comparison_sp.pdf", format="pdf", dpi=300)
plt.close()
os.system('pdfcrop interpolation_comparison_sp.pdf interpolation_comparison_sp.pdf')

fig, ax = plt.subplots(5,figsize=(4,12))
fig.subplots_adjust(hspace=0.5)
ax[0].set_title("Low-resolution function")
ax[0].plot(x1, y1, 'r.')
ax[0].set(xlim=(0, 1))
ax[1].set_title("Linear upscaling 1 -> 2")
ax[1].plot(x1, y1, 'r.')
ax[1].plot(x2, ylin2, 'b.')
ax[1].plot(xref, ylin1ref, 'k-',linewidth=0.3)
ax[1].set(xlim=(0, 1))
ax[2].set_title("Linear upscaling 2 -> 3")
ax[2].plot(x2, ylin2, 'b.')
ax[2].plot(x3, ylin3, 'g.')
ax[2].plot(xref, ylin2ref, 'k-',linewidth=0.3)
ax[2].set(xlim=(0, 1))
ax[3].set_title("Linear upscaling 1 -> 3")
ax[3].plot(x1, y1, 'r.')
ax[3].plot(x3, ylindir, 'g*')
ax[3].plot(xref, ylin1ref, 'k-',linewidth=0.3)
ax[3].set(xlim=(0, 1))
ax[4].set_title("Comparison at scale 3")
ax[4].plot(x3, ylin3, 'g.')
ax[4].plot(x3, ylindir, 'g*')
ax[4].set(xlim=(0, 1))
plt.savefig("interpolation_comparison_lin.pdf", format="pdf", dpi=300)
plt.close()
os.system('pdfcrop interpolation_comparison_lin.pdf interpolation_comparison_lin.pdf')

fig, ax = plt.subplots(5,figsize=(4,12))
fig.subplots_adjust(hspace=0.5)
ax[0].set_title("Low-resolution function")
ax[0].plot(x1, y1, 'r.')
ax[0].set(xlim=(0, 1))
ax[1].set_title("Nearest neighbor upscaling 1 -> 2")
ax[1].plot(x1, y1, 'r.')
ax[1].plot(x2, ynn2, 'b.')
ax[1].plot(xref, ynn1ref, 'k-',linewidth=0.3)
ax[1].set(xlim=(0, 1))
ax[2].set_title("Nearest neighbor upscaling 2 -> 3")
ax[2].plot(x2, ynn2, 'b.')
ax[2].plot(x3, ynn3, 'g.')
ax[2].plot(xref, ynn2ref, 'k-',linewidth=0.3)
ax[2].set(xlim=(0, 1))
ax[3].set_title("Nearest neighbor upscaling 1 -> 3")
ax[3].plot(x1, y1, 'r.')
ax[3].plot(x3, ynndir, 'g*')
ax[3].plot(xref, ynn1ref, 'k-',linewidth=0.3)
ax[3].set(xlim=(0, 1))
ax[4].set_title("Comparison at scale 3")
ax[4].plot(x3, ynn3, 'g.')
ax[4].plot(x3, ynndir, 'g*')
ax[4].set(xlim=(0, 1))
plt.savefig("interpolation_comparison_nn.pdf", format="pdf", dpi=300)
plt.close()
os.system('pdfcrop interpolation_comparison_nn.pdf interpolation_comparison_nn.pdf')
