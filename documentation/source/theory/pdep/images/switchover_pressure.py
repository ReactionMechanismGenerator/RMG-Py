#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

Tlist = numpy.arange(300.0, 2201.0, 100.0, numpy.float64)

Plist = numpy.zeros((len(Tlist),5), numpy.float64)

#                          300   400   500   600   700   800   900  1000  1100  1200  1300  1400  1500  1600  1700  1800  1900  2000  2100  2200
# m = 27
Plist[:,0] = numpy.array([-1.1, -0.9, -0.6, -0.4, -0.2,  0.0,  0.2,  0.5,  0.8,  1.1,  1.4,  1.7,  1.9,  2.1,  2.3,  2.5,  2.7,  2.9,  3.0,  3.1], numpy.float64)
# m = 39
Plist[:,1] = numpy.array([-1.6, -1.4, -1.2, -1.0, -0.7, -0.4, -0.1,  0.2,  0.5,  0.8,  1.1,  1.4,  1.7,  1.9,  2.1,  2.3,  2.5,  2.7,  2.9,  3.0], numpy.float64)
# m = 60
Plist[:,2] = numpy.array([-2.2, -2.0, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3,  0.1,  0.5,  0.8,  1.1,  1.4,  1.6,  1.9,  2.1,  2.3,  2.5,  2.7,  2.9], numpy.float64)
# m = 102
Plist[:,3] = numpy.array([-2.9, -2.7, -2.4, -2.1, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.5,  0.8,  1.1,  1.4,  1.7,  2.0,  2.2,  2.4,  2.6,  2.8], numpy.float64)
# m = 1000
Plist[:,4] = numpy.array([-4.1, -3.8, -3.5, -3.1, -2.7, -2.2, -1.7, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.1,  1.4,  1.7,  2.0,  2.2,  2.4,  2.6], numpy.float64)



poly = numpy.polyfit(Tlist, Plist, deg=4)


Tlist2 = numpy.arange(300.0, 2201.0, 10.0, numpy.float64)
Plist2 = numpy.zeros((len(Tlist2),5), numpy.float64)
for i in range(5):
    Plist2[:,i] = numpy.polyval(poly[:,i], Tlist2)

Plist = 10**Plist * 101325
Plist2 = 10**Plist2 * 101325

import pylab
fig = pylab.figure(figsize=(6,5))
ax = fig.add_subplot(1,1,1)

linespec = ['-b', '-g', '-r', '-m', '-k']
for i in range(5):
    pylab.semilogy(Tlist2, Plist2[:,i] / 1e5, linespec[i])
pylab.xlim(300, 2200)
pylab.ylim(1e-5, 1e4)
pylab.xlabel('Temperature (K)')
pylab.ylabel('Pressure (bar)')

import matplotlib.ticker
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=500.0))
ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=100.0))

fig.subplots_adjust(left=0.14, bottom=0.12, right=0.95, top=0.95, wspace=0.20, hspace=0.20)

pylab.annotate('m = 27', xy=(Tlist2[60],Plist2[60,0] / 1e5), xytext=(500,1e1), arrowprops=dict(arrowstyle="->", color='b'), color='b')
pylab.annotate('m = 39', xy=(Tlist2[100],Plist2[100,1] / 1e5), xytext=(1000,1e2), arrowprops=dict(arrowstyle="->", color='g'), color='g')
pylab.annotate('m = 60', xy=(Tlist2[120],Plist2[120,2] / 1e5), xytext=(1500,1e0), arrowprops=dict(arrowstyle="->", color='r'), color='r')
pylab.annotate('m = 102', xy=(Tlist2[80],Plist2[80,3] / 1e5), xytext=(1000,1e-2), arrowprops=dict(arrowstyle="->", color='m'), color='m')
pylab.annotate('m = 1000', xy=(Tlist2[40],Plist2[40,4] / 1e5), xytext=(500,1e-4), arrowprops=dict(arrowstyle="->", color='k'), color='k')

pylab.savefig('switchover_pressure.pdf')
pylab.savefig('switchover_pressure.png')

pylab.show()

