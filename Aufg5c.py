# script test_readColumnData.py
# -*- coding: utf-8 -*-

import numpy as np
import PhyPraKit as ppk

# -----example Code illustrating usage --------------------
if __name__ == "__main__":
  import sys, numpy as np, matplotlib.pyplot as plt
  from PhyPraKit import odFit, labxParser, readColumnData
  from scipy import interpolate
  from scipy import signal

  def Imp(f,r,l,C):
      return np.sqrt(r**2+(2.*np.pi*f*l-1/(2.*np.pi*f*C))**2)
  
  def myFunction(n,As,wnull=9000,beta=1428.571):
      F = As/((wnull**2-(2*np.pi*n)**2)**2+(4*np.pi*n*beta)**2)**0.5
      return F

  fname='Aufg5_8_2Ohmb.dat'
  ncols=6

  data_array, info_dict =\
    readColumnData(fname, ncols, delimiter=' ', pr=False)

# print what we got:
  freqa=data_array[:,0] # 1st column
  U_La=data_array[:,1]
  U_Ca=data_array[:,2]
  Ueffa=data_array[:,3]
  phiva=data_array[:,5] # 5th column
  Stroma=data_array[:,4] # 4th column


  bname='Aufg5_47Ohm.dat'
  ncols=6

  data_array, info_dict =\
    readColumnData(bname, ncols, delimiter=' ', pr=False)

# print what we got:
  freqb=data_array[:,0] # 1st column
  U_Lb=data_array[:,1]
  U_Cb=data_array[:,2]
  Ueffb=data_array[:,3]
  phivb=data_array[:,5] # 5th column
  Stromb=data_array[:,4] # 4th column

  cname='Aufg5_100Ohm.dat'
  ncols=6

  data_array, info_dict =\
    readColumnData(cname, ncols, delimiter=' ', pr=False)

# print what we got:
  freqc=data_array[:,0] # 1st column
  U_Lc=data_array[:,1]
  U_Cc=data_array[:,2]
  Ueffc=data_array[:,3]
  phivc=data_array[:,5] # 5th column
  Stromc=data_array[:,4] # 4th column
  
  p0p=[500000000,9000.,1428.]
  
  Stroma=Stroma*1000.
  Stromb=Stromb*1000.
  Stromc=Stromc*1000.
  
  # ... and fit parameters 
  print "** envelope fit"
  parp, parpe, corp, chi2p =\
    odFit(myFunction, freqc,Stromc, 0., 0.1, p0=p0p)
  print "** fit of positive envelope, chi2/df= %.2g"%(chi2p/(len(freqa)-len(parp)))
  np.set_printoptions(precision=3)
  print " --> parameters:   ", parp
  np.set_printoptions(precision=2)
  print " --> uncertainties:", parpe 
  print " --> correlations:\n", corp


  
# plot data and analysis results
  
  fig=plt.figure()
  ax1=fig.add_subplot(1,1,1)
  #ax3.set_yticks([0,50,100,150])  ax3.plot(freqa, Stroma, 'b.')
  ax1.plot(freqa, U_La, 'b')
  ax1.plot(freqa, U_Ca, 'g')
  ax1.plot(freqa, Ueffa, 'k')
  ax1.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  ax1.set_ylabel('$Spannung$ U (V) ', size='large')
  ax1.grid()
 
  #plt.savefig('Aufgabe5cU08_2.pdf')
 
  fig2=plt.figure()
  ax2=fig2.add_subplot(1,1,1)
  ax2.plot(freqb, U_Lb, 'b')
  ax2.plot(freqb, U_Cb, 'g')
  ax2.plot(freqb, Ueffa, 'k')
  ax2.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  ax2.set_ylabel('$Spannung$ U (V) ', size='large')
  ax2.grid()
  
  #plt.savefig('Aufgabe5cU47.pdf')
  
  fig3=plt.figure()
  ax3=fig3.add_subplot(1,1,1)
  ax3.plot(freqc, U_Lc, 'b')
  ax3.plot(freqc, U_Cc, 'g')
  ax3.plot(freqc, Ueffc, 'k')
  ax3.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  ax3.set_ylabel('$Spannung$ U (V) ', size='large')
  ax3.grid()
  
  #plt.savefig('Aufgabe5cU100.pdf')
  
##################################################################
# ALLE Impedanz- Resonanz plot

  
# aus pylab_exampl