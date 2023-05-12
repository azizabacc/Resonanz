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


  # values for 400mA
  Aa=(0.094,0.139,0.202,0.523,0.431,0.156,0.075,0.059,0.052,0.047)
  fa=(0.21,0.32,0.41,0.51,0.61,0.7,0.8,0.86,0.92,0.97)
  
  # values for 600mA
  Ab=(0.104,0.137,0.202,0.372,0.261,0.116,0.082)
  fb=(0.23,0.33,0.43,0.51,0.61,0.7,0.8,0.86,0.92,0.97)
  
# define functional forms for envelope curves ...
  def myFunction(n,As,wnull,beta):
      F = As/((wnull**2-(2*np.pi*n)**2)**2+(4*np.pi*n*beta)**2)**0.5
      return F      
  
  p0p=[0.7, 13, 0.08] 
  
  '''
# ... and fit parameters 
  print "** envelope fit"
  parp, parpe, corp, chi2p =\
    odFit(myFunction, fb,Ab, 0., 0.1, p0=p0p)
  print "** fit of positive envelope, chi2/df= %.2g"%(chi2p/(len(fa)-len(parp)))
  np.set_printoptions(precision=3)
  print " --> parameters:   ", parp
  np.set_printoptions(precision=2)
  print " --> uncertainties:", parpe 
  print " --> correlations:\n", corp
  


  #pard, parde, cord, chi2d =\
  #  odFit(envd, td, phid, 0., 0.03, p0=p0d)
  #print "fit of negative envelope, chi2/df= %.2g"%(chi2d/(len(td)-len(pard)))
  #np.set_printoptions(precision=3)
  #print " -> parameters:   ", pard
  #np.set_printoptions(precision=2)
  #print " -> uncertainties:", parde 
  #print " -> correlations:\n", cord


  '''

# plot data and analysis results
  fig=plt.figure()
  ax3=fig.add_subplot(1,1,1)
  #ax3.set_yticks([0,50,100,150])  ax3.plot(freqa, Stroma, 'b.')
  ax3.plot(fb, Ab, 'r.')
  x=np.linspace(0,1,100)
  #ax3.plot(x,myFunction(x, *parp),'g')
  ax3.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  ax3.set_ylabel('$Amplitude$ (rad) ', size='large')
  #ax3.set_ylim(0,110)
  ax3.grid()
  
  
  
  
  
  #plt.savefig('Aufgabe4_600mA.pdf')
  



plt.show()

