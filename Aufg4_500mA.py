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


  # vlalues for 300mA
  Aa=(0.08,0.09,0.115,0.13,0.175,0.26,0.48,0.90,1.045,1.175,1.14,1.08,0.915,0.29,0.16)
  fa=(0.264,0.304,0.353,0.393,0.44,0.49,0.525,0.542,0.55,0.56,0.56,0.563,0.568,0.618,0.67)
  
  # values for 500mA
  Ab=(0.08,0.09,0.105,0.13,0.17,0.23,0.345,0.43,0.43,0.435,0.43,0.43,0.41,0.24,0.135)
  fb=(0.24,0.309,0.352,0.393,0.439,0.483,0.527,0.544,0.554,0.555,0.556,0.563,0.570,0.613,0.66)

# define functional forms for envelope curves ...
  def myFunction(n,As,wnull,beta):
      F = As/((wnull-(2*np.pi*n)**2)**2+(4*np.pi*n*beta)**2)**0.5
      return F      
  
  p0p=[0.7, 13, 0.08] 
  
  
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




# plot data and analysis results
  fig=plt.figure()
  ax3=fig.add_subplot(1,1,1)
  #ax3.set_yticks([0,50,100,150])  ax3.plot(freqa, Stroma, 'b.')
  ax3.plot(fb, Ab, 'r.')
  x=np.linspace(0,1,100)
  ax3.plot(x,myFunction(x, *parp),'g')
  ax3.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  ax3.set_ylabel('$Amplitude$ (rad) ', size='large')
  #ax3.set_ylim(0,110)
  ax3.grid()
  
  
  
  
  
  #plt.savefig('Aufgabe4_500mA.pdf')
  



  plt.show()

