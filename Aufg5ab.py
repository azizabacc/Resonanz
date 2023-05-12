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



  fname='Aufg5_8_2Ohm.dat'
  ncols=6

  data_array, info_dict =\
    readColumnData(fname, ncols, delimiter=' ', pr=False)

# print what we got:
  freqa=data_array[:,0] # 1st column
  phiva=data_array[:,5] # 5th column
  Stroma=data_array[:,4] # 4th column


  bname='Aufg5_47Ohm.dat'
  ncols=6

  data_array, info_dict =\
    readColumnData(bname, ncols, delimiter=' ', pr=False)

# print what we got:
  freqb=data_array[:,0] # 1st column
  phivb=data_array[:,5] # 5th column
  Stromb=data_array[:,4] # 4th column

  cname='Aufg5_100Ohm.dat'
  ncols=6

  data_array, info_dict =\
    readColumnData(cname, ncols, delimiter=' ', pr=False)

# print what we got:
  freqc=data_array[:,0] # 1st column
  phivc=data_array[:,5] # 5th column
  Stromc=data_array[:,4] # 4th column
  
  Stroma=Stroma*1000.
  Stromb=Stromb*1000.
  Stromc=Stromc*1000.

# define functional forms for envelope curves ...
  def myFunction(n,As,wnull,beta):
      F = As/((wnull-n**2)**2+(2*n*beta)**2)**0.5
      return F      
  
  p0p=[1., 1500., 500.] 
  
  
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
  ax3.plot(freqa, Stroma, 'b.')
  ax3.plot(freqa, Stroma, 'b')
  
  ax3.plot(freqb, Stromb, 'r.')
  ax3.plot(freqb, Stromb, 'r')
  ax3.plot(freqc, Stromc, 'g.')
  ax3.plot(freqc, Stromc, 'g')
  
  x=np.linspace(0,3000,3000)
  #ax3.plot(x,myFunction(x, *parp))
  #ax3.plot(x,myFunction(x,100,8000,8000),'g')
  ax3.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  ax3.set_ylabel('$Amplitude$ I (mA) ', size='large')
  ax3.set_ylim(0,110)
  ax3.grid()
  
  plt.savefig('Aufgabe5a.pdf')
  
  fig2=plt.figure()
  ax2=fig2.add_subplot(1,1,1)
  ax2.plot(freqa, phiva, 'b.')
  ax2.plot(freqa, phiva, 'b')
  
  ax2.plot(freqb, phivb, 'r.')
  ax2.plot(freqb, phivb, 'r')
  ax2.plot(freqc, phivc, 'g.')
  ax2.plot(freqc, phivc, 'g')
  
  x=np.linspace(0,3000,3000)
  #ax3.plot(x,myFunction(x, *parp))
  #ax3.plot(x,myFunction(x,100,8000,8000),'g')
  ax2.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  ax2.set_ylabel('$Phasenverschiebung$ (grad) ', size='large')
  ax2.grid()
  
  
  #plt.savefig('Aufgabe5ab.pdf')
  



#plt.show()

