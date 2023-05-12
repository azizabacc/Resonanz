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
  
  
  # define functional forms for envelope curves ...
  def myFunction(n,faktor,As,wnull,beta):
      F =faktor*(As/((wnull-(2*np.pi*n)**2)**2+(2*2*np.pi*n*beta)**2)**0.5)
      return F     
  
  
  
  x=np.linspace(0,1000,10000)
  plt.plot(x,myFunction(x,1, 20,9000,1),'b')
  plt.plot(x,myFunction(x,1, 20,9000,5),'g')
  plt.plot(x,myFunction(x,1, 20,9000,10),'r')
  plt.plot(x,myFunction(x,1, 20,9000,10),'y')
  plt.show()