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



  fname='Aufg2_400mA.dat'
  ncols=2

  data_array, info_dict =\
    readColumnData(fname, ncols, delimiter=' ', pr=False)

# print what we got:
  t=data_array[:,0] # 1st column
  phi=data_array[:,1] # 2nd column

#  print "Title=", info_dict['*TITLE']
  print "t= ", t
  print "phi= ", phi

    
  print "keywords found:"
  for key in info_dict:
    if (info_dict[key]!=None): print key,':', info_dict[key]

# special line for t-omega diagramm
  teins=t
# some filtering and smoothing
  print "** filtering:"
  # * use, if there is an unwanted offset
  print "  offset correction"
  phi_c = ppk.offsetFilter(phi)          

  # * use, if sampling rate is too high  
  if len(phi)>2000: 
    print "  resampling"
    phi, t = ppk.resample(phi, t, n=int(len(phi)/1000)) # average n samples into 1




# calculate fourier spectrum 
  print "** Fourier Spectrum"
 # freq, amp = FourierSpectrum(t, phi, fmax=1.)
  freq, amp = ppk.Fourier_fft(t, phi)  # use fast algorithm
  frequency = freq[np.where(amp==max(amp))]
  print " --> Frequenz: ",frequency

# run a peak finder
  # first, determine width of typical peaks and dips 
  width=  0.5* len(t)  / (t[-1] - t[0]) / frequency
#  use convoluted template filter 
  peakind = ppk.convolutionPeakfinder(phi, width, th=0.53)
  dipind  = ppk.convolutionPeakfinder(-phi, width, th=0.53)
  if len(peakind) > 5:
    print " --> %i peaks and %i dips found"%(len(peakind), len(dipind))
    tp, phip = np.array(t[peakind]), np.array(phi[peakind])
    td, phid =np.array(t[dipind]), np.array(phi[dipind])
  else:
    print "*!!* not enough peaks found for envelope fit"
    print "     tune peakfinder parameters!"
    sys.exit(1)
    print "*!!* not enough peaks found for envelope fit"
    print "     tune peakfinder parameters!"
    sys.exit(1)
    
    
      
# extra peaks and dips found per eye
  nfoundp=6
  nfoundd=0
  tpf=(12.85,14.6,16.4,18.2,19.95,21.75)
  phipf=(0.30,0.22,0.15,0.11,0.08,0.04)
  #tdf=(175.09,176.86,178.69,180.51,182.31,184.06,185.86,187.66,189.46,191.21,192.96,194.81,196.61,198.41,200.16,201.96,203.76,205.87,207.56,209.79,211,213,215,217,219)
  #phidf=(-0.2,-0.17,-0.17,-0.15,-0.13,-0.13,-0.12,-0.1,-0.09,-0.08,-0.07,-0.05,-0.04,-0.04,-0.03,-0.02,-0.01,0,0,0,0,0,0,0,0)

# new dips and peaks
  tpn=np.zeros(len(tp)+nfoundp)
  tdn=np.zeros(len(td)+nfoundd)
  phipn=np.zeros(len(phip)+nfoundp)
  phidn=np.zeros(len(phid)+nfoundd)
  i=0
  while i<len(tp):
      tpn[i]=tp[i]
      phipn[i]=phip[i]
      i=i+1
  else:
      j=0
      while j<nfoundp:
          tpn[len(tp)+j]=tpf[j]
          phipn[len(phip)+j]=phipf[j]
          j=j+1

  i=0
  while i<len(td):
      tdn[i]=td[i]
      phidn[i]=phid[i]
      i=i+1
  else:
      j=0
      while j<nfoundd:
          tdn[len(td)+j]=tdf[j]
          phidn[len(phid)+j]=phidf[j]
          j=j+1
  # extra line because of wrong first number
  tpn[0]=0.15
  phipn[0]=2.82
  
  #print tpn
  #print tdn
  

# cubic spline interpolation with scipy.interpolate
  print "** spline interpolation"
  cs_phi=interpolate.UnivariateSpline(t, phi, s=0)
  cs_omega=cs_phi.derivative()  

# def Ekin
  def Ekin(s,T=1.39e-03):
      return 0.5*T*s**2


# define functional forms for envelope curves ...
  def env_exp_p(t, A=1., tau=75.):
    return A * np.exp(-t/tau) 

  def env_exp_d(t, A=2., tau=75.):
    return -A * np.exp(-t/tau) 

  def env_quad_p(t, a=1., b=1., c=1.):
    return a*t**2+b*t+c

  def env_quad_d(t, a=-1., b=1., c=1.):
    return a*t**2+b*t+c

  def env_linexp_p(t, A=2., tau=75., a=0.001):
    return  A*(np.exp(-t/tau) + a) 

  def env_linexp_d(t, A=2., tau=75., a=0.001):
    return  -A*(np.exp(-t/tau) + a)


# select functions for envelope ...
#  !!! Achtung: nichlineare Anpassungen benötigen Startwerte,
#           also Schätzungen der Parameterwerte
  envp=env_linexp_p     # exp + linear
  p0p=[2., 100., 0.001] # initial parameters 
  envd=env_linexp_d
  p0d=[2., 100, -0.01]  
  
#  envp=env_exp_p     # exp + linear
#  p0p=[2., 100.] # initial parameters 
#  envd=env_exp_d
#  p0d=[-2., 100.]  

  
 
# ... and fit parameters 
  print "** envelope fit"
  parp, parpe, corp, chi2p =\
    odFit(envp, tpn, phipn, 0., 0.1, p0=p0p)
  print "** fit of positive envelope, chi2/df= %.2g"%(chi2p/(len(tp)-len(parp)))
  np.set_printoptions(precision=3)
  print " --> parameters:   ", parp
  np.set_printoptions(precision=2)
  print " --> uncertainties:", parpe 
  print " --> correlations:\n", corp

  pard, parde, cord, chi2d =\
    odFit(envd, td, phid, 0., 0.03, p0=p0d)
  print "fit of negative envelope, chi2/df= %.2g"%(chi2d/(len(td)-len(pard)))
  np.set_printoptions(precision=3)
  print " -> parameters:   ", pard
  np.set_printoptions(precision=2)
  print " -> uncertainties:", parde 
  print " -> correlations:\n", cord

# Finde schnittpunkt der Einhuellenden
  
  begr=1
  i=0
  xend=np.linspace(0., 300,1000)
  while begr>0.001:
      begr=np.sqrt((envp(xend[i], *parp)-envd(xend[i], *pard))**2)
      i=i+1
  xend=xend[i]


# plot data and analysis results
  fig=plt.figure("Amplitude")
  #fig.suptitle('Aufgabe1 freie Schwingung', size='x-large', color='b')
  fig.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.93,
                    wspace=None, hspace=.25)#

  ax1=fig.add_subplot(3,1,1)
  plt.plot(t, phi)
  plt.plot(tpn, phipn, 'rx', alpha=0.5, label='peak')
  plt.plot(td, phid, 'gx', alpha=0.5, label='dip')
  x=np.linspace(0., xend, 100)
  plt.plot(x, envp(x, *parp), 
    'r-', linewidth=1, label=u'positive Einhüllende')
  plt.plot(x, envd(x, *pard), 
    'g-', linewidth=1, label=u'negative Einhüllende')
  plt.xlabel('$Zeit$ $t$ (s)', size='large')
  plt.ylabel('$Winkel$  $\phi$', size='large')
  plt.legend(loc='best', numpoints=1, prop={'size':10})
  plt.grid()

#  ax2=fig.add_subplot(3,1,2)
#  ax2.plot(teins , omega)
#  ax2.set_xlabel('$Zeit$ $t$ (s)', size='large')
#  ax2.set_ylabel('$Winkelgeschwindigkeit$  $\omega$ (1/s)', size='large')
#  ax2.grid()


  #plt.savefig('Aufgabe2400uebmAueb.pdf')
  plt.show()

  print 'tpn', phipn