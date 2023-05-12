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



  fname='Aufg1.dat'
  ncols=3

  data_array, info_dict =\
    readColumnData(fname, ncols, delimiter=' ', pr=False)

# print what we got:
  t=data_array[:,0] # 1st column
  phi=data_array[:,1] # 2nd column
  omega=data_array[:,2] # 3rd column

  print "Title=", info_dict['*TITLE']
  print "t= ", t
  print "phi= ", phi

    
  print "keywords found:"
  for key in info_dict:
    if (info_dict[key]!=None): print key,':', info_dict[key]

# special line for t-omega diagramm
  teins=t
  phieins=phi
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

# extra peaks and dips found per eye
  nfoundp=12+8
  nfoundd=17+8
  tpf=(184.96,186.76,188.56,190.31,192.11,193.91,195.66,197.51,199.26,201.06,202.86,204.61,206.5,208.67,210.46,212,214,216,218,220)
  phipf=(0.13,0.12,0.1,0.09,0.08,0.07,0.05,0.04,0.04,0.03,0.02,0.01,0,0,0,0,0,0,0,0)
  tdf=(175.09,176.86,178.69,180.51,182.31,184.06,185.86,187.66,189.46,191.21,192.96,194.81,196.61,198.41,200.16,201.96,203.76,205.87,207.56,209.79,211,213,215,217,219)
  phidf=(-0.2,-0.17,-0.17,-0.15,-0.13,-0.13,-0.12,-0.1,-0.09,-0.08,-0.07,-0.05,-0.04,-0.04,-0.03,-0.02,-0.01,0,0,0,0,0,0,0,0)

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
    return  A*(np.exp(-t/tau) - a*t) 

  def env_linexp_d(t, A=2., tau=75., a=0.001):
    return  -A*(np.exp(-t/tau) - a*t)


# select functions for envelope ...
#  !!! Achtung: nichlineare Anpassungen benötigen Startwerte,
#           also Schätzungen der Parameterwerte
  envp=env_linexp_p     # exp + linear
  p0p=[2., 100., 0.001] # initial parameters 
  envd=env_linexp_d
  p0d=[-2., 100., 0.001]  
  
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
    odFit(envd, tdn, phidn, 0., 0.03, p0=p0d)
  print "fit of negative envelope, chi2/df= %.2g"%(chi2d/(len(td)-len(pard)))
  np.set_printoptions(precision=3)
  print " -> parameters:   ", pard
  np.set_printoptions(precision=2)
  print " -> uncertainties:", parde 
  print " -> correlations:\n", cord

# plot data and analysis results
  fig=plt.figure("Amplitude", figsize=(7.5, 10.))
  #fig.suptitle('Aufgabe1 freie Schwingung', size='x-large', color='b')
  fig.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.93,
                    wspace=None, hspace=.25)#

  ax1=fig.add_subplot(3,1,1)
  #ax1.plot(t, phi)
  ax1.plot(teins,phieins)
  #ax1.plot(tpn, phipn, 'rx', alpha=0.5, label='peak')
  #ax1.plot(tdn, phidn, 'gx', alpha=0.5, label='dip')
  x=np.linspace(0., 200, 100)
  ax1.plot(x, envp(x, *parp), 
    'r-', linewidth=1, label=u'positive Einhüllende')
  ax1.plot(x, envd(x, *pard), 
    'g-', linewidth=1, label=u'negative Einhüllende')
  ax1.set_xlabel('$Zeit$ $t$ (s)', size='large')
  ax1.set_ylabel('$Winkel$  $\phi$', size='large')
  ax1.legend(loc='best', numpoints=1, prop={'size':10})
  ax1.grid()

  ax2=fig.add_subplot(3,1,2)
  ax2.plot(teins , omega)
  ax2.set_xlabel('$Zeit$ $t$ (s)', size='large')
  ax2.set_ylabel('$Winkelgeschwindigkeit$  $\omega$ (1/s)', size='large')
  ax2.grid()

  #ax3=fig.add_subplot(3,1,3)
  #ax3.plot(freq, amp, 'b.')
  #ax3.set_xlabel('$Frequenz$ $f$ (Hz)', size='large')
  #ax3.set_ylabel('$Amplitude$', size='large')
  #ax3.set_yscale('log')
  #ax3.grid()

  ax4=fig.add_subplot(3,1,3)
  x=np.linspace(0., 220, 100)
  ax4.plot(teins, Ekin(omega,T=1.39e-03), 'b-', alpha=0.5, linewidth=1)
  ax4.set_xlabel('$Zeit$ (s)', size='large')
  ax4.set_ylabel('$Energie$ (J)', size='large')
  
  plt.savefig('Aufgabe1a.pdf')

# Phasenraumdiagramme mit interpolierten Daten
  fig2=plt.figure("PhaseSpace", figsize=(10, 5.))
  fig2.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.93,
                    wspace=.25, hspace=None)
  axa=fig2.add_subplot(1,2,1)
  tplt=np.linspace(t[0], t[-1], 50000)
  axa.plot(cs_phi(tplt), cs_omega(tplt), 'b-', alpha=0.5, linewidth=1)
  axa.set_xlabel('$\phi$', size='large')
  axa.set_ylabel('$\omega$ (1/s)', size='large')
  axa.set_title("Phasenraum-Diagramm")
  
  axb=fig2.add_subplot(1,2,2)
  tplt=np.linspace(t[int(0.75*len(t))], t[-1], 20000)
  axb.plot(cs_phi(tplt), cs_omega(tplt), 'b-', alpha=0.5, linewidth=1)
  axb.set_xlabel('$\phi$', size='large')
  axb.set_ylabel('$\omega$ (1/s)', size='large')
  axb.set_title("Phasenraum-Diagramm (Ausklingphase)")  
  plt.show()
  #plt.savefig('Aufgabe1b.pdf')
  
  
  #besti