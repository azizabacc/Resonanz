# script test_readColumnData.py
# -*- coding: utf-8 -*-
''' test data input from text file with module PhyPraKit.readColumnData
..  author:: Guenter Quast <g.quast@kit.edu>
'''

from PhyPraKit import readColumnData,odFit
import numpy as np, matplotlib.pyplot as plt

fname='Aufg1.dat'
ncols=2

data_array, info_dict =\
  readColumnData(fname, ncols, delimiter=' ', pr=False)

# print what we got:
x=data_array[:,0] # 1st column
y=data_array[:,1] # 2nd column
omega=data_array[:,2] # 3rd column

print "Title=", info_dict['*TITLE']
print "x= ", x
print "y= ", y

    
print "keywords found:"
for key in info_dict:
  if (info_dict[key]!=None): print key,':', info_dict[key]



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
  #p0p=[2., 100., 0.001] # initial parameters 
  envd=env_linexp_d
  #p0d=[-2., 100., 0.001]  




# plot data and analysis results
  fig=plt.figure("Amplitude", figsize=(7.5, 10.))
  fig.suptitle('Script: Beispiel_Drehpendel.py', size='x-large', color='b')
  fig.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.93,
                    wspace=None, hspace=.25)#

  ax1=fig.add_subplot(3,1,1)
  ax1.plot(x,y)
  
  
 
 #a=np.linspace(0.0001,0.00)
  ''' 
  for i in xrange(1,4):
       ax1.plot(x, envp(x, 2.8,150-i*25,0.0004), 
         'r-', linewidth=1, label=u'positive Einhüllende')
       ax1.plot(x, envd(x, 2.8,150-i*25,0.0005), 
         'g-', linewidth=1, label=u'negative Einhüllende')
  ''' 
  ax1.plot(x, envp(x, 2.8,70,0.0006), 
    'r-', linewidth=1, label=u'positive Einhüllende')
  ax1.plot(x, envd(x, 2.8,100,0.0004), 
    'g-', linewidth=1, label=u'negative Einhüllende')
  ax1.set_xlabel('$Zeit$ $t$ (s)', size='large')
  ax1.set_ylabel('$Winkel$  $\phi$ (rad)', size='large')
  #ax1.legend(loc='best', numpoints=1, prop={'size':10})
  ax1.grid()
'''
  ax2=fig.add_subplot(3,1,2)
  ax2.plot(x , omega)
  ax2.set_xlabel('$Zeit$ $t$ (s)', size='large')
  ax2.set_ylabel('$Winkelgeschwindigkeit$  $\omega$ (1/s)', size='large')

  ax3=fig.add_subplot(3,1,3)
  ax3.plot(y, omega, 'b.')
  ax3.set_ylabel('$Winkel$  $\phi$', size='large')
  ax3.set_xlabel('$Winkelgeschwindigkeit$  $\omega$ (1/s)', size='large')
  ax3.set_yscale('log')
  ax3.grid()

'''
'''
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
'''
'''
# make a plot
fig=plt.figure(1, figsize=(5.,5.))
plt.plot(x , y)
# show systematic errors
#exs = x * np.float32(info_dict['*xRelCor'])
#eys = y * np.float32(info_dict['*yRelCor'])
#plt.errorbar(x , y, xerr=exs, yerr=eys, fmt='r.')
# and finally, the axis labels


plt.xlabel(info_dict['*xLabel'] +' / ' + info_dict['*xUnit'], size='x-large')
plt.ylabel(info_dict['*yLabel'] +' / ' + info_dict['*yUnit'], size='x-large')
'''
plt.show()

