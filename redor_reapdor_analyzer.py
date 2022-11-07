import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import curve_fit


# load bruker processed file and convert it to NMRPipe
dic, data = ng.bruker.read_pdata('C:/Users/igordas/OneDrive/Dipolar_Techniques/20221027_setup_redor/5/pdata/1/')

u = ng.bruker.guess_udic(dic, data)

# separate S0 and S and delete lines with 0's
S0 = data[1::2]
S0 = S0[~np.all(S0 == 0, axis=1)]

S = data[::2]
S = S[~np.all(S == 0, axis=1)]

#-------------------------------------------------
# Check spectrum for region to get maximum - comment it after selection!!
#fig = plt.figure(dpi=300) 
#ax = fig.add_subplot(111)
#ax.plot(S0[0])
#ax.set_xlim(1750, 2250)
#--------------------------------------------------------------

# calculate deltaS based on max values
delta = np.zeros(len(S))

for i in range(len(delta)):
    delta[i] = (max(S0[i][1750:2250])-max(S[i][1750:2250]))/max(S0[i][1750:2250])
    
# creates the NTr axis
vclist = np.linspace(1,16,16)
mas = 15000 # MAS frequency /Hz
period = 1/mas
ntr = period*vclist

#parabolic fitting
def first_order(x,a):
    y = (4/(3*(np.pi**2)))*(x**2)*a
    return y

def second_order(x,b):
    y = (16/15)*(x**2)*(b**2) - (128/315)*(x**4)*(b**4)
    return y

par1, _ = curve_fit(first_order, ntr[0:3], delta[0:3])

m = par1

par2, _ = curve_fit(second_order, ntr[0:7], delta[0:7])

D = par2

x1 = np.linspace(0,0.0003)
fit_1st = np.zeros(len(x1))

for j in range(len(fit_1st)):
    fit_1st[j] = first_order(x1[j],m)
    
x2 = np.linspace(0, 0.0008)
fit_2nd = np.zeros(len(x2))

for k in range(len(fit_2nd)):
    fit_2nd[k] = second_order(x2[k], D)

    
#calculation of distance for 1st order approximation
gamma_I = 27.116e6 #gyromagnetic ratio for 15N taken from wikipedia
gamma_S = 67.2828e6 #gyromagnetic ratio for 13C taken from wikipedia
S = 0.5 #13C spin

r6 = (4/15)*((1e-7)**2)*S*(S+1)*(gamma_I**2)*(gamma_S**2)*(sc.hbar**2)/m

r_1st = (r6**(1/6))*1e10 #distance in angstrom

dist_1 = round(r_1st[0],2)

#calculation of distance for 2nd order approximation
r3 = (1e-7)*(gamma_I)*(gamma_S)*(sc.hbar)/(2*np.pi*D)

r_2nd = (r3**(1/3))*1e10

dist_2 = round(r_2nd[0],2)

# plot the spectrum
fig = plt.figure(dpi=300) 
ax = fig.add_subplot(111)
ax.plot(ntr*1000,delta,'k.',label='Experimental Data')
ax.plot(x1*1000,fit_1st,'r-',label='1st order appr.')
ax.plot(x2*1000,fit_2nd,'b-',label='2nd order appr.')
# decorate axes
# ax.set_title("Protein 1D Spectrum")
ax.set_xlabel("NT$_r$ /ms", fontsize=16)
ax.set_ylabel('$\Delta S/S_0$', fontsize=16)
# ax.set_yticks([])
ax.set_xlim(0, 1.2) #set the limit for x-axis if needed
ax.set_ylim(0,1)
ax.text(0.2,0.5,'r = {} $\AA$'.format(dist_1))
ax.text(0.73,0.5,'r = {} $\AA$'.format(dist_2))
ax.set_title('$^{13}$C{$^{15}$N} REDOR - Glycine',fontsize = 18)
ax.legend()