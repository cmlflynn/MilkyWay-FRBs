import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import sys

def power_law(N,E0,alpha):
    # create simulated counts N in a power law distribution
    # random numbers in [0,1] for use in generating energy distribution
    x = np.random.rand(N)
    # generate power law from the random numbers x
    # transform random numbers to power law if the
    # requested law isn't flat (to within 0.001)
    if np.abs(alpha)>0.001:
        x = x**(1/alpha)
    else:
        print("Power-law is flat")
        sys.exit()
    # set up energies relative to minimum energy E0
    E = E0*x
#    E = E0*x/x
    return E

def dN_dE_plot(E,bins):    
    # create a log-log histogram of the energies E
    y, x = np.histogram(np.log10(E),bins=bins)
    # a small number to add to the FRB counts so that the
    # logarithm function can handle things
    epsilon = 1e-5
    y = np.log10(y+epsilon)
    # create bin centers from bin edges array
    dx = x[1]-x[0]
    x = x[:-1]+dx/2.0
    return x, y

def get_SN(E,w,telescope):
    print(telescope)

    if telescope == 'STARE2':
        # return S/N for STARE2 telescope from FRB energy and width
        k = 1.38e-23 # Boltzmann constant in m^2 kg s^-2 K^-1
        f = 0.36 # aperture efficiency, chosen to yield measured SEFD
        area = 0.2*0.2  # m^2
        gain = area*f/2.0/k*1E-26
        Tsys = 100.0 # K
        SEFD = Tsys/gain
        w_sec = w*1E-3 # width in sec
        f0 = 1.4E9 # Hz
        E_detected = E/(4.0*np.pi*(dist*3.09e+19)**2) # J/m^2
        power_detected = E_detected/w_sec
        flux_density = 1e+26*power_detected / f0 # Jy?
        BW = 188E6 # Hz
        Np = 2 # polarisations
        SN = flux_density*(1.0/SEFD)*np.sqrt(BW*Np*w_sec)

    if telescope == 'GReX':
        # return S/N for GReX telescope from FRB energy and width
        k = 1.38e-23 # Boltzmann constant in m^2 kg s^-2 K^-1
        f = 0.36 # aperture efficiency, chosen to yield measured SEFD
        area = 0.2*0.2  # m^2
        gain = area*f/2.0/k*1E-26
        # STARE2 --> GReX
        # gain is 1/3rd of STARE2
        gain = gain / 3.0
        Tsys = 25.0 # K
        SEFD = Tsys/gain
        w_sec = w*1E-3 # width in sec
        f0 = 1.4E9 # Hz
        E_detected = E/(4.0*np.pi*(dist*3.09e+19)**2) # J/m^2
        power_detected = E_detected/w_sec
        flux_density = 1e+26*power_detected / f0 # Jy?
        BW = 188E6 # Hz
        Np = 2 # polarisations
        SN = flux_density*(1.0/SEFD)*np.sqrt(BW*Np*w_sec)

    if telescope != "STARE2" and telescope != "GReX":
        print("Invalid telescope :"+telescope)
        sys.exit()

    print("SEFD = "+str(np.around(SEFD/1e6,2))+" MJy")
        
    return SN

############################################################################

nargs = len(sys.argv) - 1
if (nargs==3):
    filename = sys.argv[1]
    E0 = np.float(sys.argv[2])
    alpha = np.float(sys.argv[3])
else:
    print("Needs an archive file name, E0, and alpha e.g.:")
    print("python analyze_frbs.py frbsample.csv 1.0E23 -0.25")
    sys.exit()

# read the mock FRB catalog 
x,y,z,dist,dm,w = np.loadtxt(filename,usecols=(0,1,2,3,4,5),
                  dtype='float',delimiter=',',unpack=True)

# count how many FRBs were read in 
Nfrbs = len(x)

# assign energies for each FRB
E = power_law(Nfrbs,E0,alpha)

###########################################################################

font = {'family' : 'normal',
#                'weight' : 'bold',
                'size'   : 16}

plt.rc('font', **font)

# maximum search width for FRBs
w_max = 100 # ms

# mask for width limit of FRB search
wmask = w<w_max

# get S/N for each FRB from energies and widths for STARE2
SN = get_SN(E,w,telescope='STARE2')
SN_narrow = SN[wmask]
dist_narrow = dist[wmask]
snmask = SN_narrow > 10
bins = np.arange(0.0,30.0,1)
# scatter plot of FRB distances versus distance from Sun for STARE2
stare2_dist = dist_narrow[snmask]
sky = np.random.rand(len(stare2_dist))
# simulate whether FRB is in the STARE2 sky
# Galactic plane is on the sky for about 8 hours out of every 24 hours, 1/3rd of the time
mask = sky<1.0/3.0
stare2_dist = stare2_dist[mask]
# simulate whether FRB is visible (in Northern sky)
# there is roughly a 30/70 split for STARE2 between Northern and Southern sources
north_south = np.random.rand(len(stare2_dist))
mask = north_south < 0.3
stare2_dist = stare2_dist[mask]
plt.hist(stare2_dist,bins=bins,color='b',label="STARE2",width=0.9,alpha=0.5)

# now for GReX
SN_grex = get_SN(E,w,telescope='GReX')
SN_narrow_grex = SN_grex[wmask]
snmask_grex = SN_narrow_grex > 10
dist_narrow = dist[wmask]
grex_dist = dist_narrow[snmask_grex]
sky = np.random.rand(len(grex_dist))
mask = sky < 0.5
grex_dist = grex_dist[mask]
# histogram of FRB distances versus distance from Sun for GReX
plt.hist(grex_dist,bins=bins,color='g',label="GReX",width=0.9,alpha=0.5)

plt.xlabel("dist [kpc]")
plt.ylabel("N")
plt.title("STARE2 versus GReX")
plt.legend()

print("GReX FRBs : "+str(np.sum(grex_dist)))
print("STARE2 FRBs : "+str(np.sum(stare2_dist)))

print("ratio = "+str(np.sum(grex_dist)/np.sum(stare2_dist)))

plt.show()


# set up plot size
plt.figure(figsize=(20,12))

plt.rc('font', **font)

# number of subplots in x and y
nx = 4
ny = 2

# some other stuff for plotting
Rmax = 30 # kpc
Rsun = 8.5 # kpc
Zmax = 5.0 # kpc

# maximum search width for FRBs
w_max = 100 # ms

# mask for width limit of FRB search
wmask = w<w_max

# save data for FRBs that are narrow enough to be detectible
# we'll apply the S/N cut later
x_narrow = x[wmask]
y_narrow = y[wmask]
z_narrow = z[wmask]
SN_narrow = SN[wmask]
dist_narrow = dist[wmask]
dm_narrow = dm[wmask]
w_narrow = w[wmask]
E_narrow = E[wmask]

# S/N mask, to be applied after the pulse width mask
snmask = SN_narrow > 10

# power law energy distribution for FRBs
nsubplot = 1
ax = plt.subplot(ny,nx,nsubplot)
dE = 0.25 # Joules
bins = np.arange(np.log10(E0)-dE,40.0,dE)
# bin all the FRBs 
Ebin, N_Ebin = dN_dE_plot(E,bins)
plt.bar(Ebin,N_Ebin,color='b',width=0.4,label="All FRBs")
# bin the FRBs that are narrow enough to detect
Ebin, N_Ebin = dN_dE_plot(E[wmask],bins)
plt.bar(Ebin,N_Ebin,color='g',width=0.4,label="w<"+str(w_max)+"ms")
# bin the narrow FRBs that have S/N>10
Ebin, N_Ebin = dN_dE_plot(E_narrow[snmask],bins)
plt.bar(Ebin,N_Ebin,color='r',width=0.4,label="S/N>10")
plt.ylim(0,)
plt.xlabel("log10(E) [J]")
plt.ylabel("log10(N)")
plt.title(filename+"\nLF model : E0 = "+str(E0)+" J, alpha = "+str(alpha))
plt.text(0.9, 0.9, '(a)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#plt.legend()
nsubplot += 1

# scatter plot of FRB widths versus distance from Sun
ax = plt.subplot(ny,nx,nsubplot)
#plt.loglog(dist,w,'b,',label="All FRBs")
#plt.loglog(dist[wmask],w[wmask],'g,',label="w<"+str(w_max)+"ms")
#plt.loglog(dist_narrow[snmask],w_narrow[snmask],'r.',ms=1,label="S/N>10")
plt.plot(np.log10(dist),np.log10(w),'b,',label="All FRBs")
plt.plot(np.log10(dist[wmask]),np.log10(w[wmask]),'g,',label="w<"+str(w_max)+"ms")
plt.plot(np.log10(dist_narrow[snmask]),np.log10(w_narrow[snmask]),'r.',ms=1,label="S/N>10")
plt.xlabel("dist [kpc]")
plt.ylabel("width [ms]")
plt.text(0.1, 0.9, '(b)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#plt.legend()
nsubplot += 1

# scatter plot of FRB widths versus their DMs
ax = plt.subplot(ny,nx,nsubplot)
#plt.loglog(dm,w,'b,',label="All FRBs")
#plt.loglog(dm[wmask],w[wmask],'g,',label="w<"+str(w_max)+"ms")
#plt.loglog(dm_narrow[snmask],w_narrow[snmask],'r.',ms=1,label="S/N>10")
plt.plot(np.log10(dm),np.log10(w),'b,',label="All FRBs")
plt.plot(np.log10(dm[wmask]),np.log10(w[wmask]),'g,',label="w<"+str(w_max)+"ms")
plt.plot(np.log10(dm_narrow[snmask]),np.log10(w_narrow[snmask]),'r.',ms=1,label="S/N>10")
plt.xlabel("DM [pc/cc]")
plt.ylabel("width [ms]")
plt.text(0.1, 0.9, '(c)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#plt.legend()
nsubplot += 1

# view of FRBs in Milky Way plane, ie from above
ax = plt.subplot(ny,nx,nsubplot)
plt.plot(x,y,'b,',label="FRBs")
plt.plot(x[wmask],y[wmask],'g,',label="w<"+str(w_max)+"ms")
plt.plot(x_narrow[snmask],y_narrow[snmask],'r.',ms=1,label="S/N>10")
plt.xlim(-Rmax,Rmax)
plt.ylim(-Rmax,Rmax)
plt.plot(Rsun, 0.0, 'yo',label ="Sun")
plt.plot(0.0, 0.0, 'ko', label="Galactic Center")
plt.xlabel("X [kpc]")
plt.ylabel("Y [kpc]")
plt.text(0.9, 0.9, '(d)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#plt.legend()
nsubplot += 1

uselog = False
uselog = True

# histograms of FRB distances by selected types
ax = plt.subplot(ny,nx,nsubplot)
bins = np.arange(0.0,30.0,0.25)
plt.hist(dist,bins=bins,color='b',label="All FRBs",log=uselog)
plt.hist(dist[wmask],bins=bins,color='g',label="w<"+str(w_max)+"ms",log=uselog)
plt.hist(dist_narrow[snmask],bins=bins,color='r',label="S/N>10",log=uselog)
plt.xlabel("dist [kpc]")
plt.ylabel("N")
plt.text(0.9, 0.9, '(e)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#plt.legend()
nsubplot += 1

# histograms of FRB DMs by selected types
ax = plt.subplot(ny,nx,nsubplot)
bins = np.arange(0.0,3000.0,50.0)
plt.hist(dm,bins=bins,color='b',label="All FRBs",log=uselog)
plt.hist(dm[wmask],bins=bins,color='g',label="w<"+str(w_max)+"ms",log=uselog)
plt.hist(dm_narrow[snmask],bins=bins,color='r',label="S/N>10",log=uselog)
plt.xlabel("DM [pc/cc]")
plt.ylabel("N")
plt.text(0.9, 0.9, '(f)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#plt.legend()
nsubplot += 1

# histograms FRB S/Ns by selected types
ax = plt.subplot(ny,nx,nsubplot)
snbins = np.arange(-2,4,0.1)
plt.hist(np.log10(SN),bins=snbins,log=True,color='b',label="All FRBs")
plt.hist(np.log10(SN[wmask]),bins=snbins,log=True,color='g',label="w<"+str(w_max)+"ms")
plt.hist(np.log10(SN_narrow[snmask]),bins=snbins,log=True,color='r',label="S/N>10")
plt.xlabel("log(S/N)")
plt.ylabel("N")
plt.text(0.9, 0.9, '(g)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
nsubplot += 1

# side view of FRBs in Milky Way
ax = plt.subplot(ny,nx,nsubplot)
plt.plot(x,z,'b,',label="FRBs")
plt.plot(x[wmask],z[wmask],'g,',label="w<"+str(w_max)+"ms")
plt.plot(x_narrow[snmask],z_narrow[snmask],'r.',ms=1,label="S/N>10")
plt.xlim(-Rmax,Rmax)
plt.ylim(-Zmax,Zmax)
plt.plot(Rsun, 0.0, 'yo',label ="Sun")
plt.plot(0.0, 0.0, 'ko', label="Galactic Center")
plt.xlabel("X [kpc]")
plt.ylabel("Z [kpc]")
plt.text(0.9, 0.9, '(h)', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#plt.legend()

plt.tight_layout()

plt.show()

plt.loglog(dm_narrow[snmask],E_narrow[snmask],'b.',label="Mock FRBs")
plt.xlabel("DM [pc/cc]")
plt.ylabel("E [J]")
plt.loglog(333.0,1e27,'y*',label="SGR1935",ms=20)
plt.title(filename+"\nLF model : E0 = "+str(E0)+" J, alpha = "+str(alpha))
plt.show()
