import numpy as np
import matplotlib.pyplot as plt
import sys
import pygedm

nargs = len(sys.argv) - 1
if (nargs==3):
    model = sys.argv[1]
    if model!="ne2001" and model!="ymw16":
        print("Invalid ISM model")
        print("Valid models are ne2001 and ymw16")
        sys.exit()
    h_z = np.float(sys.argv[2])
    filename = sys.argv[3]
else:
    print("Needs an ISM model, a scaleheight, and an output filename, e.g.:")
    print("python galactic_frbs.py ne2001 0.05 frbsample.csv")
    print("for h_z = 50 pc (0.05 kpc)")
    print("models are ne2001 and ymw16")
    sys.exit()

# generete this many FRBs -- not all of them survive the masking though!
# about 30,000 survive the exponential disk mask
Nfrbs = 10000000 // 4

# parameters of the simulation
rgcmax = 30  # kpc, maximum distance of an FRB from Galactic center

# this creates an exponentially declining density distribution
# of FRBs in the Galactic disk, scale length 4 kpc, 2-D distribution at this stage
rgc = (np.random.random(Nfrbs) ** 0.5) * rgcmax
test = np.random.random(Nfrbs)
mask = test < np.exp(-rgc / 4.0)  # 4 kpc is the scale-length of the exponential

# create the exponential distribution 
rgc = rgc[mask]

# count how many FRBs were made
Nfrbs = len(rgc)

# some Galactic parameters
Rsun = 8.5  # kpc
Rmax = 30.0  # kpc

# if we are using the NE2001 model, then just generate x,y positions
# for the FRBs made above in the old stellar disk. If using the YMW16
# model, start from scratch but for the exact same number of FRBs as
# the exp model, but generate new positions based on the electron
# density n_e map for YMW16 at Z=0. This mocks up a young disk
# population instead.

if model == "ne2001":
    print("Using exponential disk model to generate FRB positions")
    # now make the x,y,z values for the events generated above
    phi = np.random.random(len(rgc)) * 2.0 * np.pi
    x_frb = rgc * np.cos(phi)
    y_frb = rgc * np.sin(phi)

if model == "ymw16":
    print("Using YMW16 n_e model at the mid-plane to generate FRB positions")
    # load a 100x100 grid of n_e densities at Z=0 in the YMW16 model
    d_2016 = np.load("YMW2016.nemodel.npy")
    # central pixel in map set to zero (center of Galaxy gets no FRBs)
    d_2016[50, 50] = 0.0
    # normalise map to a maximum value of 1
    d_2016_norm = d_2016 / np.max(d_2016)
    x_frb = []
    y_frb = []
    frb_count = 0
    # generate mock FRBs similarly to electron density map
    while frb_count < Nfrbs:
        # generate a random x,y position in the galaxy
        x1 = np.random.rand() * 100
        y1 = np.random.rand() * 100
        ix = int(x1)
        iy = int(100 - y1)
        # if this is not in the exact middle of the galaxy
        if not (ix == 50 and iy == 50):
            # generate a random number in range [0,1]
            rannum = np.random.rand()
            # use the normalised map and this number to decide if
            # an FRB should be at this position
            if d_2016_norm[iy, ix] > rannum:
                x_frb = np.append(x_frb, (x1 - 50) * 0.4)
                y_frb = np.append(y_frb, (y1 - 50) * 0.4)
                frb_count += 1
                
# the total number of events
Nfrb = len(x_frb)
print("Total number of simulated FRBs = "+str(Nfrb))

# distribute as a Gaussian in Z with h_z scaleheight
z_frb = np.random.normal(0.0, h_z, size=len(x_frb))

# compute galactic coordinates, gl and gb
r2d = 180.0/np.pi
gl = np.arctan2(x_frb-Rsun,y_frb)
gl = gl + np.pi/2.0
for i in range(len(gl)):
    if gl[i] > np.pi:
        gl[i] = gl[i] - 2.0*np.pi
gl = gl*r2d
dist = np.sqrt((x_frb-Rsun)**2+y_frb**2+z_frb**2)
gb = np.arctan2(z_frb,dist)*r2d

# work out the DM for each simulated FRB
frb_dm = np.zeros(len(x_frb))
for i in range(len(x_frb)):
    dm, tau = pygedm.dist_to_dm(float(gl[i]), float(gb[i]), float(dist[i])*1000.0, method=model)
    frb_dm[i] = dm.value

f0 = 1.4 # GHz
# figure out the scattering of each FRB based on its DM and the Bhat data
# log (base10) of the scattering time, using the Bhat relation (Bhat 2006)
logtscatt = (-6.46 + 0.154*np.log10(frb_dm) + 1.1*np.log10(frb_dm)**2 - 3.9*np.log10(f0))

# convert log of the scattering time to the scattering time
t_scatt = 10.0**logtscatt

# intrinsic width of FRBs -- we need to assume some sort of lower limit
t_intrinsic = 0.61 / 1e5
t_dm = 8.3*(122e-3)*frb_dm/f0**3

#convert this to ms
t_dm_ms = t_dm*0.001
t_pulse = np.sqrt(t_intrinsic**2 + t_scatt**2 + t_dm_ms**2)

# we also add a factor of 10 scatter around this DM-scatter
# relation, since that's what the pulsars show
t_pulse = np.log10(t_pulse) + np.random.normal(0, 1.0, size=len(x_frb))
t_pulse = 10.0**t_pulse

f = open(filename,'w')
f.write("#x,y,z,dist,DM,width\n")
for i in range(len(x_frb)):
    if (t_pulse[i]<1e10):
        line = str(x_frb[i])+","+str(y_frb[i])+","+str(z_frb[i])+","
        line += str(dist[i])+","+str(frb_dm[i])+","+str(t_pulse[i])
        f.write(line+"\n")
f.close()
print("Finished. Output in "+filename)

sys.exit()
    
