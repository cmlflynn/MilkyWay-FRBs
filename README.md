# MilkyWay-FRBs


There are two codes: 

galactic_frbs.py -- which creates a database of FRB sources in X,Y,Z around the Milky Way, with DMs from an ISM model, and a width based on the Bhat curve (relates DM to pulse width / scattering)

analyze_frbs.py -- which plots the results of such a database. 

To run:

python galactic_frbs.py ymw16 0.05 1.4 frbsample_ymw16_0.05_1.4.csv

This uses YMW16 as the ISM mode, a scaleheight of 0.05 kpc in Z, a survey frequency of 1.4 GHz (for the Bhat law), and writes the results to the CSV file. 

Example output:

#x,y,z,dist,DM,width
3.5220033845868617,1.0918745087119022,0.04382651599596037,5.096524404861009,634.9056396484375,273.9694149697946
-5.745411513317328,8.956414266633045,-0.0342995714854575,16.827069922000856,2351.2783203125,5348.314650032203
-8.29202207043251,9.698977288811717,0.09359871688225069,19.392032549028478,2277.095458984375,277399.1958871533

X,Y,Z coordinates are in kpc, dist is the distance from the Sun (kpc), DM is in pc/cc (derived from ISM model) and the width is derived from the Bhat relation (ms)

Some widths can be very large -- no width larger than 100 ms has significance as it is essentially undetectible in any present or planned survey.

To analyse the FRBs:

python analyze_frbs.py frbsample_ymw16_0.05_1.4_.csv 1e25 -1.5

This reads in the sources in the CSV file, assigns a luminosity to each from a power law LF with E0=1E25 erg and powerlaw-slope -1.5.

All sources then have a S/N computed at Earth, and S/N>10 implies a detection. 

It will produce a plot which is a comparison of the FRBs detectible in STARE2 and GReX, and plots of the various observational and other parameters of each FRB as seen by STARE2. 




