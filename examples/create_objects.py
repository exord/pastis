from matplotlib import pyplot as plt
plt.ion()
import pastis

# Initialise pastis(takes some time). This reads the BT-Settl spectra and
# Darmouth evolution tracks, as well as some photometric filters.
pastis.initialize()

# Import the AstroClasses module (must be done once pastis is initialised)
from pastis import AstroClasses as ac

# Instaniate a star whose Teff, logg and z are known.
# The star is assumed to be at 10 pc, with zero redenning.
star = ac.Target(teff=5700, z=0, logg=4.4, ebmv=0.0, dist=10)

# Once the star is instantaniated, you can get its magnitude, for example
# we iterate over all bands loaded at initialization (there are many many more!)
for band in pastis.photometry.NormalisedFilters:
    # This is linked with the CoRoT filters, which do not have a fixed
    # zero point. We skip them.
    if band not in pastis.photometry.ZeroMag:
        continue
    print('Magnitude in {} band: {:.2f}.'.format(band, star.get_mag(band)))

# We can plot the spectrum
plt.semilogx(pastis.photometry.ww, star.get_spectrum(), '-b',
             label='Teff=5700; Dist=10; E(B-V)=0')
plt.xlim(0, 1e4)

# Now let's see how changing some parameters affect the spectrum
# First make the stellar Teff smaller
star.teff = 4500

# Plot the new spectrum (the force argument is new! download the new version if you see the same spectra twice)
plt.semilogx(pastis.photometry.ww, star.get_spectrum(force=True), '-r',
             label='Teff=4500; Dist=10; E(B-V)=0')

# Now let's include some redenning on the first spectrum
star.teff = 5700
star.ebmv = 0.05

plt.semilogx(pastis.photometry.ww, star.get_spectrum(force=True), '-g',
             label='Teff=5700; Dist=10; E(B-V)=0.05')



