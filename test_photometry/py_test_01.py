from astropy.io import fits
fits_image_filename = fits.open('C:/Users/Jake/MPhys/test_photometry/PAL5_3p6um.fits')

hdul = fits.open(fits_image_filename)
hdul.info()