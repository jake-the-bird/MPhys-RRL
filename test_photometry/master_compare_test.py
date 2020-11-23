import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import glob
import time
import copy
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.nddata import NDData
from astropy.visualization import simple_norm
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.time import Time
from photutils import DAOStarFinder
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from photutils.psf import IterativelySubtractedPSFPhotometry as ISPSF
from photutils.psf import extract_stars
from photutils.psf.groupstars import DAOGroup
from photutils import EPSFBuilder
from photutils.background import MMMBackground
from photutils.utils import make_random_cmap


# could add a check to see if apertures need to be plotted (in case i want to plot just the fits)
def plotting(region, ap, an, cmap, choice, savename):
    plt.imshow(region, cmap=cmap, origin='lower', norm=LogNorm(), interpolation='nearest')
    
    if choice == 'ap_only':
        ap.plot(color='blue', lw=.5, alpha=1.)
    elif choice == 'both_aps':
        ap.plot(color='blue', lw=.5, alpha=.5)
        an.plot(color='red', lw=.5, alpha=.5)
        
    plt.colorbar(fraction = 0.05)
    #plt.grid(True)
    plt.minorticks_on()
    plt.grid(b=True, which='major', lw=.5, alpha=.4, color='black')
    plt.grid(b=True, which='minor', lw=.5, alpha=.2, color='black')
    
    if savename != None:
        plt.savefig('images/'+savename, dpi=500)
    
    plt.gcf().set_size_inches(10, 6)
    plt.show()
    #plt.close()
    
def print_table_nicely(table):   # note: 'table' must be an astropy Table type
    for col in table.colnames:
        table[col].info.format = '%.8g'
    print(table)
    
def star_find(data, sigma, fwhm, std, roundness, sharphi):
    daofind = DAOStarFinder(threshold=sigma*std, fwhm=fwhm, roundlo=-roundness, roundhi=roundness, sharphi=sharphi)
    sources = daofind(data)
    return sources
    
def ap_phot(data, ap, an, bkg_method):
    all_apers = [ap, an]
    table = aperture_photometry(data, all_apers)
    
    if bkg_method == 'mean':
        bkg_mean = table['aperture_sum_1'] / an.area
        table['bkg_sum_mean'] = bkg_mean * ap.area
        table['star_bkgsub'] = table['aperture_sum_0'] - table['bkg_sum_mean']
        
        an_err = bkg_mean * ap.area / np.sqrt(an.area)
        ap_err = bkg_mean * np.sqrt(ap.area)
        poisson_err = np.sqrt(table['star_bkgsub'])
        table['counts_err'] = np.sqrt(an_err**2 + ap_err**2 + poisson_err**2)
        #table['perc_err'] = table['counts_err'] / table['star_bkgsub']
    
    elif bkg_method == 'median':
        annulus_masks = an.to_mask(method='center')
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_mask, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_mask)
            
        bkg_median = np.array(bkg_median)
        table['bkg_sum_median'] = bkg_median * ap.area
        table['star_bkgsub'] = table['aperture_sum_0'] - table['bkg_sum_median']
        
        an_err = bkg_median * ap.area / np.sqrt(an.area)
        ap_err = bkg_median * np.sqrt(ap.area)
        poisson_err = np.sqrt(table['star_bkgsub'])
        table['counts_err'] = np.sqrt(an_err**2 + ap_err**2 + poisson_err**2)
        #table['perc_err'] = table['counts_err'] / table['star_bkgsub']
        
    else:
        print('\n\nPlease give a valid bkg_method kthx\n\n')
            
    return table

def apparent_magnitude_err(counts, apcorr, zmag_err, counts_err, apcorr_err):
    mag_err = np.sqrt(zmag_err**2 + (2.5 * np.sqrt((counts_err/counts)**2 + (apcorr_err/apcorr)**2) / np.log(10))**2)
    return mag_err

def make_catalog(x, y, header):
    crd = np.transpose((x, y))
    w = wcs.WCS(header)
    world = w.wcs_pix2world(crd, 0)
    ra = world[:,0]
    dec = world[:,1]
    cat = SkyCoord(ra, dec, frame='icrs', unit='deg')
    
    return cat, ra, dec


run_all = False  # set to False to run only the first epoch
base_dir = 'C:/Users/Jake/MPhys-code/MPhys-RRL/test_photometry/data/PAL5/'
channel = '3p6um'
#channel = '4p5um'
method = 'aperture'
#method = 'PSF'

star_matching = 'no'

sigma_level = 6.
sigma_level_PSF = 50.
FWHM = 5.
r_ap = 6.
r_in = 6.
r_out = 14.
roundness = 0.5
sharphi = 0.9
sharphi_PSF = 0.7
#sharplo = 0.65

if channel == '3p6um':
    zmag       = 18.80    # given in/calculated using IRAC handbook, section 4.8
    zmag_err   =  0.02    # calculated from zmag = 2.5log10(F0/C), F0 = 280.9 +/- 4.1 (from IRAC, table 4.1)
    apcorr     =  1.1233  # aperture correction for 6, 6-14 pix apertures in channel 1; given on IRAC website
    apcorr_err =  0.0225  # taking accuracy to be ~2%, as mentioned in IRAC, section 4.10
elif channel == '4p5um':
    zmag       = 18.32
    zmag_err   =  0.02
    apcorr     =  1.1336
    apcorr_err =  0.0227
else:
    print('Please select a valid channel')


file_m = base_dir+'PAL5/PAL5_'+channel+'.fits'
with fits.open(file_m) as hdu_list:        
    print(file_m)
    image_m = hdu_list[0].data
    hdr_m = hdu_list[0].header
    exptime = hdr_m['EXPTIME']
    fluxconv = hdr_m['FLUXCONV']
    conv = exptime / fluxconv
    data_m = image_m * conv

mean_m, median_m, std_m = sigma_clipped_stats(data_m, sigma=sigma_level)

sources_m = star_find(data_m, sigma_level, FWHM, std_m, roundness, sharphi)
pos_m = np.transpose((sources_m['xcentroid'], sources_m['ycentroid']))
apertures_m = CircularAperture(pos_m, r_ap)
annuli_m = CircularAnnulus(pos_m, r_in, r_out)
print('No. of stars detected: {}'.format(len(sources_m)))
plotting(data_m, apertures_m, annuli_m, 'Greys', choice='both_aps', savename=None)

sources_m['id'].name = 'id_master'

cat_m, ra_m, dec_m = make_catalog(sources_m['xcentroid'], sources_m['ycentroid'], hdr_m)
sources_m['RA'] = ra_m
sources_m['dec'] = dec_m


print('Photometry method selected: {}'.format(method))

###   INITIALISING IMAGE WE'LL BE COMPARING OTHER EPOCHS TO WHEN STAR MATCHING   ###

comp_epoch = 1  # epoch to compare all others against when matching stars   
epoch = 1  # counter to be incremented every epoch
LC_mag = []
LC_err = []
LC_time = []

for filename in glob.glob(base_dir+'*/PAL5__e[0-9]_'+channel+'.fits', recursive=True) + glob.glob(base_dir+'*/PAL5__e[0-9][0-9]_'+channel+'.fits', recursive=True):
    
    ###   OPENING FITS FILE AND CONVERTING TO COUNTS   ###
    
    with fits.open(filename) as hdu_list:        
        print(filename)
        image_data = hdu_list[0].data
        hdr = hdu_list[0].header
        exptime = hdr['EXPTIME']
        fluxconv = hdr['FLUXCONV']
        conv = exptime / fluxconv
        print('EXPTIME: {0}\nFLUXCONV: {1}'.format(exptime, fluxconv))
        data = image_data * conv
    
    file_corr = base_dir+'PAL5__e'+str(epoch)+'/PAL5__e'+str(epoch)+'_correction_'+channel+'.fits'
    with fits.open(file_corr) as hdu_list:
        data_corr = hdu_list[0].data
        
    mean, median, std = sigma_clipped_stats(data, sigma=sigma_level)
    print('Mean: {0}\nMedian: {1}\nStd dev: {2}'.format(mean, median, std))
    
    
    
    if method == 'aperture':
    
        ###   LOCATING STARS AND PLOTTING THEM   ###
        # could maybe put most of this in its own starfind function? including aperture and position bits

        sources = star_find(data, sigma_level, FWHM, std, roundness, sharphi)
        #sources = star_find(data, 300., FWHM, std, roundness, sharphi)
        pos = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(pos, r_ap)
        annuli = CircularAnnulus(pos, r_in, r_out)
        print('No. of stars detected: {0}'.format(len(sources)))
        plotting(data, apertures, annuli, 'Greys', choice='both_aps', savename=None)
        
        #starlist = np.array(sources['xcentroid', 'ycentroid'])
        #np.savetxt('outputs/starlist_02.txt', starlist, delimiter=',')

        ###   DOING APERTURE PHOTOMETRY   ###

        phot_table = ap_phot(data, apertures, annuli, 'median')
        
        phot_table['apparent_mag'] = float('NaN')
        phot_table['mag_err'] = float('NaN')
        for i in range(len(phot_table)):
            locorr = data_corr[int(phot_table['ycenter'][i].value)][int(phot_table['xcenter'][i].value)]  # location-dependent correction at star's location
            if phot_table['star_bkgsub'][i] >= 0:
                phot_table['apparent_mag'][i] = zmag - 2.5 * math.log10(apcorr * locorr * phot_table['star_bkgsub'][i] / conv)
                phot_table['mag_err'][i] = apparent_magnitude_err(phot_table['star_bkgsub'][i], apcorr, zmag_err, phot_table['counts_err'][i], apcorr_err)
        
        #phot_table.write('outputs/table_'+channel+'_e'+str(epoch)+'.txt', format='csv', overwrite=True)
        #print_table_nicely(phot_table)
        print_table_nicely(phot_table['id', 'xcenter', 'ycenter', 'star_bkgsub', 'counts_err', 'apparent_mag', 'mag_err'])
        
        ###   COMPARING AGAINST MASTER LIST   ###
        
        cat_ep, ra_ep, dec_ep = make_catalog(phot_table['xcenter'], phot_table['ycenter'], hdr)
        phot_table['RA'] = ra_ep
        phot_table['dec'] = dec_ep
        
        start = time.perf_counter()
        
        for ma in range(len(cat_m)):
            print(ma)
            seps = []
            for ep in range(len(cat_ep)):
                seps.append(cat_ep[ep].separation(cat_m[ma]).value)
            id_good = np.asarray(seps).argmin()

        print('Ma: {0}\tEp: {1}\tId: {2}\tSep: {3}'.format(ma, ep, id_good, seps[id_good]))
        
        
        # this is waaaay too slow over two massive tables :(
        # plus, how would i handle multiple epoch stars matching the same master star??
        
        
        print('Time taken to calculate all separations: {}s'.format(time.perf_counter() - start))
        print(len(seps))
            
        
        
        
        ###   STAR MATCHING   ###
        
        if star_matching == 'yes':
        
            LC_time.append(Time(hdr['DATE_OBS'], format='isot', scale='utc').mjd)

            if epoch == comp_epoch:  # this may only work in this format if comp_epoch = 1 lol
                comp_table = copy.copy(phot_table)
                cat_comp, ra_comp, dec_comp = make_catalog(comp_table['xcenter'], comp_table['ycenter'], hdr)
                comp_table['RA'] = ra_comp
                comp_table['dec'] = dec_comp

                #comp_table['id', 'RA', 'dec', 'xcenter', 'ycenter', 'apparent_mag', 'mag_err'].write('outputs/comp_table_'+channel+'_e'+str(epoch)+'.txt', format='csv', overwrite=True, delimiter=' ')
            else:
                match_table = copy.copy(phot_table)
                cat_match, ra, dec = make_catalog(match_table['xcenter'], match_table['ycenter'], hdr, master)
                match_table['RA'] = ra
                match_table['dec'] = dec

                #match_table['id', 'RA', 'dec', 'xcenter', 'ycenter', 'apparent_mag', 'mag_err'].write('outputs/match_table_'+channel+'_e'+str(epoch)+'.txt', format='csv', overwrite=True, delimiter=' ')

                idx, d2d, d3d = cat_match.match_to_catalog_sky(cat_comp)

                # Selection criteria:
                max_sep = 0.25 * u.arcsec
                selection = (d2d > max_sep)
                match_index = idx
                match_index[selection] = -99.
                ind = ((match_index >= 0))  # keeps only positive indices (ie. removes -99s)
                print('Number of common stars:', sum(ind))

                match_table = match_table[ind]
                comp_table_new = comp_table[match_index][ind]  # creating a different object so comp_table is untouched for the next iteration
                delta_mag = comp_table_new['apparent_mag'] - match_table['apparent_mag']
                print_table_nicely(match_table['id', 'xcenter', 'ycenter', 'RA', 'dec', 'apparent_mag', 'mag_err'])
                print_table_nicely(comp_table_new['id', 'xcenter', 'ycenter', 'RA', 'dec', 'apparent_mag', 'mag_err'])

                plt.plot(delta_mag, match_table['apparent_mag'], 'b+')
                #plt.plot(delta_mag, mag1, 'b+')
                plt.xlabel('Difference in magnitude between epochs')
                plt.ylabel('Magnitude')
                plt.grid()
                plt.gca().invert_yaxis()
                plt.gcf().set_size_inches(10, 8)
                plt.show()

                #LC_mag.append(match_table['apparent_mag'][0])  # doesn't currently include first epoch mag

                c1 = SkyCoord(match_table['RA'], match_table['dec'], frame='icrs', unit='deg')
                c2 = SkyCoord(comp_table_new['RA'], comp_table_new['dec'], frame='icrs', unit='deg')
                sep = c1.separation(c2)
                print('Largest separation between matching RA and dec coords: {} arcsec'.format(max(sep.arcsecond)))
        
        
        
    
    elif method == 'PSF':
        
        ###   DETECTING MORE STARS   ###
        
        daofind_PSF = DAOStarFinder(threshold=sigma_level*std, fwhm=FWHM, roundlo=-roundness, roundhi=roundness, sharphi=0.85)
        sources_PSF = daofind_PSF(data)
        
        #sources_PSF = star_find(data, sigma_level, FWHM, std, roundness, sharphi)
        pos_PSF = np.transpose((sources_PSF['xcentroid'], sources_PSF['ycentroid']))
        ap_PSF = CircularAperture(pos_PSF, r_ap)
        plotting(data, ap_PSF, an=None, cmap='Greys', choice='ap_only', savename='psf_detect_test_00.png')
        print('Number of PSF stars: {}'.format(len(sources_PSF)))
        
        ###   GROUPING STARS BASED ON PROXIMITY TO NEIGHBOURS   ###
        
        sources_PSF['xcentroid'].name = 'x_0'
        sources_PSF['ycentroid'].name = 'y_0'
        daogroup = DAOGroup(crit_separation=2.5*FWHM)
        star_groups = daogroup(sources_PSF)
        star_groups = star_groups.group_by('group_id')
        print_table_nicely(star_groups)
        
        ncolors = max(star_groups['group_id'])
        cmap = make_random_cmap(ncolors=ncolors, seed=1612)
        plt.imshow(data, origin='lower', norm=LogNorm(), interpolation='nearest', cmap='Greys')
        for i, group in enumerate(star_groups.groups):
            pos_group = np.transpose([group['x_0'], group['y_0']])
            ap_group = CircularAperture(pos_group, r_ap)
            ap_group.plot(color=cmap.colors[i], lw=1.)
        plt.gcf().set_size_inches(15, 9)
        plt.show()
        
        ###   DOING PSF PHOTOMETRY   ###
        
        data_nonans = np.nan_to_num(data, nan=0.00001, copy=True) # changing NaN values to a float so PSF fit doesn't crash
        
        #bkg_estimation = MMMBackground()
        #fitter = LevMarLSQFitter()
        fit_rad = 5
        
        epsf.x_0.fixed = True
        epsf.y_0.fixed = True
        init_pos = Table(names=['x_0', 'y_0'], data=[sources_PSF['x_0'], sources_PSF['y_0']])
        
        start = time.perf_counter()
        PSF_photometry = ISPSF(finder = daofind_PSF,
                              group_maker = daogroup,
                              bkg_estimator = MMMBackground(),
                              psf_model = epsf,
                              fitter = LevMarLSQFitter(),
                              fitshape = 2*fit_rad+1,
                              niters = 10,
                              aperture_radius = 6.)
        PSF_table = PSF_photometry(image=data_nonans, init_guesses=init_pos)
        residual_image = PSF_photometry.get_residual_image()
        print('Time taken to fit PSF model: {}s'.format(time.perf_counter() - start))
        
        plotting(residual_image, ap=None, an=None, cmap='viridis', choice=None, savename=None)
        
        PSF_table['apparent_mag'] = float('NaN')
        PSF_table['mag_err'] = float('NaN')
        for i in range(len(PSF_table)):
            locorr = data_corr[int(PSF_table['y_fit'][i])][int(PSF_table['x_fit'][i])]  # location-dependent correction at star's location
            if PSF_table['flux_fit'][i] >= 0:
                PSF_table['apparent_mag'][i] = zmag - 2.5 * math.log10(apcorr * locorr * PSF_table['flux_fit'][i] / conv)
                PSF_table['mag_err'][i] = apparent_magnitude_err(PSF_table['flux_fit'][i], apcorr, zmag_err, PSF_table['flux_unc'][i], apcorr_err)
        
        
        cat, ra, dec = make_catalog(PSF_table['x_fit'], PSF_table['y_fit'], hdr)
        PSF_table['RA'] = ra
        PSF_table['dec'] = dec
        
        #print_table_nicely(PSF_table['id', 'group_id', 'iter_detected', 'x_fit', 'y_fit', 'apparent_mag', 'mag_err'])
        print_table_nicely(PSF_table['id', 'x_fit', 'y_fit', 'RA', 'dec', 'apparent_mag', 'mag_err'])
        #PSF_table['id', 'RA', 'dec', 'x_fit', 'y_fit', 'apparent_mag', 'mag_err'].write('outputs/matching/'+method+'_'+channel+'_matching_e'+str(epoch)+'.txt', format='csv', overwrite=True, delimiter=' ')
        
        LC_time.append(Time(hdr['DATE_OBS'], format='isot', scale='utc').mjd)
        
        
    print('\n\n\n')
    
    if run_all == False:
        break
    
    #if epoch == 2:
    #    break
    
    epoch += 1