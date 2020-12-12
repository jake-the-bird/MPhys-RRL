### Script generously provided by Vicky Scowcroft ###
### Converted to Python 3.x and edited by Jake Bird, December 1 2020 ###




import numpy as np
import matplotlib.pyplot as plt
import sys
import gloess_fits as gf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
#matplotlib.rc('text',usetex=True)
from matplotlib import rcParams
#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Garamond']

def runGloess(mag_ch1, err_ch1, mag_ch2, err_ch2, LC_time, period, starname, xir1, xir2, wantcolour):

    ## Converting the gloess fourtran/pgplot code to python/matplotlib
    ## June 15 2012

    ## Version 1.0
    ## last edit - June 19 2012

    ## Next thing to add:
    ##Print fits to an output text file


    ## Open the input data file and read the info

    #input = sys.argv[1]
    #input = 'c:/Users/Jake/MPhys-code/MPhys-RRL/test_photometry/gloess-master/RRLyr.gloess_in'
    #input = 'c:/Users/Jake/MPhys-code/MPhys-RRL/test_photometry/gloess-master/UVOct.gloess_in'
    #counter = 0

    ## Want to know whether the IRAC data is phased or not. 
    ## If it is phased, must reduce the uncertainty by another factor of sqrt(N)

    #nir1 = sum(ir1<50)
    #nir2= sum(ir2<50)

    phased = 0   # our data is NOT phased

    ir1 = np.array(mag_ch1)
    #nir1 = len(ir1)
    nir1 = sum(ir1<50)
    eir1 = np.array(err_ch1)
    #xir1 = 0.10

    ir2 = np.array(mag_ch2)
    #nir2 = len(ir2)
    nir2 = sum(ir2<50)
    eir2 = np.array(err_ch2)
    #xir2 = 0.10

    mjd = np.array(LC_time)

    # Phases don't need to be done individually by band - only depends on P
    phase = (mjd / period) - np.floor(mjd / period)
    #phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))

    # Usage:  fit_one_band(data,err,phases,n,smooth):
    maxvals = []
    minvals = []

    if nir1 > 0:
        maxvals.append(np.amax(ir1[ir1<50])-1.4)
        minvals.append(np.amin(ir1[ir1<50])-1.4)
    if nir2 > 0:
        maxvals.append(np.amax(ir2[ir2<50])-1.8)
        minvals.append(np.amin(ir2[ir2<50])-1.8)


    maxvals = np.array(maxvals)
    minvals = np.array(minvals)

    max = np.max(maxvals)
    min = np.min(minvals)
    print(starname, '---- Period =', period, 'days')
    print('------------------------------------------------------')

    # Set up names for output files

    #fitname = starname + '.glo_fits'
    #avname = str(starname) + '.glo_avs'

    #avsout = open(avname,'w')
    #fitout = open(fitname,'w')

    maxlim = max + 0.5
    minlim = min - 0.5

    plt.clf()

    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)
    #plt.figure(figsize=(16.0,10.0))

    gs = gridspec.GridSpec(3, 4)
    ax1 = plt.subplot(gs[:, 0:2])
    ax2 = plt.subplot(gs[0, 2:4])
    ax3 = plt.subplot(gs[1, 2:4])
    ax4 = plt.subplot(gs[2, 2:4])
    ax1.axis([1,3.5,(maxlim),(minlim)])
    titlestring = str(starname) + ', P = ' + str(period) + ' days'
    #print titlestring
    plt.suptitle(titlestring, fontsize=20)

    ax1.set_ylabel('Magnitude')
    ax1.set_xlabel('Phase $\phi$')

    ## Fitting and plotting for each band
    print(nir1, nir2)

    if nir1 > 0:
        ir11, ir1x, yir1, yeir1, xphaseir1 = gf.fit_one_band(ir1,eir1,phase,nir1,xir1)
        ax1.plot(ir1x,ir11-1.4,'k-')
        ax1.plot(xphaseir1,yir1-1.4,color='MediumVioletRed',marker='o',ls='None', label='$[3.6]-1.4$')
    ## for RRLyrae WISE plots:
    #	ax1.plot(ir1x,ir11+1.,'k-')
    # 	ax1.plot(xphaseir1,yir1+1.,color='Turquoise',marker='o',ls='None', label='W1+1.0')
        aveir1, adevir1, sdevir1, varir1, skewir1, kurtosisir1, ampir1 = gf.moment(ir11[200:300],100)
        if phased == 1:
            factor = np.sqrt(nir1)
        if phased == 0:
            factor = 1 
        if nir1 > 1:
            #print('<[3.6]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f} N ir1 = {3}'.format(aveir1, sdevir1/factor, ampir1,nir1), file = avsout)
            print('<[3.6]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir1, sdevir1/factor, ampir1))
        if nir1 == 1:
            #print('[3.6] = {0:.3f} --- single point'.format(aveir1), file=avsout)
            print('[3.6] = {0:.3f} --- single point'.format(aveir1))
            
        ax2.axis([1,3.5,(np.average(ir11[200:300]) + 0.4),(np.average(ir11[200:300]) - 0.4)])
        ax2.yaxis.tick_right()
        ax2.plot(ir1x,ir11,'k-')
        ax2.plot(xphaseir1,yir1,color='MediumVioletRed',marker='o',ls='None', label='$[3.6]$')
        ax2.annotate('$[3.6]$', xy=(0.04, 0.8375), xycoords='axes fraction', fontsize=16)
    else:
        aveir1 = float("NaN")
        ampir1 = float("NaN")
        sdevir1 = float("NaN")

    if nir2 > 0:
        ir21, ir2x, yir2, yeir2, xphaseir2 = gf.fit_one_band(ir2,eir2,phase,nir2,xir2)
        ax1.plot(ir2x,ir21-1.8,'k-')
        ax1.plot(xphaseir2,yir2-1.8,color='DeepPink',marker='o',ls='None', label='$[4.5]-1.8$')
    ## For RRLyrae WISE plots:
    #	ax1.plot(ir2x,ir21,'k-')
    # 	ax1.plot(xphaseir2,yir2,color='Gold',marker='o',ls='None', label='W2')
        aveir2, adevir2, sdevir2, varir2, skewir2, kurtosisir2, ampir2= gf.moment(ir21[200:300],100)
        if phased == 1:
            factor = np.sqrt(nir2)
        if phased == 0:
            factor = 1
        if nir2 > 1:
            #print('<[4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f} N IR2 = {3}'.format(aveir2, sdevir2/factor, ampir2,nir2), file = avsout)
            print('<[4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir2, sdevir2/factor, ampir2))
        if nir2 == 1:
            #print('[4.5] = {0:.3f} --- single point'.format(aveir2), file = avsout)
            print('[4.5] = {0:.3f} --- single point'.format(aveir2))
            
        ax3.axis([1,3.5,(np.average(ir21[200:300]) + 0.4),(np.average(ir21[200:300]) - 0.4)])
        ax3.yaxis.tick_right()
        ax3.plot(ir2x,ir21,'k-')
        ax3.plot(xphaseir2,yir2,color='DeepPink',marker='o',ls='None', label='$[3.6]$')
        ax3.annotate('$[4.5]$', xy=(0.04, 0.8375), xycoords='axes fraction',fontsize=16)
    else:
        aveir2 = float("NaN")
        ampir2 = float("NaN")
        sdevir2 = float("NaN")


    handles, labels = ax1.get_legend_handles_labels() 
    #ax1.legend(handles[::-1],labels[::-1],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., numpoints=1)
    ax1.legend(handles[::-1],labels[::-1],loc=4, numpoints=1,prop={'size':10})

    #plt.setp(ax1.get_xticklabels(),visible=False)

    if (wantcolour == 'yes' or wantcolour == True):
        if nir1 == nir2:
            ### Define the colour curve
            colour_curve = ir11 - ir21
            ## Define the colour points
            ch1_points = yir1[yir1<99]
            ch2_points = yir2[yir2<99]
            colour_points = ch1_points - ch2_points
            colour_phases = xphaseir1[yir1<99]

            colour_points = np.concatenate((colour_points,colour_points,colour_points,colour_points,colour_points))
            colour_phases = np.concatenate((colour_phases,(colour_phases+1.),(colour_phases+2.),(colour_phases+3.),(colour_phases+4.)))


            avecol, adevcol, sdevcol, varcol, skewcol, kurtosiscol, ampcol = gf.moment(colour_curve[200:300],100)

            #print('<[3.6] - [4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(avecol, sdevcol/factor, ampcol), file = avsout)
            print('<[3.6] - [4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(avecol, sdevcol/factor, ampcol))

            print(np.average(ir11[200:300]) + 0.3)
            print(np.average(ir11[200:300]) - 0.3)


            #divider = make_axes_locatable(ax1)
            #axcol = divider.append_axes("bottom",1.2,pad=0.1,sharex=ax1)
            myaxis2 = [1,3.5,-0.2,0.2]
            ax4.axis(myaxis2)
            ax4.yaxis.tick_right()
            ax4.yaxis.set_major_locator(plt.FixedLocator([-0.1,0,0.1]))
            ax4.plot(ir1x,colour_curve,'k-')
            ax4.plot(colour_phases,colour_points,color='Black',marker='o',ls='None', label='$[3.6]-[4.5]$')

            ax4.set_xlabel('Phase $\phi$')
            #ax4.annotate('$[3.6] - [4.5]$', xy=(1.1, 0.135), xycoords='data')
            ax4.annotate('$[3.6] - [4.5]$', xy=(0.04, 0.8375), xycoords='axes fraction',fontsize=16)

            ax4.hlines(0,1,3.5,'k','dashdot')
        else:
            print("Channels don't have the same number of measurements, so colour curve couldn't be plotted")

    plt.setp(ax2.get_xticklabels(),visible=False)
    plt.setp(ax3.get_xticklabels(),visible=False)

    plotname = str(starname)+'.eps'
    #plt.savefig(plotname, transparent='True')
    plt.savefig(plotname)

    #avsout.close()
    plt.show()
    #fitout.close()

    return aveir1, ampir1, sdevir1, aveir2, ampir2, sdevir2