# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

# __::((odmkScipyFIR.py))::__

# Python ODMK timing sequencer module

# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header end-------------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

import os
import sys
import numpy as np
#from numpy import cos, sin, pi, absolute, arange
from scipy.signal import kaiserord, lfilter, firwin, freqz
from pylab import figure, clf, plot, xlabel, ylabel, xlim, ylim, title, grid, axes, show
import matplotlib.pyplot as plt



rootDir = 'C:/odmkDev/odmkCode/odmkPython/'
audioSrcDir = 'C:/odmkDev/odmkCode/odmkPython/audio/multiWavSrcDir/'
audioOutDir = 'C:/odmkDev/odmkCode/odmkPython/audio/wavout/'

defaultTxtOutDir = 'C:/odmkDev/odmkCode/odmkPython/DSP/firTestOutDir/'


sys.path.insert(0, rootDir+'util')
from odmkClear import *
#from odmkPlotUtil import *
import odmkPlotUtil as odmkplt

sys.path.insert(1, rootDir+'audio')
import odmkWavIO as waveio

sys.path.insert(2, rootDir+'DSP')
import odmkClocks as clks
import odmkSigGen1 as sigGen

# temp python debugger - use >>>pdb.set_trace() to set break
# import pdb


# // *---------------------------------------------------------------------* //
plt.close('all')
clear_all()

# // *---------------------------------------------------------------------* //

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : function definitions
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# #########################################################################
# begin : file output
# #########################################################################

# // *-----------------------------------------------------------------* //
# // *---TXT write simple periodic waveforms (sigLength # samples)
# // *-----------------------------------------------------------------* //

def sig2txt(sigIn, nChan, outNm, outDir='None'):
    ''' writes data to TXT file
        signal output name = outNm (expects string) '''

    if outDir != 'None':
        try:
            if isinstance(outDir, str):
                txtOutDir = outDir
                os.makedirs(txtOutDir, exist_ok=True)
        except NameError:
            print('Error: outNm must be a string')
    else:
        txtOutDir = defaultTxtOutDir

    try:
        if isinstance(outNm, str):
            sigOutFull = txtOutDir+outNm
    except NameError:
        print('Error: outNm must be a string')

    # writes data to .TXT file:
    outputFile = open(sigOutFull, 'w', newline='')

    if nChan == 0:
        print('ERROR: Number of Channels must be >= 1')
    elif nChan == 1:
        for i in range(len(sigIn)):
            outputFile.write(str(sigIn[i]) + '\n')
    else:
        for i in range(len(sigIn[0])):
            lineTmp = ""
            for j in range(len(sigIn) - 1):
                strTmp = str(sigIn[j, i]) + str('    ')
                lineTmp = lineTmp + strTmp
            lineTmp = lineTmp + str(sigIn[len(sigIn) - 1, i]) + '\n'
            outputFile.write(lineTmp)            

    outputFile.close()

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# end : function definitions
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# // *---------------------------------------------------------------------* //


#------------------------------------------------
# Create a signal for demonstration.
#------------------------------------------------

sample_rate = 48000.0
nsamples = 1024

#sample_rate = 100.0
#nsamples = 400
#t = arange(nsamples) / sample_rate
#x = cos(2*pi*0.5*t) + 0.2*sin(2*pi*2.5*t+0.1) + \
#        0.2*sin(2*pi*15.3*t) + 0.1*sin(2*pi*16.7*t + 0.1) + \
#            0.1*sin(2*pi*23.45*t+.8)


#t = arange(nsamples) / sample_rate
#x = cos(2*pi*0.5*t) + 0.2*sin(2*pi*2.5*t+0.1) + \
#        0.2*sin(2*pi*15.3*t) + 0.1*sin(2*pi*16.7*t + 0.1) + \
#            0.1*sin(2*pi*23.45*t+.8)
            
F_1KHz = 560.
A_1KHz = 1.0
 
F_15KHz = 9000.
A_15KHz = 0.5
 
t = np.arange(nsamples) / sample_rate
x = A_1KHz * np.sin(2*np.pi*F_1KHz*t) + A_15KHz*np.sin(2*np.pi*F_15KHz*t)

#x = 3*cos(2*pi*t)+cos(2*pi*3*t)+2*cos(2*pi*5*t)

# write to output file
outNm = 'input.dat'
sig2txt(x, 1, outNm)


#------------------------------------------------
# Create a FIR filter and apply it to x.
#------------------------------------------------

# The Nyquist rate of the signal.
nyq_rate = sample_rate / 2.0

# The cutoff frequency of the filter.
cutoff_hz = 6000.

numCoeff = 47

# *---can't control # of taps with this method---*
# The desired width of the transition from pass to stop,
# relative to the Nyquist rate.  We'll design the filter
# with a 5 Hz transition width.
#width = 5.0/nyq_rate
#width = 1000.0/nyq_rate
# The desired attenuation in the stop band, in dB.
#ripple_db = 60.0
# Compute the order and Kaiser parameter for the FIR filter.
#N, beta = kaiserord(ripple_db, width)
# Use firwin with a Kaiser window to create a lowpass FIR filter.
#taps = firwin(N, cutoff_hz/nyq_rate, window=('kaiser', beta))


## Use firwin to create a lowpass FIR filter
fir_coeff = firwin(numCoeff, cutoff_hz/nyq_rate)

finalTaps = ""
for j in range(int(len(fir_coeff)/4)):
    print("%.16f" % fir_coeff[j*4]+",    %.16f" % fir_coeff[j*4+1]+",    %.16f" % fir_coeff[j*4+2]+",    %.16f" % fir_coeff[j*4+3]+",    ")
for k in range(len(fir_coeff)%4):
    finalTaps = finalTaps+"%.16f" % fir_coeff[len(fir_coeff)-(len(fir_coeff)%4-k)]+",    "
print(finalTaps)

# write to output file
outNm = 'fir_coeff.dat'
sig2txt(fir_coeff, 1, outNm)


# Use lfilter to filter x with the FIR filter.
filtered_x = lfilter(fir_coeff, 1.0, x)

# write to output file
outNm = 'ref_res.dat'
sig2txt(filtered_x, 1, outNm)

#------------------------------------------------
# Plot the FIR filter coefficients.
#------------------------------------------------

figure(1)
plot(fir_coeff, 'bo-', linewidth=2)
title('Filter Coefficients (%d taps)' % numCoeff)
grid(True)

#------------------------------------------------
# Plot the magnitude response of the filter.
#------------------------------------------------

figure(2)
clf()
w, h = freqz(fir_coeff, worN=8000)
plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
xlabel('Frequency (Hz)')
ylabel('Gain')
title('Frequency Response')
ylim(-0.05, 1.05)
grid(True)

# Upper inset plot.
ax1 = axes([0.43, 0.55, .37, .25])
plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
xlim(0,7000.0)
ylim(0.9965, 1.007)
grid(True)

# Lower inset plot
ax2 = axes([0.43, 0.23, .37, .25])
plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
xlim(7000.0, 15000.0)
ylim(0.0, 0.0125)
grid(True)

#------------------------------------------------
# Plot the original and filtered signals.
#------------------------------------------------

# The phase delay of the filtered signal.
delay = 0.5 * (numCoeff-1) / sample_rate

figure(3)
# Plot the original signal.
plot(t, x)
# Plot the filtered signal, shifted to compensate for the phase delay.
plot(t-delay, filtered_x, 'r-')
# Plot just the "good" part of the filtered signal.  The first N-1
# samples are "corrupted" by the initial conditions.
plot(t[numCoeff-1:]-delay, filtered_x[numCoeff-1:], 'g', linewidth=4)

xlabel('t')
grid(True)

show()




#from numpy import sin, arange, pi
#from scipy.signal import lfilter, firwin
#from pylab import figure, plot, grid, show
# 
##------------------------------------------------
## Create a signal for demonstration.
##------------------------------------------------
## 320 samples of (1000Hz + 15000 Hz) at 48 kHz
#sample_rate = 48000.
#nsamples = 320
# 
#F_1KHz = 1000.
#A_1KHz = 1.0
# 
#F_15KHz = 15000.
#A_15KHz = 0.5
# 
#t = arange(nsamples) / sample_rate
#signal = A_1KHz * sin(2*pi*F_1KHz*t) + A_15KHz*sin(2*pi*F_15KHz*t)
# 
##------------------------------------------------
## Create a FIR filter and apply it to signal.
##------------------------------------------------
## The Nyquist rate of the signal.
#nyq_rate = sample_rate / 2.
# 
## The cutoff frequency of the filter: 6KHz
#cutoff_hz = 6000.0
# 
## Length of the filter (number of coefficients, i.e. the filter order + 1)
#numtaps = 29
# 
## Use firwin to create a lowpass FIR filter
#fir_coeff = firwin(numtaps, cutoff_hz/nyq_rate)
# 
## Use lfilter to filter the signal with the FIR filter
#filtered_signal = lfilter(fir_coeff, 1.0, signal)
# 
##------------------------------------------------
## Plot the original and filtered signals.
##------------------------------------------------
# 
## The first N-1 samples are "corrupted" by the initial conditions
#warmup = numtaps - 1
# 
## The phase delay of the filtered signal
#delay = (warmup / 2) / sample_rate
# 
#figure(1)
## Plot the original signal
#plot(t, signal)
# 
## Plot the filtered signal, shifted to compensate for the phase delay
#plot(t-delay, filtered_signal, 'r-')
# 
## Plot just the "good" part of the filtered signal.  The first N-1
## samples are "corrupted" by the initial conditions.
#plot(t[warmup:]-delay, filtered_signal[warmup:], 'g', linewidth=4)
# 
#grid(True)
# 
#show()
# 
##------------------------------------------------
## Print values
##------------------------------------------------
#def print_values(label, values):
#    var = "float32_t %s[%d]" % (label, len(values))
#    print("%-30s = {%s}" % (var, ', '.join(["%+.10f" % x for x in values])))
# 
#print_values('signal', signal)
#print_values('fir_coeff', fir_coeff)
#print_values('filtered_signal', filtered_signal)