# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkPVx1_tb.py))::__
#
# Python testbench for odmkClocks, odmkSigGen1 objects
# required lib:
# odmkClocks ; odmkSigGen1
#
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header end-------------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

import os
import csv
import wave
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from odmkClear import *

import odmkClocks as clks
import odmkSigGen1 as sigGen

# temp python debugger - use >>>pdb.set_trace() to set break
import pdb

# // *---------------------------------------------------------------------* //
clear_all()

# // *---------------------------------------------------------------------* //

print('// //////////////////////////////////////////////////////////////// //')
print('// *--------------------------------------------------------------* //')
print('// *---::ODMK Signal Generator 1::---*')
print('// *--------------------------------------------------------------* //')
print('// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ //')


# // *---------------------------------------------------------------------* //

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : function definitions
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# // *---------------------------------------------------------------------* //
# // *--Math Functions--*
# // *---------------------------------------------------------------------* //

def cyclicZn(n):
    ''' calculates the Zn roots of unity '''
    cZn = np.zeros((n, 1))*(0+0j)    # column vector of zero complex values
    for k in range(n):
        # z(k) = e^(((k)*2*pi*1j)/n)        # Define cyclic group Zn points
        cZn[k] = np.cos(((k)*2*np.pi)/n) + np.sin(((k)*2*np.pi)/n)*1j   # Euler's identity

    return cZn


# // *---------------------------------------------------------------------* //
# // *--Plot Functions--*
# // *---------------------------------------------------------------------* //

def odmkPlot1D(fnum, sig, xLin, pltTitle, pltXlabel, pltYlabel, lncolor='red', lnstyle='-', lnwidth=1.00, pltGrid=True, pltBgColor='black'):
    ''' ODMK 1D Matplotlib plot
        required inputs:
        fnum => unique plot number
        sig => signal to plot
        xLin => linear space to define x-axis (0 to max x-axis length-1)
        pltTitle => text string for plot title
        pltXlabel => text string for x-axis
        pltYlabel => text string for y-axis
        optional inputs:
        lncolor => line color (default = red ; html color names, html color codes??)
        lnstyle => line style (default = plain line ; * ; o ; etc..)
        lnwidth => line width
        pltGrid => use grid : default = True ; <True;False>
        pltBgColor => backgroud color (default = black) '''

    # Input signal
    plt.figure(num=fnum, facecolor='silver', edgecolor='k')
    # check if xLin is < than or = to sig
    if len(xLin) > len(sig):
        print('ERROR: length of xLin x-axis longer than signal length')
        return 1
    elif len(xLin) == len(sig):
        odmkMatPlt = plt.plot(xLin, sig)
    else:
        odmkMatPlt = plt.plot(xLin, sig[0:len(xLen)])

    plt.setp(odmkMatPlt, color=lncolor, ls=lnstyle, linewidth=lnwidth)
    plt.xlabel(pltXlabel)
    plt.ylabel(pltYlabel)
    plt.title(pltTitle)
    plt.grid(color='c', linestyle=':', linewidth=.5)
    plt.grid(pltGrid)
    # plt.xticks(np.linspace(0, Fs/2, 10))
    ax = plt.gca()
    ax.set_axis_bgcolor(pltBgColor)

    return 0
    
def odmkMultiPlot1D(fnum, sigArray, xLin, pltTitle, pltXlabel, pltYlabel, colorMp='gnuplot', lnstyle='-', lnwidth=1.00, pltGrid=True, pltBgColor='black'):
    ''' ODMK 1D Matplotlib multi-plot
        required inputs:
        fnum => unique plot number
        sig => signal to plot : 2D Numpy array
        xLin => linear space to define x-axis (0 to max x-axis length-1)
        pltTitle => text string for plot title
        pltXlabel => text string for x-axis
        pltYlabel => text string for y-axis
        optional inputs:
        lncolor => line color (default = red ; html color names, html color codes??)
        lnstyle => line style (default = plain line ; * ; o ; etc..)
        lnwidth => line width
        pltGrid => use grid : default = True ; <True;False>
        pltBgColor => backgroud color (default = black) '''

    # define the color map
    try:
        cmap = plt.cm.get_cmap(colorMp)
    except ValueError as e:
        print('ValueError: ', e)
    colors = cmap(np.linspace(0.0, 1.0, len(sigArray[0, :])))

    # Input signal
    plt.figure(num=fnum, facecolor='silver', edgecolor='k')
    # check if xLin is < than or = to sig
    if len(xLin) > len(sigArray[:, 0]):
        print('ERROR: length of xLin x-axis longer than signal length')
        return 1
    else:
        if len(xLin) == len(sigArray[:, 0]):
            # odmkMatPlt = []
            for i in range(len(sinArray[0, :])):
                plt.plot(xLin, sigArray[:, i], color=colors[i], ls=lnstyle, linewidth=lnwidth)
        else:
            # odmkMatPlt = []
            for i in range(len(sinArray[0, :])):
                plt.plot(xLin, sigArray[0:len(xLin), i], color=colors[i], ls=lnstyle, linewidth=lnwidth)

        plt.xlabel(pltXlabel)
        plt.ylabel(pltYlabel)
        plt.title(pltTitle)
        plt.grid(color='c', linestyle=':', linewidth=.5)
        plt.grid(pltGrid)
        # plt.xticks(np.linspace(0, Fs/2, 10))
        ax = plt.gca()
        ax.set_axis_bgcolor(pltBgColor)

    return 0


# // *---------------------------------------------------------------------* //
# // *--Phase Vocoder Functions--*
# // *---------------------------------------------------------------------* //

def odmkPVA(PVAin, Awin, NFFT, Ra, Fs):
    ''' odmk phase vocoder Analysis
        PVAin -> wave file / signal source in
        Awin -> phase vocoder analysis window (default = cos)
        PVAout -> phase vocoder output array <rows = bins, cols = frames>
        NFFT -> fft length
        Ra -> Analysis sample hop
        fs -> sampling frequency '''

    sigLength = len(PVAin)

    # int posin, posout, i, k, mod;
    # float *sigframe, *specframe, *lastph;
    # float fac, scal, phi, mag, delta

    sigframe = np.array(NFFT)
    specframe = np.array(NFFT)

    lastph = np.zeros(NFFT/2)

    fac = Fs / (Ra * 2*np.pi)
    scal = 2*np.pi * Ra/NFFT

    inIdx = 0
    outIdx = 0
    # for(posin=posout=0; posin < input_size; posin+=hopsize):
    for j in range(0, sigLength, Ra):
        mod = inIdx % NFFT
        # window & rotate a signal frame
        for i in range(NFFT):
            if (inIdx+i < sigLength):
                sigframe[(i + mod) % NFFT] = PVAin[inIdx + i] * Awin[i]
            else:
                sigframe[(i + mod) % NFFT] = 0

        # transform it
        specframe = sp.fft(sigframe)
        # rectangular to polar
        # specMag = np.abs(specframe)
        # specPhase = np.arctan2(specframe.imag, specframe.real)

        # convert to PV output
        k = 1
        h = 2
        # for(h=2,k=1; h < fftsize; h+=2, k++){
        for h in range(0, NFFT, 2):

            mag = sp.sqrt(specframe[h]*specframe[h] + specframe[h+1]*specframe[h+1])
            phi = np.atan2(specframe[h+1], specframe[h])
            # phase diffs
            delta = phi - lastph[k]
            lastph[k] = phi

            # unwrap the difference, so it lies between -pi and pi
            while (delta > np.pi):
                delta -= 2 * np.pi
            while (delta < -np.pi):
                delta += 2 * np.pi

            # construct the amplitude-frequency pairs
            specframe[h] = mag
            specframe[h+1] = (delta + k * scal) * fac

            k += 1

        # output it
        # for(i=0; i < fftsize; i++, posout++)
        for g in range(NFFT):
            output[outIdx] = specframe[g]
            outIdx += 1

    return PVAout


#int pvs(float* input, float* window, float* output,
#          int input_size, int fftsize, int hopsize, float sr){
#
#int posin, posout, k, i, output_size, mod;
#float *sigframe, *specframe, *lastph;
#float fac, scal, phi, mag, delta;
#
#sigframe = new float[fftsize];
#specframe = new float[fftsize];
#lastph = new float[fftsize/2];
#memset(lastph, 0, sizeof(float)*fftsize/2);
#
#output_size = input_size*hopsize/fftsize;
#
#fac = (float) (hopsize*twopi/sr);
#scal = sr/fftsize;
#
#for(posout=posin=0; posout < output_size; posout+=hopsize){ 
#
#   // load in a spectral frame from input 
#   for(i=0; i < fftsize; i++, posin++)
#        specframe[i] = input[posin];
#	
# // convert from PV input to DFT coordinates
# for(i=2,k=1; i < fftsize; i+=2, k++){
#   delta = (specframe[i+1] - k*scal)*fac;
#   phi = lastph[k]+delta;
#   lastph[k] = phi;
#   mag = specframe[i];
#  
#  specframe[i] = (float) (mag*cos(phi));
#  specframe[i+1] = (float) (mag*sin(phi)); 
#  
#}
#   // inverse-transform it
#   ifft(specframe, sigframe, fftsize);
#
#   // unrotate and window it and overlap-add it
#   mod = posout%fftsize;
#   for(i=0; i < fftsize; i++)
#       if(posout+i < output_size)
#          output[posout+i] += sigframe[(i+mod)%fftsize]*window[i];
#}
#delete[] sigframe;
#delete[] specframe;
#delete[] lastph;
#
#return output_size;
#}


#int pva(float *input, float *window, float *output, 
#        int input_size, int fftsize, int hopsize, float sr){
#
#int posin, posout, i, k, mod;
#float *sigframe, *specframe, *lastph;
#float fac, scal, phi, mag, delta, pi = (float)twopi/2;
#
#sigframe = new float[fftsize];
#specframe = new float[fftsize];
#lastph = new float[fftsize/2];
#memset(lastph, 0, sizeof(float)*fftsize/2);
#
#fac = (float) (sr/(hopsize*twopi));
#scal = (float) (twopi*hopsize/fftsize);
#
#for(posin=posout=0; posin < input_size; posin+=hopsize){
#      mod = posin%fftsize;
#	// window & rotate a signal frame
#      for(i=0; i < fftsize; i++) 
#          if(posin+i < input_size)
#            sigframe[(i+mod)%fftsize]
#                     = input[posin+i]*window[i];
#           else sigframe[(i+mod)%fftsize] = 0;
#
#      // transform it
#      fft(sigframe, specframe, fftsize);
#
#      // convert to PV output
#      for(i=2,k=1; i < fftsize; i+=2, k++){
#
#      // rectangular to polar
#      mag = (float) sqrt(specframe[i]*specframe[i] + 
#                        specframe[i+1]*specframe[i+1]);  
#      phi = (float) atan2(specframe[i+1], specframe[i]);
#      // phase diffs
#      delta = phi - lastph[k];
#      lastph[k] = phi;
#         
#      // unwrap the difference, so it lies between -pi and pi
#      while(delta > pi) delta -= (float) twopi;
#      while(delta < -pi) delta += (float) twopi;
#
#      // construct the amplitude-frequency pairs
#      specframe[i] = mag;
#	  specframe[i+1] = (delta + k*scal)*fac;
#
#      }
#      // output it
#      for(i=0; i < fftsize; i++, posout++)
#			  output[posout] = specframe[i];
#		  
#}
#delete[] sigframe;
#delete[] specframe;
#delete[] lastph;
#
#return posout;
#}
#
#int pvs(float* input, float* window, float* output,
#          int input_size, int fftsize, int hopsize, float sr){
#
#int posin, posout, k, i, output_size, mod;
#float *sigframe, *specframe, *lastph;
#float fac, scal, phi, mag, delta;
#
#sigframe = new float[fftsize];
#specframe = new float[fftsize];
#lastph = new float[fftsize/2];
#memset(lastph, 0, sizeof(float)*fftsize/2);
#
#output_size = input_size*hopsize/fftsize;
#
#fac = (float) (hopsize*twopi/sr);
#scal = sr/fftsize;
#
#for(posout=posin=0; posout < output_size; posout+=hopsize){ 
#
#   // load in a spectral frame from input 
#   for(i=0; i < fftsize; i++, posin++)
#        specframe[i] = input[posin];
#	
# // convert from PV input to DFT coordinates
# for(i=2,k=1; i < fftsize; i+=2, k++){
#   delta = (specframe[i+1] - k*scal)*fac;
#   phi = lastph[k]+delta;
#   lastph[k] = phi;
#   mag = specframe[i];
#  
#  specframe[i] = (float) (mag*cos(phi));
#  specframe[i+1] = (float) (mag*sin(phi)); 
#  
#}
#   // inverse-transform it
#   ifft(specframe, sigframe, fftsize);
#
#   // unrotate and window it and overlap-add it
#   mod = posout%fftsize;
#   for(i=0; i < fftsize; i++)
#       if(posout+i < output_size)
#          output[posout+i] += sigframe[(i+mod)%fftsize]*window[i];
#}
#delete[] sigframe;
#delete[] specframe;
#delete[] lastph;
#
#return output_size;
#}



# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# end : function definitions
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : Loading & Formatting img and sound files
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


rootDir = 'C:/usr/eschei/odmkPython/odmk/werk/'


print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Set Master Parameters for output::---*')
print('// *--------------------------------------------------------------* //')


# // *---------------------------------------------------------------------* //
# // *--Primary parameters--*
# // *---------------------------------------------------------------------* //

# length of x in seconds:
xLength = 1

# audio sample rate:
fs = 44100.0

# sample period
T = 1.0 / fs

# audio sample bit width
bWidth = 24

# video frames per second:
framesPerSec = 30.0

bpm = 133.0

# time signature: 0 = 4/4; 1 = 3/4
timeSig = 0


print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Instantiate clock & signal Generator objects::---*')
print('// *--------------------------------------------------------------* //')

tbClocks = clks.odmkClocks(xLength, fs, bpm, framesPerSec)

numSamples = tbClocks.totalSamples

tbSigGen = sigGen.odmkSigGen1(numSamples, fs)

# // *---------------------------------------------------------------------* //

tbQtrBar = tbClocks.clkQtrBar()
tbQtrBar5 = tbClocks.clkQtrBar(nBar=5)


print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Generate source waveforms::---*')
print('// *--------------------------------------------------------------* //')


print('\n::Mono Sine waves::')
print('generated mono sin signals @ 2.5K and 5K Hz')
# generate simple mono sin waves
testFreq1 = 2500.0
testFreq2 = 5000.0

sin2_5K = tbSigGen.monosin(testFreq1)
sin5K = tbSigGen.monosin(testFreq2)

# // *---------------------------------------------------------------------* //

print('\n::Multi Sine source::')
print('generated array of sin signals "sinArray"')
testFreqs = [666.0, 777.7, 2300.0, 6000.0, 15600.0]
numFreq = len(testFreqs)

print('Frequency Array (Hz):')
print(testFreqs)

sinArray = np.array([])
for freq in testFreqs:
    sinArray = np.concatenate((sinArray, tbSigGen.monosin(freq)))
sinArray = sinArray.reshape((numFreq, numSamples))
sinArray = sinArray.transpose()

# // *---------------------------------------------------------------------* //

# generate a set of orthogonal frequencies

print('\n::Orthogonal Multi Sine source::')
print('generated array of orthogonal sin signals "orthoSinArray"')

# for n freqs, use 2n+1 => skip DC and negative freqs!
# ex. for cyclicZn(15), we want to use czn[1, 2, 3, ... 7]

numOrthoFreq = 7
czn = cyclicZn(2*numOrthoFreq + 1)

orthoFreqArray = np.array([])
for c in range(1, numOrthoFreq+1):
    cznph = np.arctan2(czn[c].imag, czn[c].real)
    cznFreq = (fs*cznph)/(2*np.pi)
    orthoFreqArray = np.append(orthoFreqArray, cznFreq)

print('Orthogonal Frequency Array (Hz):')
print(orthoFreqArray)

orthoSinArray = np.array([])
for freq in orthoFreqArray:
    orthoSinArray = np.concatenate((orthoSinArray, tbSigGen.monosin(freq)))
orthoSinArray = orthoSinArray.reshape((numOrthoFreq, numSamples))
orthoSinArray = orthoSinArray.transpose()


# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Check source waveform spectrum::---*')
print('// *--------------------------------------------------------------* //')

# FFT length
N = 2048

# // *---------------------------------------------------------------------* //

y1 = sin2_5K[0:N]
y2 = sin5K[0:N]

# forward FFT
y1_FFT = sp.fft(y1)
y1_Mag = np.abs(y1_FFT)
y1_Phase = np.arctan2(y1_FFT.imag, y1_FFT.real)
# scale and format FFT out for plotting
y1_FFTscale = 2.0/N * np.abs(y1_FFT[0:N/2])

y2_FFT = sp.fft(y2)
y2_Mag = np.abs(y2_FFT)
y2_Phase = np.arctan2(y2_FFT.imag, y2_FFT.real)
# scale and format FFT out for plotting
y2_FFTscale = 2.0/N * np.abs(y2_FFT[0:N/2])

# inverse FFT
y1_IFFT = sp.ifft(y1_FFT)

y2_IFFT = sp.ifft(y2_FFT)

# check
yDiff = y2_IFFT - y2

# // *---------------------------------------------------------------------* //

# ::sinArray::

yArray = np.array([])
yMagArray = np.array([])
yPhaseArray = np.array([])
yScaleArray = np.array([])
# for h in range(len(sinArray[0, :])):
for h in range(numFreq):    
    yFFT = sp.fft(sinArray[0:N, h])
    yArray = np.concatenate((yArray, yFFT))
    yScaleArray = np.concatenate((yScaleArray, 2.0/N * np.abs(yFFT[0:N/2])))
#    yMagArray = np.concatenate((yMagArray, np.abs(yFFT)))    
#    yPhaseArray = np.concatenate((yPhaseArray, np.arctan2(yFFT.imag, yFFT.real)))

yArray = yArray.reshape((numFreq, N))
yArray = yArray.transpose()

yScaleArray = yScaleArray.reshape((numFreq, N/2))
yScaleArray = yScaleArray.transpose()

# yMagArray = yMagArray.reshape((numFreqs, N))
# yMagArray = yMagArray.transpose()

# yPhaseArray = yPhaseArray.reshape((numFreqs, N))
# yPhaseArray = yPhaseArray.transpose()

# // *---------------------------------------------------------------------* //

# ::orthoSinArray::

yOrthoArray = np.array([])
yOrthoMagArray = np.array([])
yOrthoPhaseArray = np.array([])
yOrthoScaleArray = np.array([])
# for h in range(len(sinArray[0, :])):
for h in range(numOrthoFreq):
    yOrthoFFT = sp.fft(orthoSinArray[0:N, h])
    yOrthoArray = np.concatenate((yOrthoArray, yOrthoFFT))
    yOrthoScaleArray = np.concatenate((yOrthoScaleArray, 2.0/N * np.abs(yOrthoFFT[0:N/2])))

yOrthoArray = yOrthoArray.reshape((numOrthoFreq, N))
yOrthoArray = yOrthoArray.transpose()

yOrthoScaleArray = yOrthoScaleArray.reshape((numOrthoFreq, N/2))
yOrthoScaleArray = yOrthoScaleArray.transpose()

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Phase Vocoder Begin::---*')
print('// *--------------------------------------------------------------* //')

sigLength = 4096

NFFT = 64

Ra = 32

PVAin = y2[0:sigLength]

Awin = np.blackman(NFFT)


PVAout = odmkPVA(PVAin, Awin, NFFT, Ra, fs)





# // *---------------------------------------------------------------------* //

# // *---------------------------------------------------------------------* //
# // *---Mono FFT plots---*
# // *---------------------------------------------------------------------* //

# define a sub-range for wave plot visibility
tLen = 50

fnum = 1
pltTitle = 'Input Signal y1 (first '+str(tLen)+' samples)'
pltXlabel = 'y1: '+str(testFreq1)+' Hz'
pltYlabel = 'Magnitude'


sig = y1[0:tLen]
# define a linear space from 0 to 1/2 Fs for x-axis:
xaxis = np.linspace(0, tLen, tLen)


odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)


fnum = 2
pltTitle = 'Scipy FFT Mag: y1_FFTscale '+str(testFreq1)+' Hz'
pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
pltYlabel = 'Magnitude (scaled by 2/N)'

# sig <= direct

# define a linear space from 0 to 1/2 Fs for x-axis:
xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)

odmkPlot1D(fnum, y1_FFTscale, xfnyq, pltTitle, pltXlabel, pltYlabel)

# // *---------------------------------------------------------------------* //
# // *---plot a single sin out of array---*
# // *---------------------------------------------------------------------* //

#nn = 0
#
## FFT length
#N = 2048
#
#sinAtst = sinArray[0:N, nn]
#sinAtst_frq = testFreqs[nn]
#
## forward FFT
#sinAtst_FFT = sp.fft(sinAtst)
#sinAtst_Mag = np.abs(sinAtst_FFT)
#sinAtst_Phase = np.arctan2(sinAtst_FFT.imag, sinAtst_FFT.real)
## scale and format FFT out for plotting
#sinAtst_FFTscale = 2.0/N * np.abs(sinAtst_FFT[0:N/2])
#
## inverse FFT
#sinAtst_IFFT = sp.ifft(sinAtst_FFT)
#
#fnum = 300
#pltTitle = 'Input Signal sinAtst (first '+str(tLen)+' samples)'
#pltXlabel = 'sinAtst: '+str(sinAtst_frq)+' Hz'
#pltYlabel = 'Magnitude'
#
#sig = sinAtst[0:tLen]
#
## define a linear space from 0 to 1/2 Fs for x-axis:
#xaxis = np.linspace(0, tLen, tLen)
#
#odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)
#
#
#fnum = 301
#pltTitle = 'Scipy FFT Mag: sinAtst '+str(sinAtst_frq)+' Hz'
#pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
#pltYlabel = 'Magnitude (scaled by 2/N)'
#
## define a linear space from 0 to 1/2 Fs for x-axis:
#xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
#
#odmkPlot1D(fnum, sinAtst_FFTscale, xfnyq, pltTitle, pltXlabel, pltYlabel)


# // *---------------------------------------------------------------------* //
# // *---Multi Plot - source signal array vs. FFT MAG out array---*
# // *---------------------------------------------------------------------* //

fnum = 3
pltTitle = 'Input Signals: sinArray (first '+str(tLen)+' samples)'
pltXlabel = 'sinArray time-domain wav'
pltYlabel = 'Magnitude'

# define a linear space from 0 to 1/2 Fs for x-axis:
xaxis = np.linspace(0, tLen, tLen)

odmkMultiPlot1D(fnum, sinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='gnuplot')


fnum = 4
pltTitle = 'FFT Mag: yScaleArray multi-osc '
pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
pltYlabel = 'Magnitude (scaled by 2/N)'

# define a linear space from 0 to 1/2 Fs for x-axis:
xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)

odmkMultiPlot1D(fnum, yScaleArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='gnuplot')


# // *---------------------------------------------------------------------* //
# // *---Orthogonal Sine Plot - source signal array vs. FFT MAG out array---*
# // *---------------------------------------------------------------------* //

fnum = 5
pltTitle = 'Input Signals: orthoSinArray (first '+str(tLen)+' samples)'
pltXlabel = 'orthoSinArray time-domain wav'
pltYlabel = 'Magnitude'

# define a linear space from 0 to 1/2 Fs for x-axis:
xaxis = np.linspace(0, tLen, tLen)

odmkMultiPlot1D(fnum, orthoSinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='gnuplot')


fnum = 6
pltTitle = 'FFT Mag: yOrthoScaleArray multi-osc '
pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
pltYlabel = 'Magnitude (scaled by 2/N)'

# define a linear space from 0 to 1/2 Fs for x-axis:
xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)

odmkMultiPlot1D(fnum, yOrthoScaleArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')

# // *---------------------------------------------------------------------* //

plt.show()

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::done::---*')
print('// *--------------------------------------------------------------* //')

# // *---------------------------------------------------------------------* //
