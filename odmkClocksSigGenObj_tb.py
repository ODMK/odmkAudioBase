# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkClocksSigGenObj_tb.py))::__
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
print('// *---::FFT processing::---*')
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

odmkMultiPlot1D(fnum, sinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='gist_stern')


fnum = 4
pltTitle = 'FFT Mag: yScaleArray multi-osc '
pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
pltYlabel = 'Magnitude (scaled by 2/N)'

# define a linear space from 0 to 1/2 Fs for x-axis:
xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)

odmkMultiPlot1D(fnum, yScaleArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='gist_stern')


# // *---------------------------------------------------------------------* //
# // *---Orthogonal Sine Plot - source signal array vs. FFT MAG out array---*
# // *---------------------------------------------------------------------* //

fnum = 5
pltTitle = 'Input Signals: orthoSinArray (first '+str(tLen)+' samples)'
pltXlabel = 'orthoSinArray time-domain wav'
pltYlabel = 'Magnitude'

# define a linear space from 0 to 1/2 Fs for x-axis:
xaxis = np.linspace(0, tLen, tLen)

odmkMultiPlot1D(fnum, orthoSinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')


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


#tLen = 200
#
## Input signal
#fig1 = plt.figure(num=1, facecolor='silver', edgecolor='k')
#odmkSrcplt1 = plt.plot(y[0:tLen])
#plt.setp(odmkSrcplt1, color='red', ls='-', linewidth=1.00)
#plt.xlabel('monosin5K: '+str(testFreq)+' Hz')
#plt.ylabel('Magnitude')
#plt.title('Input Signal (first '+str(tLen)+' samples)')
#plt.grid(color='c', linestyle=':', linewidth=.5)
#plt.grid(True)
## plt.xticks(np.linspace(0, Fs/2, 10))
#ax = plt.gca()
#ax.set_axis_bgcolor('black')
#
## define a linear space from 0 to 1/2 Fs for x-axis:
#xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
#
## FFT Magnitude out plot (0-fs/2)
#fig2 = plt.figure(num=2, facecolor='silver', edgecolor='k')
#odmkFFTplt1 = plt.plot(xfnyq, yfscale)
#plt.setp(odmkFFTplt1, color='red', ls='-', linewidth=1.00)
#plt.xlabel('Frequency: 0 - '+str(fs / 2)+' Hz')
#plt.ylabel('Magnitude (scaled by 2/N)')
#plt.title('Scipy FFT: Fs = '+str(fs)+' N = '+str(N))
#plt.grid(color='c', linestyle=':', linewidth=.5)
#plt.grid(True)
## plt.xticks(np.linspace(0, Fs/2, 10))
#ax = plt.gca()
#ax.set_axis_bgcolor('black')
