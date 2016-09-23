# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkSigGenObj_tb.py))::__
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

rootDir = 'C:/odmkDev/odmkCode/odmkPython/'
audioScrDir = 'C:/odmkDev/odmkCode/odmkPython/audio/wavsrc/'
audioOutDir = 'C:/odmkDev/odmkCode/odmkPython/audio/wavout/'

import sys

#sys.path.insert(0, 'C:/odmkDev/odmkCode/odmkPython/util')
sys.path.insert(0, rootDir+'util')
from odmkClear import *
#from odmkPlotUtil import *
import odmkPlotUtil as odmkplt

#sys.path.insert(1, 'C:/odmkDev/odmkCode/odmkPython/DSP')
sys.path.insert(1, rootDir+'audio')
import native_wav as wavio

#sys.path.insert(1, 'C:/odmkDev/odmkCode/odmkPython/DSP')
sys.path.insert(2, rootDir+'DSP')
import odmkClocks as clks
import odmkSigGen1 as sigGen

# temp python debugger - use >>>pdb.set_trace() to set break
#import pdb

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
fs = 48000.0

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

tbclkDownBeats = tbClocks.clkDownBeats()


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

cos2_5K = tbSigGen.monocos(testFreq1)
cos5K = tbSigGen.monocos(testFreq2)

tri2_5K = tbSigGen.monotri(testFreq1)
tri5K = tbSigGen.monotri(testFreq2)

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

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::FFT processing::---*')
print('// *--------------------------------------------------------------* //')

# FFT length
N = 2048

print('\nPerforming FFT with length = ' + str(N))

# // *---------------------------------------------------------------------* //

y1 = sin2_5K[0:N]
y2 = sin5K[0:N]

# forward FFT
y1_FFT = sp.fft(y1)
y1_Mag = np.abs(y1_FFT)
y1_Phase = np.arctan2(y1_FFT.imag, y1_FFT.real)
# scale and format FFT out for plotting
y1_FFTscale = 2.0/N * np.abs(y1_FFT[0:int(N/2)])

# inverse FFT
y1_IFFT = sp.ifft(y1_FFT)


y2_FFT = sp.fft(y2)
y2_Mag = np.abs(y2_FFT)
y2_Phase = np.arctan2(y2_FFT.imag, y2_FFT.real)
# scale and format FFT out for plotting
y2_FFTscale = 2.0/N * np.abs(y2_FFT[0:int(N/2)])

# inverse FFT
y2_IFFT = sp.ifft(y2_FFT)

# check
yDiff = y2_IFFT - y2

y3tri = tri2_5K[0:N]
y3tri_FFT = sp.fft(y3tri)
y3tri_Mag = np.abs(y3tri_FFT)
y3tri_Phase = np.arctan2(y3tri_FFT.imag, y3tri_FFT.real)
# scale and format FFT out for plotting
y3tri_FFTscale = 2.0/N * np.abs(y3tri_FFT[0:int(N/2)])

# // *---------------------------------------------------------------------* //

# ::sinArray::

yArray = np.array([])
yMagArray = np.array([])
yPhaseArray = np.array([])
yScaleArray = np.array([])
# for h in range(len(sinArray[0, :])):
for h in range(numFreq):    
    yFFT = sp.fft(sinArray[h, 0:N])
    yArray = np.concatenate((yArray, yFFT))
    yScaleArray = np.concatenate((yScaleArray, 2.0/N * np.abs(yFFT[0:int(N/2)])))

yArray = yArray.reshape((numFreq, N))

yScaleArray = yScaleArray.reshape((numFreq, int(N/2)))

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
    yOrthoFFT = sp.fft(orthoSinArray[h, 0:N])
    yOrthoArray = np.concatenate((yOrthoArray, yOrthoFFT))
    yOrthoScaleArray = np.concatenate((yOrthoScaleArray, 2.0/N * np.abs(yOrthoFFT[0:int(N/2)])))

yOrthoArray = yOrthoArray.reshape((numOrthoFreq, N))

yOrthoScaleArray = yOrthoScaleArray.reshape(numOrthoFreq, (int(N/2)))


# // *---------------------------------------------------------------------* //

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Plotting::---*')
print('// *--------------------------------------------------------------* //')


odmkPlots = 0
if odmkPlots == 1:
    

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
    
    
    odmkplt.odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)
    
    
    fnum = 2
    pltTitle = 'FFT Mag: y1_FFTscale '+str(testFreq1)+' Hz'
    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
    pltYlabel = 'Magnitude (scaled by 2/N)'
    
    # sig <= direct
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    odmkplt.odmkPlot1D(fnum, y1_FFTscale, xfnyq, pltTitle, pltXlabel, pltYlabel)
    
    
    
    # // *---------------------------------------------------------------------* //
    # // *---Multi Plot - source signal array vs. FFT MAG out array---*
    # // *---------------------------------------------------------------------* //

    fnum = 3
    pltTitle = 'Input Signals: sinArray (first '+str(tLen)+' samples)'
    pltXlabel = 'sinArray time-domain wav'
    pltYlabel = 'Magnitude'
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    
    odmkplt.odmkMultiPlot1D(fnum, sinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')
    
    
    fnum = 4
    pltTitle = 'FFT Mag: yScaleArray multi-osc '
    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
    pltYlabel = 'Magnitude (scaled by 2/N)'
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    odmkplt.odmkMultiPlot1D(fnum, yScaleArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')
    
    
    # // *---------------------------------------------------------------------* //
    # // *---Orthogonal Sine Plot - source signal array vs. FFT MAG out array---*
    # // *---------------------------------------------------------------------* //
    
    fnum = 5
    pltTitle = 'Input Signals: orthoSinArray (first '+str(tLen)+' samples)'
    pltXlabel = 'orthoSinArray time-domain wav'
    pltYlabel = 'Magnitude'
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    
    odmkplt.odmkMultiPlot1D(fnum, orthoSinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')
    
    
    fnum = 6
    pltTitle = 'FFT Mag: yOrthoScaleArray multi-osc '
    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
    pltYlabel = 'Magnitude (scaled by 2/N)'
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    odmkplt.odmkMultiPlot1D(fnum, yOrthoScaleArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')
   
    # // *-----------------------------------------------------------------* //
        
    
    # define a sub-range for wave plot visibility
    tLen = 500
    
    fnum = 7
    pltTitle = 'Input Signal tri2_5K (first '+str(tLen)+' samples)'
    pltXlabel = 'tri2_5K: '+str(testFreq1)+' Hz'
    pltYlabel = 'Magnitude'
    
    sig = tri2_5K[0:tLen]
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    
    odmkplt.odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)
    
    # // *-----------------------------------------------------------------* //    
    
    fnum = 8
    pltTitle = 'FFT Mag: y3tri_FFTscale '+str(testFreq1)+' Hz'
    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
    pltYlabel = 'Magnitude (scaled by 2/N)'
    
    # sig <= direct
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    odmkplt.odmkPlot1D(fnum, y3tri_FFTscale, xfnyq, pltTitle, pltXlabel, pltYlabel)    

    # // *-----------------------------------------------------------------* //


else:    # comment-off/on: toggle plots below  
    print('\n')
    print('// *---::No Plotting / Debugging::---*')

    plt.show()

# // *---------------------------------------------------------------------* //

odmkTxtFileOut = 1
if odmkTxtFileOut == 1:

    print('\n')
    print('// *--------------------------------------------------------------* //')
    print('// *---::Write Signal to .txt File::---*')
    print('// *--------------------------------------------------------------* //')

    # rootDir = 'C:/odmkDev/odmkCode/odmkPython/'
    outputDir = rootDir+'DSP/werk/'

    # // *-----------------------------------------------    
    # write a 1D sine signal to .txt
    sigIn = sin2_5K
    
    outNmTXT = 'sin2_5K_src.txt'
    nChan = 1
    
    tbSigGen.sig2txt(sigIn, nChan, outNmTXT, outDir=outputDir)
    print('\nwrote data to file: '+outputDir+outNmTXT)


    # // *-----------------------------------------------    
    # write multi-array of 1D sine signals to .txt
    sigIn = yOrthoScaleArray
    
    outNmTXT = 'yOrthoScaleArray.txt'
    nChan = len(yOrthoScaleArray)
    
    tbSigGen.sig2txt(sigIn, nChan, outNmTXT, outDir=outputDir)
    print('\nwrote data to file: '+outputDir+outNmTXT)


# // *---------------------------------------------------------------------* //

odmkCsvFileOut = 0
if odmkCsvFileOut == 1:

    print('\n')
    print('// *--------------------------------------------------------------* //')
    print('// *---::Write Signal to .csv File::---*')
    print('// *--------------------------------------------------------------* //')
    
    sigIn = sin2_5K
    
    outNmCSV = 'sin2_5K_src.csv'
    
    # rootDir = 'C:/odmkDev/odmkCode/odmkPython/'
    outputDir = rootDir+'DSP/werk/'
    
    tbSigGen.sig2csv(sigIn, outNmCSV, outDir=outputDir)


# // *---------------------------------------------------------------------* //


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
