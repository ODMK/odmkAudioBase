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
plt.close('all')
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


# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : Primary parameters definitions
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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


# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : object definition
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Instantiate clock & signal Generator objects::---*')
print('// *--------------------------------------------------------------* //')

tbClocks = clks.odmkClocks(xLength, fs, bpm, framesPerSec)

numSamples = tbClocks.totalSamples

tbSigGen = sigGen.odmkSigGen1(numSamples, fs)

# // *---------------------------------------------------------------------* //

tbclkDownBeats = tbClocks.clkDownBeats()


# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : Generate source waveforms
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


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

odmkTestFreqArray5_1 = [666.0, 777.7, 2300.0, 6000.0, 15600.0]
odmkTestFreqArray5_2 = [444.0, 1776.0]
odmkTestFreqArray5_3 = [3200.0, 6400.0, 9600.0, 12800.0, 16000.0, 19200.0,  22400.0]
multiSinV1 = tbSigGen.multiSin(odmkTestFreqArray5_2)

numOrtFreqs = 7
nCzn = cyclicZn(2*numOrtFreqs + 1)

nOrthogonalArray = np.array([])
for c in range(1, numOrtFreqs+1):
    nCznPh = np.arctan2(nCzn[c].imag, nCzn[c].real)
    nOrthogonalArray = (fs*nCznPh)/(2*np.pi)


# // *---------------------------------------------------------------------* //

print('\n::Multi Sine source::')
print('generated array of sin signals "sinArray"')
# testFreqs = [666.0, 777.7, 2300.0, 6000.0, 15600.0]
testFreqs = odmkTestFreqArray5_1
numFreqSinArray = len(testFreqs)

print('Frequency Array (Hz):')
print(testFreqs)

sinArray = np.array([])
for freq in testFreqs:
    sinArray = np.concatenate((sinArray, tbSigGen.monosin(freq)))
sinArray = sinArray.reshape((numFreqSinArray, numSamples))

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

# generate a composite signal of an array of frequencies

print('\n::Orthogonal Composite Multi Sine source::')
print('generated a Composite array of sin signals "orthoSinComp1"')

# for n freqs, use 2n+1 => skip DC and negative freqs!
# ex. for cyclicZn(15), we want to use czn[1, 2, 3, ... 7]

compFreqArray = odmkTestFreqArray5_3

sinComp1 = tbSigGen.multiSin(compFreqArray)

print('Generated composite sine signal: sinComp1')


# // *---------------------------------------------------------------------* //

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : FFT processing
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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

## ::sinArray Probes::
#ttLen = 200
#sig = sinArray[3, 0:ttLen]
## define a linear space from 0 to 1/2 Fs for x-axis:
#xaxis = np.linspace(0, ttLen, ttLen)
#
#fnum = 300
#pltTitle = 'SigGen output: sinArray[0] (first '+str(ttLen)+' samples)'
#pltXlabel = 'time'
#pltYlabel = 'Magnitude (scaled by ???)'
#
#odmkplt.odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)
#
#
#yProbeArray1 = sp.fft(sinArray[3, 0:N])
#yProbeScale1 = 2.0/N * np.abs(yProbeArray1[0:int(N/2)])
## define a linear space from 0 to 1/2 Fs for x-axis:
#xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
#
#fnum = 301
#pltTitle = 'FFT Mag: y3tri_FFTscale '+str(testFreq1)+' Hz'
#pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
#pltYlabel = 'Magnitude (scaled by 2/N)'
#
#odmkplt.odmkPlot1D(fnum, yProbeScale1, xfnyq, pltTitle, pltXlabel, pltYlabel)


# ::sinArray::
numFreqs = numFreqSinArray    # defined above for gen of sinArray

yArray = np.array([])
yMagArray = np.array([])
yPhaseArray = np.array([])
yScaleArray = np.array([])
# for h in range(len(sinArray[0, :])):
for h in range(numFreqs):    
    yFFT = sp.fft(sinArray[h, 0:N])
    yArray = np.concatenate((yArray, yFFT))
    yScaleArray = np.concatenate((yScaleArray, 2.0/N * np.abs(yFFT[0:int(N/2)])))
yArray = yArray.reshape((numFreqs, N))
yScaleArray = yScaleArray.reshape((numFreqs, int(N/2)))

# yMagArray = yMagArray.reshape((numFreqs, N))
# yMagArray = yMagArray.transpose()

# yPhaseArray = yPhaseArray.reshape((numFreqs, N))
# yPhaseArray = yPhaseArray.transpose()

# // *---------------------------------------------------------------------* //

# ::orthoSinArray::
numFreqs = numOrthoFreq

yOrthoArray = np.array([])
yOrthoMagArray = np.array([])
yOrthoPhaseArray = np.array([])
yOrthoScaleArray = np.array([])
# for h in range(len(sinArray[0, :])):
for h in range(numFreqs):
    yOrthoFFT = sp.fft(orthoSinArray[h, 0:N])
    yOrthoArray = np.concatenate((yOrthoArray, yOrthoFFT))
    yOrthoScaleArray = np.concatenate((yOrthoScaleArray, 2.0/N * np.abs(yOrthoFFT[0:int(N/2)])))
yOrthoArray = yOrthoArray.reshape((numFreqs, N))
yOrthoScaleArray = yOrthoScaleArray.reshape(numFreqs, (int(N/2)))


# // *---------------------------------------------------------------------* //

# ::compositeSinOut1::
ySinComp1 = sinComp1[0:N]

# forward FFT
sinComp1_FFT = sp.fft(ySinComp1)
sinComp1_Mag = np.abs(sinComp1_FFT)
sinComp1_Phase = np.arctan2(sinComp1_FFT.imag, sinComp1_FFT.real)
# scale and format FFT out for plotting
sinComp1_FFTscale = 2.0/N * np.abs(sinComp1_FFT[0:int(N/2)])

# inverse FFT
sinComp1_IFFT = sp.ifft(sinComp1_FFT)


# // *---------------------------------------------------------------------* //

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : Plotting
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# // *---------------------------------------------------------------------* 

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Plotting::---*')
print('// *--------------------------------------------------------------* //')


# // *---------------------------------------------------------------------* //

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


    # // *---------------------------------------------------------------------* //
    # // *---Composite Sine Plot - source signal array vs. FFT MAG out array---*
    # // *---------------------------------------------------------------------* //

   # // *-----------------------------------------------------------------* //


    # odmkTestFreqArray5_3
    tLen = 200
    
    sig = sinComp1[0:tLen]
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    

    fnum = 9
    pltTitle = 'SigGen output: sinComp1 Composite waveform (first '+str(tLen)+' samples)'
    pltXlabel = 'time'
    pltYlabel = 'Magnitude (scaled by ???)'
    
    odmkplt.odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)    

    # // *-----------------------------------------------------------------* //

    fnum = 10
    pltTitle = 'FFT Mag: sinComp1_FFTscale Composite waveform'
    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
    pltYlabel = 'Magnitude (scaled by 2/N)'
    
    # sig <= direct
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    odmkplt.odmkPlot1D(fnum, sinComp1_FFTscale, xfnyq, pltTitle, pltXlabel, pltYlabel)    

    # // *-----------------------------------------------------------------* //




    plt.show()

else:    # comment-off/on: toggle plots below  
    print('\n')
    print('// *---::Plotting Bypassed::---*')


# // *---------------------------------------------------------------------* //

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : Write Signal to .txt File
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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
    # write a 1D sine signal to .txt
    sigIn = sinComp1
    
    outNmTXT = 'sinComp1.txt'
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
