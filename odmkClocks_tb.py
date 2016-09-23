# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkClocks_tb.py))::__
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
    
    
# // *---------------------------------------------------------------------* //
# // *--Plot Functions--*
# // *---------------------------------------------------------------------* //    

#def odmkPlot1D(fnum, sig, xLin, pltTitle, pltXlabel, pltYlabel, lncolor='red', lnstyle='-', lnwidth=1.00, pltGrid=True, pltBgColor='black'):
#    ''' ODMK 1D Matplotlib plot
#        required inputs:
#        fnum => unique plot number
#        sig => signal to plot
#        xLin => linear space to define x-axis (0 to max x-axis length-1)
#        pltTitle => text string for plot title
#        pltXlabel => text string for x-axis
#        pltYlabel => text string for y-axis
#        optional inputs:
#        lncolor => line color (default = red ; html color names, html color codes??)
#        lnstyle => line style (default = plain line ; * ; o ; etc..)
#        lnwidth => line width
#        pltGrid => use grid : default = True ; <True;False>
#        pltBgColor => backgroud color (default = black) '''
#
#    # Input signal
#    plt.figure(num=fnum, facecolor='silver', edgecolor='k')
#    # check if xLin is < than or = to sig
#    if len(xLin) > len(sig):
#        print('ERROR: length of xLin x-axis longer than signal length')
#        return 1
#    elif len(xLin) == len(sig):
#        odmkMatPlt = plt.plot(xLin, sig)
#    else:
#        odmkMatPlt = plt.plot(xLin, sig[0:len(xLin)])
#
#    plt.setp(odmkMatPlt, color=lncolor, ls=lnstyle, linewidth=lnwidth)
#    plt.xlabel(pltXlabel)
#    plt.ylabel(pltYlabel)
#    plt.title(pltTitle)
#    plt.grid(color='c', linestyle=':', linewidth=.5)
#    plt.grid(pltGrid)
#    # plt.xticks(np.linspace(0, Fs/2, 10))
#    ax = plt.gca()
#    ax.set_axis_bgcolor(pltBgColor)
#
#    return 0
#    
#def odmkMultiPlot1D(fnum, sigArray, xLin, pltTitle, pltXlabel, pltYlabel, colorMp='gnuplot', lnstyle='-', lnwidth=1.00, pltGrid=True, pltBgColor='black'):
#    ''' ODMK 1D Matplotlib multi-plot
#        required inputs:
#        fnum => unique plot number
#        sig => signal to plot : 2D Numpy array
#        xLin => linear space to define x-axis (0 to max x-axis length-1)
#        pltTitle => text string for plot title
#        pltXlabel => text string for x-axis
#        pltYlabel => text string for y-axis
#        optional inputs:
#        lncolor => line color (default = red ; html color names, html color codes??)
#        lnstyle => line style (default = plain line ; * ; o ; etc..)
#        lnwidth => line width
#        pltGrid => use grid : default = True ; <True;False>
#        pltBgColor => backgroud color (default = black) '''
#
#    # define the color map
#    try:
#        cmap = plt.cm.get_cmap(colorMp)
#    except ValueError as e:
#        print('ValueError: ', e)
#    colors = cmap(np.linspace(0.0, 1.0, len(sigArray[0, :])))
#
#    # Input signal
#    plt.figure(num=fnum, facecolor='silver', edgecolor='k')
#    # check if xLin is < than or = to sig
#    if len(xLin) > len(sigArray[:, 0]):
#        print('ERROR: length of xLin x-axis longer than signal length')
#        return 1
#    else:
#        if len(xLin) == len(sigArray[:, 0]):
#            # odmkMatPlt = []
#            for i in range(len(sinArray[0, :])):
#                plt.plot(xLin, sigArray[:, i], color=colors[i], ls=lnstyle, linewidth=lnwidth)
#        else:
#            # odmkMatPlt = []
#            for i in range(len(sinArray[0, :])):
#                plt.plot(xLin, sigArray[0:len(xLin), i], color=colors[i], ls=lnstyle, linewidth=lnwidth)
#
#        plt.xlabel(pltXlabel)
#        plt.ylabel(pltYlabel)
#        plt.title(pltTitle)
#        plt.grid(color='c', linestyle=':', linewidth=.5)
#        plt.grid(pltGrid)
#        # plt.xticks(np.linspace(0, Fs/2, 10))
#        ax = plt.gca()
#        ax.set_axis_bgcolor(pltBgColor)
#
#    return 0    

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

# audio sample rate:
fs = 48000.0

# sample period
T = 1.0 / fs

# audio sample bit width
bWidth = 24

# video frames per second:
framesPerSec = 30.0

bpm = 133.0

# time signature: 4 = 4/4; 3 = 3/4
timeSig = 4

# calculate length of 1　bar in seconds:
spb = 60.0 / bpm
secondsPerBar = timeSig * spb

# length of x in seconds:
#xLength = 1
xLength = secondsPerBar

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Instantiate clock & signal Generator objects::---*')
print('// *--------------------------------------------------------------* //')

# // *---------------------------------------------------------------------* //

tbClocks = clks.odmkClocks(xLength, fs, bpm, framesPerSec)


tbxLength = tbClocks.xLength
tbFs = tbClocks.fs
tbBpm = tbClocks.bpm
tbframesPerSec = tbClocks.framesPerSec

# *---Audio secondary Parameters---*

# set bps - beats per second
tbBps = tbClocks.bps
# set spb - Seconds per beat
tbSpb = tbClocks.spb
# set samplesPerBeat
tbsamplesPerBeat = tbClocks.samplesPerBeat
# set samples per bar / 1 bar = tsig beats (ex. 4/4: 4*samplesPerBeat)
tbsamplesPerBar = tbClocks.samplesPerBar

# set totalSamples - Total audio samples in x
tbtotalSamples = tbClocks.totalSamples
# set totalBeats - total beats in x
tbtotalBeats = tbClocks.totalBeats


# *---Video secondary Parameters---*

# set samplesPerFrame - Num audio samples per video frame
tbsamplesPerFrame = tbClocks.samplesPerFrame
# set framesPerBeat - Num video frames per beat
tbframesPerBeat = tbClocks.framesPerBeat

# set totalFrames - Total video frames in x
TBtotalFrames = tbClocks.totalFrames


tbclkDownBeats = tbClocks.clkDownBeats()

tbclkDownFrames = tbClocks.clkDownFrames()

tbclkQtrBeat = tbClocks.clkQtrBeat()

nbar = 3
tbclkQtrBeatBar = tbClocks.clkQtrBeatBar(nbar)

n = 7
tbclkDivNBeat = tbClocks.clkDivNBeat(n)


print('\nCreated odmkClocks object "tbClock" with the following parameters:')
print('\nAn odmkClocks object has been instanced with:')
print('xLength = '+str(tbxLength))
print('fs = '+str(tbFs))
print('bpm = '+str(tbBpm))
print('framesPerSec = '+str(tbframesPerSec))
print('beatsPerSecond = '+str(tbBps))
print('secondsPerBeat = '+str(tbSpb))
print('samplesPerBeat = '+str(tbsamplesPerBeat))
print('samplesPerBar = '+str(tbsamplesPerBar))
print('totalSamples = '+str(tbtotalSamples))
print('totalBeats = '+str(tbtotalBeats))
print('samplesPerFrame = '+str(tbsamplesPerFrame))
print('framesPerBeat = '+str(tbframesPerBeat))
print('totalFrames = '+str(TBtotalFrames))


# // *---------------------------------------------------------------------* //

odmkPlots = 0;
if odmkPlots == 1:
    
    print('\n')
    print('// *--------------------------------------------------------------* //')
    print('// *---::Plotting::---*')
    print('// *--------------------------------------------------------------* //')

    # // *---------------------------------------------------------------------* //
    # // *---tbclkDivNBeat---*
    # // *---------------------------------------------------------------------* //
    
    # define a sub-range for wave plot visibility
    tLen = tbtotalSamples
    
    fnum = 1
    pltTitle = 'Input Signal tbclkDivNBeat (first '+str(tLen)+' samples)'
    pltXlabel = 'tbclkDivNBeat: '+str(n)+' beats per bar'
    pltYlabel = 'Magnitude'

    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    
    odmkplt.odmkPlot1D(fnum, tbclkDivNBeat[0:tLen], xaxis, pltTitle, pltXlabel, pltYlabel)
    
    
    
#    # // *---------------------------------------------------------------------* //
#    # // *---Multi Plot - source signal array vs. FFT MAG out array---*
#    # // *---------------------------------------------------------------------* //
#
#    fnum = 3
#    pltTitle = 'Input Signals: sinArray (first '+str(tLen)+' samples)'
#    pltXlabel = 'sinArray time-domain wav'
#    pltYlabel = 'Magnitude'
#    
#    # define a linear space from 0 to 1/2 Fs for x-axis:
#    xaxis = np.linspace(0, tLen, tLen)
#    
#    odmkMultiPlot1D(fnum, sinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='gist_stern')
#    
#    
#    fnum = 4
#    pltTitle = 'FFT Mag: yScaleArray multi-osc '
#    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
#    pltYlabel = 'Magnitude (scaled by 2/N)'
#    
#    # define a linear space from 0 to 1/2 Fs for x-axis:
#    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
#    
#    odmkMultiPlot1D(fnum, yScaleArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='gist_stern')
#    
#    
#    # // *---------------------------------------------------------------------* //
#    # // *---Orthogonal Sine Plot - source signal array vs. FFT MAG out array---*
#    # // *---------------------------------------------------------------------* //
#    
#    fnum = 5
#    pltTitle = 'Input Signals: orthoSinArray (first '+str(tLen)+' samples)'
#    pltXlabel = 'orthoSinArray time-domain wav'
#    pltYlabel = 'Magnitude'
#    
#    # define a linear space from 0 to 1/2 Fs for x-axis:
#    xaxis = np.linspace(0, tLen, tLen)
#    
#    odmkMultiPlot1D(fnum, orthoSinArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')
#    
#    
#    fnum = 6
#    pltTitle = 'FFT Mag: yOrthoScaleArray multi-osc '
#    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
#    pltYlabel = 'Magnitude (scaled by 2/N)'
#    
#    # define a linear space from 0 to 1/2 Fs for x-axis:
#    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
#    
#    odmkMultiPlot1D(fnum, yOrthoScaleArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='hsv')
#    
#    # // *-----------------------------------------------------------------* //
#        
#    
#    # define a sub-range for wave plot visibility
#    tLen = 500
#    
#    fnum = 7
#    pltTitle = 'Input Signal tri2_5K (first '+str(tLen)+' samples)'
#    pltXlabel = 'tri2_5K: '+str(testFreq1)+' Hz'
#    pltYlabel = 'Magnitude'
#    
#    sig = tri2_5K[0:tLen]
#    # define a linear space from 0 to 1/2 Fs for x-axis:
#    xaxis = np.linspace(0, tLen, tLen)
#    
#    odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)
#    
#    # // *-----------------------------------------------------------------* //    
#    
#    fnum = 8
#    pltTitle = 'Scipy FFT Mag: y1_FFTscale '+str(testFreq1)+' Hz'
#    pltXlabel = 'Frequency: 0 - '+str(fs / 2)+' Hz'
#    pltYlabel = 'Magnitude (scaled by 2/N)'
#    
#    # sig <= direct
#    
#    # define a linear space from 0 to 1/2 Fs for x-axis:
#    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
#    
#    odmkPlot1D(fnum, y3tri_FFTscale, xfnyq, pltTitle, pltXlabel, pltYlabel)    
#
#    # // *-----------------------------------------------------------------* //

else:    # comment-off/on: toggle plots below
    print('\n')  
    print('// *---::No Plotting / Debugging::---*')    
    
    plt.show()

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
