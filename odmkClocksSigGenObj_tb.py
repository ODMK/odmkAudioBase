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

N = tbClocks.totalSamples

tbSigGen = sigGen.odmkSigGen1(N, fs)

# // *---------------------------------------------------------------------* //

tbQtrBar = tbClocks.clkQtrBar()
tbQtrBar5 = tbClocks.clkQtrBar(nBar=5)


print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Instantiate clock & signal Generator objects::---*')
print('// *--------------------------------------------------------------* //')

testFreq = 5000.0

monosin5K = tbSigGen.monosin(testFreq)

y = monosin5K

# forward FFT
yf = sp.fft(y)
yfMag = np.abs(yf)
yfPhase = np.arctan2(yf.imag, yf.real)

# inverse FFT
yInv = sp.ifft(yf)

yDiff = yInv - y


# scale and format FFT out for plotting
yfscale = 2.0/N * np.abs(yf[0:N/2])


# // *---------------------------------------------------------------------* //
# *---plots---*
tLen = 200

# Input signal
fig1 = plt.figure(num=1, facecolor='silver', edgecolor='k')
odmkSrcplt1 = plt.plot(y[0:tLen])
plt.setp(odmkSrcplt1, color='red', ls='-', linewidth=1.00)
plt.xlabel('monosin5K: '+str(testFreq)+' Hz')
plt.ylabel('Magnitude')
plt.title('Input Signal (first '+str(tLen)+' samples)')
plt.grid(color='c', linestyle=':', linewidth=.5)
plt.grid(True)
# plt.xticks(np.linspace(0, Fs/2, 10))
ax = plt.gca()
ax.set_axis_bgcolor('black')

# define a linear space from 0 to 1/2 Fs for x-axis:
xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)

# FFT Magnitude out plot (0-fs/2)
fig2 = plt.figure(num=2, facecolor='silver', edgecolor='k')
odmkFFTplt1 = plt.plot(xfnyq, yfscale)
plt.setp(odmkFFTplt1, color='red', ls='-', linewidth=1.00)
plt.xlabel('Frequency: 0 - '+str(fs / 2)+' Hz')
plt.ylabel('Magnitude (scaled by 2/N)')
plt.title('Scipy FFT: Fs = '+str(fs)+' N = '+str(N))
plt.grid(color='c', linestyle=':', linewidth=.5)
plt.grid(True)
# plt.xticks(np.linspace(0, Fs/2, 10))
ax = plt.gca()
ax.set_axis_bgcolor('black')


#def odmkMatPlot1D(fnum, sig, xLin, pltTitle, pltXlabel, pltYlabel, lncolor='red', lnstyle='-', lnwidth=1.00, pltGrid=True, pltBgColor='black'):
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

