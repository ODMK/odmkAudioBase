# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

# __::((odmkWavIO_tb.py))::__

# Created on Sat Aug 20 17:28:14 2016

# 24 bit vectorized .wav processing by the odorousbeast
# this is my house, it's a horrorhouse
#
# Microsoft/RIFF WAVE files use 2's-complement signed integers for PCM samples
# the byte-order of those integers is little-endian.
#
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header end-------------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

import sys
import numpy as np
import matplotlib.pyplot as plt

rootDir = 'C:/odmkDev/odmkCode/odmkPython/'
audioScrDir = 'C:/odmkDev/odmkCode/odmkPython/audio/wavsrc/'
audioOutDir = 'C:/odmkDev/odmkCode/odmkPython/audio/wavout/'


#sys.path.insert(0, 'C:/odmkDev/odmkCode/odmkPython/util')
sys.path.insert(0, rootDir+'util')
from odmkClear import *
#from odmkPlotUtil import *
import odmkPlotUtil as odmkplt

#sys.path.insert(1, 'C:/odmkDev/odmkCode/odmkPython/DSP')
sys.path.insert(1, rootDir+'audio')
import odmkWavIO as waveio

#sys.path.insert(1, 'C:/odmkDev/odmkCode/odmkPython/DSP')
sys.path.insert(2, rootDir+'DSP')
import odmkClocks as clks
import odmkSigGen1 as sigGen

# temp python debugger - use >>>pdb.set_trace() to set break
# import pdb


# // *---------------------------------------------------------------------* //
clear_all()

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


# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Read .wav file & extract parameters::---*')
print('// *--------------------------------------------------------------* //')

tbWaveIO = waveio.odmkWavIO()

wavfile_A = 'bbsineAD.wav'


wavA_stereo, wavA_param = tbWaveIO.wavRead(wavfile_A, audioScrDir)

# wavIn_param = {'fSampleRate': fSampleRate, 'fSampleBits': fSampleBits, 
#                'fChannels': fChannels, 'fNumSamples': fNumSamples}


wavA_sampleRate = wavA_param['fSampleRate']
wavA_sampleBits = wavA_param['fSampleBits']
wavA_channels = wavA_param['fChannels']
wavA_numSamples = wavA_param['fNumSamples']

print('\n::((read .wav file '+wavfile_A+' with the following parameters))::')
print('\nSample Rate = '+str(wavA_sampleRate))
print('Sample Bitwidth = '+str(wavA_sampleBits))
print('Number of Channels = '+str(wavA_channels))
print('Number of Samples = '+str(wavA_numSamples))


# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Clone .wav file & write::---*')
print('// *--------------------------------------------------------------* //')


wavClone_stereo = wavA_stereo
wavClone_param = wavA_param
fs = 44100

wavfile_clone = 'wavclone000.wav'
wavfile_dir = audioOutDir

tbWaveIO.wavWrite(wavClone_stereo, wavfile_clone, wavfile_dir, fs)

print('\n::((wrote .wav file '+wavfile_clone+' with the following parameters))::')
print('\nSample Rate = '+str(wavA_param['fSampleRate']))
print('Sample Bitwidth = '+str(wavA_param['fSampleBits']))
print('Number of Channels = '+str(wavA_param['fChannels']))
print('Number of Samples = '+str(wavA_param['fNumSamples']))


# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Read cloned .wav file & verify against original::---*')
print('// *--------------------------------------------------------------* //')



# Verify

#wavDiff = bytes(len(wavC))
#
#byteErrorCnt = 0
#for j in range(len(wavC)):
#    if wavC[j] != wavA[j]:
#        wavDiff[j] = wavC[j] - wavA[j]
#        byteErrorCnt += 1
#if byteErrorCnt == 0:
#    print('\nSuccessfully cloned: '+wavfile_A)
#else:
#    print('ERROR: cloning unsuccessful, byte error count = '+str(byteErrorCnt))


# // *---------------------------------------------------------------------* //

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

# // *---------------------------------------------------------------------* //


#print('\n')
#print('// *--------------------------------------------------------------* //')
#print('// *---::Load TXT file into OBJ "sinesrc"::---*')
#print('// *--------------------------------------------------------------* //')
#
#sinesrc = u'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\csvsrc\\sintest1.txt'
#
## reads .csv data into Python List:
#datain = []
#with open(sinesrc, mode='r') as infile:
#    for line in infile.readlines():
#        l,dumb = line.strip().split(',')
#        datain.append((l))
#
#print('\nLoaded file: ')
#
#lgth = len(list(datain))    # get length by iterating csvin obj (only way?)
#print('Length of datain = '+str(lgth))


# // *---------------------------------------------------------------------* //


# /////////////////////////////////////////////////////////////////////////////////////////////
# #############################################################################################
# begin : plotting
# #############################################################################################
# /////////////////////////////////////////////////////////////////////////////////////////////

plottt = 0
if plottt == 1:

    # ********* TB CLK and Signal Plotting *********
    
    # // *---------------------------------------------------------------------* //
    # // *---Mono FFT plots---*
    # // *---------------------------------------------------------------------* //
    
    # define a sub-range for wave plot visibility
    
    #tLen = len(wavA_stereo[:,0])
    tLen = len(wavA_stereo[0])    
    
    fnum = 1
    pltTitle = 'Input Signal <<wavA_Left>>'
    pltXlabel = 'time domain signal'
    pltYlabel = 'amplitude'

    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)

    #odmkplt.odmkPlot1D(fnum, wavA_stereo[:,0], xaxis, pltTitle, pltXlabel, pltYlabel)
    odmkplt.odmkPlot1D(fnum, wavA_stereo[0], xaxis, pltTitle, pltXlabel, pltYlabel)


    tLen = 333
    
    fnum = 2
    pltTitle = 'Input Signal <<wavA_Left>> 1st '+str(tLen)+' samples'
    pltXlabel = 'time domain signal'
    pltYlabel = 'amplitude'
    
    # sig <= direct
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    
    odmkplt.odmkPlot1D(fnum, wavA_stereo[0], xaxis, pltTitle, pltXlabel, pltYlabel)


    tLen = len(wavA_stereo[1])    
    
    fnum = 3
    pltTitle = 'Input Signal <<wavA_Right>>'
    pltXlabel = 'time domain signal'
    pltYlabel = 'amplitude'

    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)

    odmkplt.odmkPlot1D(fnum, wavA_stereo[1], xaxis, pltTitle, pltXlabel, pltYlabel)


    tLen = 333
    
    fnum = 4
    pltTitle = 'Input Signal <<wavA_Right>> 1st '+str(tLen)+' samples'
    pltXlabel = 'time domain signal'
    pltYlabel = 'amplitude'
    
    # sig <= direct
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    
    odmkplt.odmkPlot1D(fnum, wavA_stereo[1], xaxis, pltTitle, pltXlabel, pltYlabel)


    plt.show()

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::done::---*')
print('// *--------------------------------------------------------------* //')

# // *---------------------------------------------------------------------* //