# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkFFT_ftxt.py))::__
#
# Python .wav file read / write
# CSV read/write
# FFT <-> iFFT verification
# Basic Spectral Plotting
#
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header end-------------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

import os
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


rootDir = 'C:\\odmkDev\\odmkCode\\odmkPython\\'
signalSrcDir = 'C:\\odmkDev\\odmkCode\\odmkPython\\DSP\\test\\'


#sys.path.insert(0, 'C:/odmkDev/odmkCode/odmkPython/util')
sys.path.insert(0, rootDir+'util')
from odmkClear import *
import odmkPlotUtil as odmkplt

sys.path.insert(1, rootDir+'DSP')
import odmkClocks as clks
import odmkWavGen1 as wavGen


# temp python debugger - use >>>pdb.set_trace() to set break
import pdb

# // *---------------------------------------------------------------------* //
plt.close('all')
clear_all()

# // *---------------------------------------------------------------------* //

print('// //////////////////////////////////////////////////////////////// //')
print('// *--------------------------------------------------------------* //')
print('// *---::ODMK Fourier Transforms::---*')
print('// *--------------------------------------------------------------* //')
print('// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ //')


# // *---------------------------------------------------------------------* //

#print('\n')
#print('// *--------------------------------------------------------------* //')
#print('// *---::Basic Python/Scipy FFT test::---*')
#print('// *--------------------------------------------------------------* //')
#
#
## Sampling Freq
#Fs = 48000.0
## Number of sample points
#N = 4096
## sample period
#T = 1.0 / Fs
#
#freqs = [5000.0, 7000.0]
#
## create composite sin source
#x = np.linspace(0.0, N*T, N)
#y = np.sin(freqs[0] * 2.0*np.pi*x) + 0.5*np.sin(freqs[1] * 2.0*np.pi*x)
#
#print('\nCreated mixed sinusoid source signal "y"')
#
## forward FFT
#yf = sp.fft(y)
#yfMag = np.abs(yf)
#yfPhase = np.arctan2(yf.imag, yf.real)
#
## inverse FFT
#yInv = sp.ifft(yf)
#
#yDiff = yInv - y
#
#
## scale and format FFT out for plotting
#yfscale = 2.0/N * np.abs(yf[0:int(N/2)])
#
#print('\nPerformed FFT, calculate Mag and Phase, create scaled signal "yfscale"')
#
#
## // *---------------------------------------------------------------------* //
## *---plots---*
#
#tLen = 200
#
## Input signal
#fig1 = plt.figure(num=1, facecolor='silver', edgecolor='k')
#odmkSrcplt1 = plt.plot(y[0:tLen])
#plt.setp(odmkSrcplt1, color='red', ls='-', linewidth=1.00)
#plt.xlabel('mixed sine: '+str(freqs[0])+' + '+str(freqs[1])+' Hz')
#plt.ylabel('Magnitude')
#plt.title('Input Signal (first '+str(tLen)+' samples)')
#plt.grid(color='c', linestyle=':', linewidth=.5)
#plt.grid(True)
## plt.xticks(np.linspace(0, Fs/2, 10))
#ax = plt.gca()
#ax.set_axis_bgcolor('black')
#
#
#
## define a linear space from 0 to 1/2 Fs for x-axis:
#xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
#
#
## FFT Magnitude out plot (0-Fs/2)
#fig2 = plt.figure(num=2, facecolor='silver', edgecolor='k')
#odmkFFTplt1 = plt.plot(xfnyq, yfscale)
#plt.setp(odmkFFTplt1, color='red', ls='-', linewidth=1.00)
#plt.xlabel('Frequency: 0 - '+str(Fs / 2)+' Hz')
#plt.ylabel('Magnitude (scaled by 2/N)')
#plt.title('Scipy FFT: Fs = '+str(Fs)+' N = '+str(N))
#plt.grid(color='c', linestyle=':', linewidth=.5)
#plt.grid(True)
## plt.xticks(np.linspace(0, Fs/2, 10))
#ax = plt.gca()
#ax.set_axis_bgcolor('black')
#
#
## Output signal vs. Input signal
#fig3 = plt.figure(num=3, facecolor='silver', edgecolor='k')
#odmkSrcplt2 = plt.plot(y[0:tLen])
#plt.setp(odmkSrcplt2, color='red', ls='-', linewidth=1.00)
#odmkResplt1 = plt.plot(yInv[0:tLen])
#plt.setp(odmkResplt1, color='orange', ls='-', linewidth=1.00)
#plt.xlabel('mixed sine: '+str(freqs[0])+' + '+str(freqs[1])+' Hz')
#plt.ylabel('Magnitude')
#plt.title('Input (red) VS Output (orange) (first '+str(tLen)+' samples)')
#plt.grid(color='c', linestyle=':', linewidth=.5)
#plt.grid(True)
## plt.xticks(np.linspace(0, Fs/2, 10))
#ax = plt.gca()
#ax.set_axis_bgcolor('black')
#
#
## Output vs. Input difference
#fig4 = plt.figure(num=4, facecolor='silver', edgecolor='k')
#odmkDiffplt = plt.plot(yDiff)
#plt.setp(odmkDiffplt, color='orange', ls='-', linewidth=1.00)
#plt.xlabel('time domain')
#plt.ylabel('Magnitude')
#plt.title('Output VS Input Difference')
#plt.grid(color='c', linestyle=':', linewidth=.5)
#plt.grid(True)
## plt.xticks(np.linspace(0, Fs/2, 10))
#ax = plt.gca()
#ax.set_axis_bgcolor('black')

# // *---------------------------------------------------------------------* //


# // *---------------------------------------------------------------------* //

if 0:
    
    #Simple single mono-signal time & frequency plots

    print('\n')
    print('// *--------------------------------------------------------------* //')
    print('// *---:: Mono Signal - time domain plot, freq domain plot ::---*')
    print('// *--------------------------------------------------------------* //')
    
    
    # Sampling Freq (must agree with OSC source)
    Fs = 100000000.0
    #Fs = 48000.0
    
    # Number of sample points
    N = 10000
    # sample period
    T = 1.0 / Fs
    


    signalsrc = signalSrcDir+'osc4T_sin_out.txt'
    #signalsrc = signalSrcDir+'osc4T_saw_out.txt'
    #signalsrc = signalSrcDir+'osc4T_sqr_out.txt'
    #signalsrc = signalSrcDir+'osc4T_pwm_out.txt'
    #signalsrc = signalSrcDir+'osc4T_lfo_out.txt'
    
    # reads .csv data into Numpy array:
    datalist = []
    with open(signalsrc, mode='r') as infile:
        for line in infile.readlines():
            datalist.append(float(line))
    datain_txt = np.array(datalist)
    
    src_name = os.path.split(signalsrc)[1]
    # src_path = os.path.split(sinesrc)[0]
    
    print('\nLoaded file: '+src_name)
    
    lgth = len(list(datain_txt))    # get length by iterating csvin obj (only way?)
    print('Length of datain = '+str(lgth))
    
        
    # forward FFT
    datain_txt_fft = sp.fft(datain_txt[0:N])
    datain_txt_fmag = np.abs(datain_txt_fft)
    datain_txt_phase = np.arctan2(datain_txt_fft.imag, datain_txt_fft.real)
    
    # inverse FFT
    datain_txt_inv = sp.ifft(datain_txt_fft)
    
    datain_txt_diff = datain_txt_inv - datain_txt[0:N]
    
    
    # scale and format FFT out for plotting
    #datain_txt_scale = 2.0/N * np.abs(datain_txt_fft[0:int(N/2)])
    datain_txt_scale = 2.0/N * np.abs(datain_txt_fft[0:int(N/2)])
    
    
    print('\nPerformed FFT, calculate Mag and Phase, create scaled signal "datain_txt_scale"')
        
    
    # // *---------------------------------------------------------------------* //
    # *---plots---*
    
    #tLen = 16384
    tLen = 48000
    
    # Input signal
    fig101 = plt.figure(num=101, facecolor='silver', edgecolor='k')
    odmkSrcplt1 = plt.plot(datain_txt[0:tLen])     # scale for viewing peaks
    plt.setp(odmkSrcplt1, color='red', ls='-', linewidth=1.00)
    plt.xlabel('time domain signal (1st '+str(tLen)+' samples)')
    plt.ylabel('Magnitude')
    plt.title('test DATAIN (first '+str(tLen)+' samples)')
    plt.grid(color='c', linestyle=':', linewidth=.5)
    plt.grid(True)
    # plt.xticks(np.linspace(0, Fs/2, 10))
    ax = plt.gca()
    ax.set_axis_bgcolor('black')
    
    
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    
    # FFT Magnitude out plot (0-Fs/2)
    fig102 = plt.figure(num=102, facecolor='silver', edgecolor='k')
    odmkFFTplt1 = plt.plot(xfnyq, datain_txt_scale)
    plt.setp(odmkFFTplt1, color='red', ls='-', linewidth=1.00)
    plt.xlabel('Frequency: 0 - '+str(Fs / 2)+' Hz')
    plt.ylabel('Magnitude (scaled by 2/N)')
    plt.title('test DATAIN - FFT: Fs = '+str(Fs)+', NFFT = '+str(N))
    plt.grid(color='c', linestyle=':', linewidth=.5)
    plt.grid(True)
    # plt.xticks(np.linspace(0, Fs/2, 10))
    ax = plt.gca()
    ax.set_axis_bgcolor('black')


# // *---------------------------------------------------------------------* //


# // *---------------------------------------------------------------------* //

if 1:

    print('\n')
    print('// *--------------------------------------------------------------* //')
    print('// *---:: Mono Signal - time domain plot, freq domain plot ::---*')
    print('// *--------------------------------------------------------------* //')
    
    
    # Sampling Freq (must agree with OSC source)
    #Fs = 100000000.0
    Fs = 48000.0
    
    # Number of sample points
    N = 10000
    # sample period
    T = 1.0 / Fs
    
    
    #signalSrc = signalSrcDir+'moogHL_refL_out.dat'
    signalSrc = signalSrcDir+'out_dly_L.dat'


    signalSrc1 = signalSrcDir+'osc4T_sin_out.txt'
    signalSrc2 = signalSrcDir+'osc4T_saw_out.txt'
    signalSrc3 = signalSrcDir+'osc4T_sqr_out.txt'
    signalSrc4 = signalSrcDir+'osc4T_pwm_out.txt'
    signalSrc5 = signalSrcDir+'osc4T_lfo_out.txt'
    
    signalSrcArray = [signalSrc1, signalSrc2, signalSrc3, signalSrc4, signalSrc5]
    
    
    # // *---------------------------------------------------------------------* //
    # *---read in data from txt file---*
    
    # reads .csv data into Numpy array:
    datalist = []
    with open(signalSrc, mode='r') as infile:
        for line in infile.readlines():
            datalist.append(float(line))
    datain_txt = np.array(datalist)
    
    src_name = os.path.split(signalSrc)[1]
    # src_path = os.path.split(sinesrc)[0]
    
    print('\nLoaded file: '+src_name)
    
    lgth = len(list(datain_txt))    # get length by iterating csvin obj (only way?)
    print('Length of datain = '+str(lgth))
    
    
#    for j in signalSrcArray:
#        datalist = []
#        with open(signalSrc, mode='r') as infile:
#            for line in infile.readlines():
#                datalist.append(float(line))
#        datain_txt = np.array(datalist)
#
#        src_name = os.path.split(signalSrc)[1]
#        # src_path = os.path.split(sinesrc)[0]
#
#        print('\nLoaded file: '+src_name)
#
#        lgth = len(list(datain_txt))    # get length by iterating csvin obj (only way?)
#        print('Length of datain = '+str(lgth))



    # // *---------------------------------------------------------------------* //
    # *---convert to freq domain & scale---*    
        
    # forward FFT
    datain_txt_fft = sp.fft(datain_txt[0:N])
    datain_txt_fmag = np.abs(datain_txt_fft)
    datain_txt_phase = np.arctan2(datain_txt_fft.imag, datain_txt_fft.real)
    
    # inverse FFT
    datain_txt_inv = sp.ifft(datain_txt_fft)
    
    datain_txt_diff = datain_txt_inv - datain_txt[0:N]
    
    
    # scale and format FFT out for plotting
    #datain_txt_scale = 2.0/N * np.abs(datain_txt_fft[0:int(N/2)])
    datain_txt_scale = 2.0/N * np.abs(datain_txt_fft[0:int(N/2)])
    
    
    print('\nPerformed FFT, calculate Mag and Phase, create scaled signal "datain_txt_scale"')
        
    
    # // *---------------------------------------------------------------------* //
    # *---plots---*
    
    tLen = N
    
    # Input signal
    fig101 = plt.figure(num=101, facecolor='silver', edgecolor='k')
    odmkSrcplt1 = plt.plot(datain_txt[0:tLen])     # scale for viewing peaks
    plt.setp(odmkSrcplt1, color='red', ls='-', linewidth=1.00)
    plt.xlabel('time domain signal (1st '+str(tLen)+' samples)')
    plt.ylabel('Magnitude')
    plt.title('test DATAIN (first '+str(tLen)+' samples)')
    plt.grid(color='c', linestyle=':', linewidth=.5)
    plt.grid(True)
    # plt.xticks(np.linspace(0, Fs/2, 10))
    ax = plt.gca()
    ax.set_axis_bgcolor('black')
    
    
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    
    # FFT Magnitude out plot (0-Fs/2)
    fig102 = plt.figure(num=102, facecolor='silver', edgecolor='k')
    odmkFFTplt1 = plt.plot(xfnyq, datain_txt_scale)
    plt.setp(odmkFFTplt1, color='red', ls='-', linewidth=1.00)
    plt.xlabel('Frequency: 0 - '+str(Fs / 2)+' Hz')
    plt.ylabel('Magnitude (scaled by 2/N)')
    plt.title('test DATAIN - FFT: Fs = '+str(Fs)+', NFFT = '+str(N))
    plt.grid(color='c', linestyle=':', linewidth=.5)
    plt.grid(True)
    # plt.xticks(np.linspace(0, Fs/2, 10))
    ax = plt.gca()
    ax.set_axis_bgcolor('black')
    

#    fnum = 3
#    pltTitle = 'Input Signals: dataArray ('+str(tLen)+' samples)'
#    pltXlabel = 'dataArray time-domain wav'
#    pltYlabel = 'Magnitude'
#    
#    # define a linear space from 0 to 1/2 Fs for x-axis:
#    xaxis = np.linspace(0, tLen, tLen)
#    
#    #pdb.set_trace()
#    
#    odmkplt.odmkMultiPlot1D(fnum, dataArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='cool')


# // *---------------------------------------------------------------------* //




# // *---------------------------------------------------------------------* //

if 0:

    print('\n')
    print('// *--------------------------------------------------------------* //')
    print('// *---:: 2 Mono Signal Overlay ::---*')
    print('// *---:: time domain plot, freq domain plot ::---*')
    print('// *--------------------------------------------------------------* //')
    
    
    # Sampling Freq (must agree with OSC source)
    #Fs = 100000000.0
    Fs = 48000.0
    
    # Number of sample points
    N = 10000
    # sample period
    T = 1.0 / Fs
    
    
    # sinesrc = u'C:\\usr\\eschei\\odmkPython\\odmk\\DSP\\csvsrc\\sintest1.csv'
    #sinesrc = signalSrcDir+'dds_osc_output.txt'
    #sinesrc = signalSrcDir+'dds_sqr_output.txt'
    #sinesrc = signalSrcDir+'pwm_output.txt'
    #sinesrc = signalSrcDir+'dds_lfo_output.txt'
    
    signalsrc1 = signalSrcDir+'in_src_L.dat'
    signalsrc2 = signalSrcDir+'out_dly_L.dat'
    
    # reads .csv data into Numpy array:
    datalist1 = []
    with open(signalsrc1, mode='r') as infile:
        for line in infile.readlines():
            datalist1.append(float(line))
    datain1 = np.array(datalist1)
    
    src_name1 = os.path.split(signalsrc1)[1]
    # src_path = os.path.split(sinesrc)[0]
    
    print('\nLoaded file: '+src_name1)
    
    lgth1 = len(list(datain1))    # get length by iterating csvin obj (only way?)
    print('Length of datain = '+str(lgth1))
    
    
    datalist2 = []
    with open(signalsrc2, mode='r') as infile:
        for line in infile.readlines():
            datalist2.append(float(line))
    datain2 = np.array(datalist2)
    
    src_name2 = os.path.split(signalsrc2)[1]
    # src_path = os.path.split(sinesrc)[0]
    
    print('\nLoaded file: '+src_name2)
    
    lgth2 = len(list(datain2))    # get length by iterating csvin obj (only way?)
    print('Length of datain = '+str(lgth2))    
    
    
    
    dataArray = np.array([])
    dataArray = np.concatenate((datain1, datain2))
    dataArray = dataArray.reshape((2, N))    
    
    
        
    # forward FFT
    datain1_fft = sp.fft(datain1[0:N])
    datain1_fmag = np.abs(datain1_fft)
    datain1_phase = np.arctan2(datain1_fft.imag, datain1_fft.real)
    
    datain2_fft = sp.fft(datain2[0:N])
    datain2_fmag = np.abs(datain2_fft)
    datain2_phase = np.arctan2(datain2_fft.imag, datain2_fft.real)
    
    # inverse FFT
    #datain1_txt_inv = sp.ifft(datain1_txt_fft)
    #datain2_txt_inv = sp.ifft(datain2_txt_fft)
    
    #datain1_txt_diff = datain1_txt_inv - datain1_txt[0:N]
    
    
    # scale and format FFT out for plotting
    #datain_txt_scale = 2.0/N * np.abs(datain_txt_fft[0:int(N/2)])
    datain1_scale = 2.0/N * np.abs(datain1_fft[0:int(N/2)])
    datain2_scale = 2.0/N * np.abs(datain2_fft[0:int(N/2)])
    
    dataFFTMagArray = np.array([])
    dataFFTMagArray = np.concatenate((datain1_scale, datain2_scale))
    dataFFTMagArray = dataFFTMagArray.reshape((2, int(N/2)))
    
    print('\nPerformed FFT, calculate Mag and Phase, create scaled signal "datain_txt_scale"')
        
    
    # // *---------------------------------------------------------------------* //
    # *---plots---*

    
    # // *---------------------------------------------------------------------* //
    # // *---Multi Plot - source signal array vs. FFT MAG out array---*
    # // *---------------------------------------------------------------------* //

    tLen = N

    fnum = 5
    pltTitle = 'Input Signals: dataArray ('+str(tLen)+' samples)'
    pltXlabel = 'dataArray time-domain wav'
    pltYlabel = 'Magnitude'
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xaxis = np.linspace(0, tLen, tLen)
    
    #pdb.set_trace()
    
    odmkplt.odmkMultiPlot1D(fnum, dataArray, xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='cool')
    
    
    fnum = 6
    pltTitle = 'FFT Mag: dataFFTMagArray '
    pltXlabel = 'Frequency: 0 - '+str(Fs / 2)+' Hz'
    pltYlabel = 'Magnitude (scaled by 2/N)'
    
    # define a linear space from 0 to 1/2 Fs for x-axis:
    xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    odmkplt.odmkMultiPlot1D(fnum, dataFFTMagArray, xfnyq, pltTitle, pltXlabel, pltYlabel, colorMp='cool')


# // *---------------------------------------------------------------------* //


# // *---------------------------------------------------------------------* //

plt.show()

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::done::---*')
print('// *--------------------------------------------------------------* //')

# // *---------------------------------------------------------------------* //