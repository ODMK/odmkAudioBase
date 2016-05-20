# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkFFTiFFTx.py))::__
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
print('// *---::ODMK Fourier Transforms::---*')
print('// *--------------------------------------------------------------* //')
print('// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ //')


# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Basic Python/Scipy FFT test::---*')
print('// *--------------------------------------------------------------* //')


# Sampling Freq
Fs = 44100.0
# Number of sample points
N = 1024
# sample period
T = 1.0 / Fs

freqs = [5000.0, 7000.0]

# create composite sin source
x = np.linspace(0.0, N*T, N)
y = np.sin(freqs[0] * 2.0*np.pi*x) + 0.5*np.sin(freqs[1] * 2.0*np.pi*x)

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
plt.xlabel('mixed sine: '+str(freqs[0])+' + '+str(freqs[1])+' Hz')
plt.ylabel('Magnitude')
plt.title('Input Signal (first '+str(tLen)+' samples)')
plt.grid(color='c', linestyle=':', linewidth=.5)
plt.grid(True)
# plt.xticks(np.linspace(0, Fs/2, 10))
ax = plt.gca()
ax.set_axis_bgcolor('black')



# define a linear space from 0 to 1/2 Fs for x-axis:
xfnyq = np.linspace(0.0, 1.0/(2.0*T), N/2)


# FFT Magnitude out plot (0-Fs/2)
fig2 = plt.figure(num=2, facecolor='silver', edgecolor='k')
odmkFFTplt1 = plt.plot(xfnyq, yfscale)
plt.setp(odmkFFTplt1, color='red', ls='-', linewidth=1.00)
plt.xlabel('Frequency: 0 - '+str(Fs / 2)+' Hz')
plt.ylabel('Magnitude (scaled by 2/N)')
plt.title('Scipy FFT: Fs = '+str(Fs)+' N = '+str(N))
plt.grid(color='c', linestyle=':', linewidth=.5)
plt.grid(True)
# plt.xticks(np.linspace(0, Fs/2, 10))
ax = plt.gca()
ax.set_axis_bgcolor('black')


# Output signal vs. Input signal
fig3 = plt.figure(num=3, facecolor='silver', edgecolor='k')
odmkSrcplt2 = plt.plot(y[0:tLen])
plt.setp(odmkSrcplt2, color='red', ls='-', linewidth=1.00)
odmkResplt1 = plt.plot(yInv[0:tLen])
plt.setp(odmkResplt1, color='orange', ls='-', linewidth=1.00)
plt.xlabel('mixed sine: '+str(freqs[0])+' + '+str(freqs[1])+' Hz')
plt.ylabel('Magnitude')
plt.title('Input (red) VS Output (orange) (first '+str(tLen)+' samples)')
plt.grid(color='c', linestyle=':', linewidth=.5)
plt.grid(True)
# plt.xticks(np.linspace(0, Fs/2, 10))
ax = plt.gca()
ax.set_axis_bgcolor('black')


# Output vs. Input difference
fig4 = plt.figure(num=4, facecolor='silver', edgecolor='k')
odmkDiffplt = plt.plot(yDiff)
plt.setp(odmkDiffplt, color='orange', ls='-', linewidth=1.00)
plt.xlabel('time domain')
plt.ylabel('Magnitude')
plt.title('Output VS Input Difference')
plt.grid(color='c', linestyle=':', linewidth=.5)
plt.grid(True)
# plt.xticks(np.linspace(0, Fs/2, 10))
ax = plt.gca()
ax.set_axis_bgcolor('black')

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Load CSV data into "datain" (numpy array, floats)::---*')
print('// *--------------------------------------------------------------* //')

sinesrc = u'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\csvsrc\\sintest1.csv'

# reads .csv data into Numpy array:
datalist = []
with open(sinesrc, mode='r') as infile:
    for line in infile.readlines():
        l, dumb = line.strip().split(',')
        datalist.append(float(l))
datain = np.array(datalist)

src_name = os.path.split(sinesrc)[1]
# src_path = os.path.split(sinesrc)[0]

print('\nLoaded file: '+src_name)

lgth = len(list(datain))    # get length by iterating csvin obj (only way?)
print('Length of datain = '+str(lgth))


## t1 = np.arange(0.0, 5.0, 0.02)
#fig1 = plt.figure(num=1, facecolor='silver', edgecolor='k')
## odmkLinePlt1 = plt.plot(t1, sinesrc)
#odmkLinePlt1 = plt.plot(datain[0:999])
## plt.setp(odmkLinePlt1, color='white', ls='-', marker='o', mfc='orange', linewidth=1.00)
#plt.setp(odmkLinePlt1, color='orange', ls='-', linewidth=1.00)
#plt.xlabel('x-axis-label')
#plt.ylabel('y-axis-label')
#plt.title('sinesrc first 1000 samples')
## plt.annotate('sinesrc',
##             xy=(2.23, .74),
##             xycoords='data', color='midnightblue', ha='center',
##             fontsize=20)
## plt.text(1.3, .54, r'$\aleph\ =  f(t): red dots$', color='black', fontsize=16)
## plt.text(1.3, .44, r'$\beth\ = -f(t): blue triangles$', color='black', fontsize=16)
## plt.axis([0, 5, -1, 1])
#plt.grid(True)
#ax = plt.gca()
#ax.set_axis_bgcolor('black')
    

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Read .wav file & extract parameters"::---*')
print('// *--------------------------------------------------------------* //')


# wavfile_in = 'C:\\odmk\\odmkPython\\audio\\wavSource\\odiusEMattack002.wav'
# wavfile_in = 'C:\\odmk\\odmkPython\\audio\\wavSource\\death_by_black_hole_vox.wav'
# wavfile_in = 'C:\\odmk\\odmkPython\\audio\\wavSource\\YT-KONGA2x02.wav'
# wavfile_in = 'C:\\odmk\\odmkPython\\audio\\wavSource\\601_90ALTOSAX_C_chop.wav'
# wavfile_in = 'C:\\odmk\\odmkPython\\audio\\wavSource\\The_Amen_Break_odmk.wav'
wavfile_in = 'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\wavsrc\\MAXIMALL_Ddrumz3.wav'
# wavfile_in = 'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\wavsrc\\rmeWav\\100_24\\100_24.wav'



fwav = wave.open(wavfile_in, 'r')

fSampleRate = fwav.getframerate()
fSampleWidth = fwav.getsampwidth()
fSampleBits = 8*fSampleWidth
fChannels = fwav.getnchannels()
fNumSamples = fwav.getnframes()
fparams = fwav.getparams()

print('\nWav Sample Rate = '+str(fSampleRate))
print('Wav Sample Bits = '+str(fSampleBits))
print('Wav num Channels = '+str(fChannels))
print('Wav num Frames = '+str(fNumSamples))
print('Wav Parameters = '+str(fparams))

audio_data = fwav.readframes(fNumSamples)

fwav.close()

# // *---------------------------------------------------------------------* //

#print('\n')
#print('// *--------------------------------------------------------------* //')
#print('// *---::Clone .wav file & write, reopen to verify parameters"::---*')
#print('// *--------------------------------------------------------------* //')
#
#
##wavfile_out = 'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\wavsrc\\werk\\wavclone000.wav'
#wavfile_out = 'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\wavsrc\\werk\\wavclone000.wav'
#
#fclone = wave.open(wavfile_out, 'w')
#
#
#fclone.setframerate(fSampleRate)    # set the frame rate
#fclone.setsampwidth(fSampleWidth)    # the sample width
#fclone.setnchannels(fChannels)    # set the number of channels
#fclone.setnframes(fNumSamples)    # set the number of frames
#
## fclone.setparams(fparams)    # set all parameters at once
#
#
#fclone.writeframes(audio_data)    # write audio frames and patch up the file header
#
#fclone.close()    # patch up the file header and close the output file

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::sinesrc processing - FFT <-> iFFT, mag & phase plots"::---*')
print('// *--------------------------------------------------------------* //')

halflength = round(len(datain) / 2)

print('\nPerform FFT & calculate Mag and Phase')
print('Create: "datainFF", "datainMAG", "datainPHASE" (numpy array, float)')
datainFFT = sp.fft(datain)
datainMAG = sp.absolute(datain)
datainPHASE = np.arctan2(datainFFT.imag, datainFFT.real)


dataout = sp.ifft(datainFFT)

print('\nCheck results: max(dataout - datain) = '+str(max(dataout - datain)))

#xf = np.linspace(0.0, 1.0/(2.0*T), N/2)



## t1 = np.arange(0.0, 5.0, 0.02)
#fig2 = plt.figure(num=2, facecolor='silver', edgecolor='k')
## odmkLinePlt1 = plt.plot(t1, sinesrc)
#odmkLinePlt1 = plt.plot(datain[0:999])
## plt.setp(odmkLinePlt1, color='white', ls='-', marker='o', mfc='orange', linewidth=1.00)
#plt.setp(odmkLinePlt1, color='orange', ls='-', linewidth=1.00)
#plt.xlabel('x-axis-label')
#plt.ylabel('y-axis-label')
#plt.title('sinesrc first 1000 samples')
## plt.annotate('sinesrc',
##             xy=(2.23, .74),
##             xycoords='data', color='midnightblue', ha='center',
##             fontsize=20)
## plt.text(1.3, .54, r'$\aleph\ =  f(t): red dots$', color='black', fontsize=16)
## plt.text(1.3, .44, r'$\beth\ = -f(t): blue triangles$', color='black', fontsize=16)
## plt.axis([0, 5, -1, 1])
#plt.grid(True)
#ax = plt.gca()
#ax.set_axis_bgcolor('black')


# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::wav processing - FFT <-> iFFT, mag & phase plots"::---*')
print('// *--------------------------------------------------------------* //')

halflength = round(len(datain) / 2)

print('\nPerform FFT & calculate Mag and Phase')
print('Create: "datainFF", "datainMAG", "datainPHASE" (numpy array, float)')
datainFFT = sp.fft(datain)
datainMAG = sp.absolute(datain)
datainPHASE = np.arctan2(datainFFT.imag, datainFFT.real)


dataout = sp.ifft(datainFFT)

print('\nCheck results: max(dataout - datain) = '+str(max(dataout - datain)))

#xf = np.linspace(0.0, 1.0/(2.0*T), N/2)



## t1 = np.arange(0.0, 5.0, 0.02)
#fig2 = plt.figure(num=2, facecolor='silver', edgecolor='k')
## odmkLinePlt1 = plt.plot(t1, sinesrc)
#odmkLinePlt1 = plt.plot(datain[0:999])
## plt.setp(odmkLinePlt1, color='white', ls='-', marker='o', mfc='orange', linewidth=1.00)
#plt.setp(odmkLinePlt1, color='orange', ls='-', linewidth=1.00)
#plt.xlabel('x-axis-label')
#plt.ylabel('y-axis-label')
#plt.title('sinesrc first 1000 samples')
## plt.annotate('sinesrc',
##             xy=(2.23, .74),
##             xycoords='data', color='midnightblue', ha='center',
##             fontsize=20)
## plt.text(1.3, .54, r'$\aleph\ =  f(t): red dots$', color='black', fontsize=16)
## plt.text(1.3, .44, r'$\beth\ = -f(t): blue triangles$', color='black', fontsize=16)
## plt.axis([0, 5, -1, 1])
#plt.grid(True)
#ax = plt.gca()
#ax.set_axis_bgcolor('black')


# // *---------------------------------------------------------------------* //

plt.show()

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::done::---*')
print('// *--------------------------------------------------------------* //')

# // *---------------------------------------------------------------------* //