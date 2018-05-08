# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkHilbertSSB.py))::__
#
# Frequency shift a signal using SSB modulation.
#
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header end-------------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

import os
import wave
import numpy as np
import scipy as sp
import scipy.signal
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


rootDir = 'C:/odmkDev/odmkCode/odmkPython/'

import sys

#sys.path.insert(0, 'C:/odmkDev/odmkCode/odmkPython/util')
sys.path.insert(0, rootDir+'util')
from odmkClear import *
#from odmkPlotUtil import *
import odmkPlotUtil as odmkplt

#sys.path.insert(1, 'C:/odmkDev/odmkCode/odmkPython/DSP')
sys.path.insert(1, rootDir+'DSP')
import odmkClocks as clks
import odmkSigGen1 as sigGen


# temp python debugger - use >>>pdb.set_trace() to set break
import pdb

# // *---------------------------------------------------------------------* //
plt.close('all')
clear_all()

# // *---------------------------------------------------------------------* //


# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : function definitions
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def nextpow2(x):
    """Return the first integer N such that 2**N >= abs(x)"""
    return int(np.ceil(np.log2(np.abs(x))))


def freq_shift(x, f_shift, dt):
    """
    Shift the specified signal by the specified frequency.
    """
    # Pad the signal with zeros to prevent the FFT invoked by the transform from
    # slowing down the computation:
    N_orig = len(x)
    N_padded = 2**nextpow2(N_orig)
    t = np.arange(0, N_padded)
    hFilter = scipy.signal.hilbert(np.hstack((x, np.zeros(N_padded-N_orig, x.dtype))))
    #pdb.set_trace()
    return (hFilter*np.exp(2j*np.pi*f_shift*dt*t))[:N_orig].real    
    
    #return (scipy.signal.hilbert(np.hstack((x, np.zeros(N_padded-N_orig, x.dtype))))*np.exp(2j*np.pi*f_shift*dt*t))[:N_orig].real


# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# end : function definitions
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


#dt = 1e-3
#fs = 1/dt

fs = 48000
dt = 1/fs

T = 1.0
t = np.arange(0, T, dt)
N = len(t)

# Construct original signal:
x = 3*np.cos(2*np.pi*t)+np.cos(2*np.pi*3*t)+2*np.cos(2*np.pi*5*t)

# Uncomment the code below to construct a more interesting signal:
# N_taps = 2500
# np.random.seed(1)
# b = sp.signal.remez(N_taps, [0.0, 14.5/fs, 15.5/fs, 19.5/fs, 20.5/fs, 0.5], [0, 1, 0])
# x = sp.signal.lfilter(b, 1.0, np.random.rand(len(t)))
    
# Frequency shift:
f_shift = 10.0

# Shift signal's frequency components by using the Hilbert transform
# to perform SSB modulation:
x_shift = freq_shift(x, f_shift, dt)



fHilbert=sp.signal.remez(99, [0.03, 0.47],[1], type='hilbert');
#freqz(b,axisFreqz=[0,np.pi,-50,10],axisPhase=[0,np.pi,-np.pi,np.pi])




# Plot results:
f = np.fft.fftfreq(N, dt)
xf = np.fft.fft(x).real
xf_shift = np.fft.fft(x_shift).real
start = 0
stop = int((25.0/(fs/2.0))*(N/2.0))
plt.clf()
plt.subplot2grid((2, 3), (0, 0), rowspan=2, colspan=2)
plt.plot(t, x, 'b', t, x_shift, 'r-')
plt.xlabel('t (s)')
plt.ylabel('x(t)')
plt.legend(('original', 'shifted'))
plt.subplot2grid((2, 3), (0, 2))
plt.stem(f[start:stop], xf[start:stop])
plt.title('Original')
plt.subplot2grid((2, 3), (1, 2))
plt.stem(f[start:stop], xf_shift[start:stop])
plt.title('Shifted')
plt.xlabel('F (Hz)')
plt.tight_layout()
plt.suptitle('Frequency Shifting Using SSB Modulation')
plt.draw()



# // *---------------------------------------------------------------------* //

# // *---------------------------------------------------------------------* //
# // *---Mono FFT plots---*
# // *---------------------------------------------------------------------* //

# define a sub-range for wave plot visibility
tLen = len(t)

fnum = 2
pltTitle = 'Input Signal x (first '+str(tLen)+' samples)'
pltXlabel = 'x'
pltYlabel = 'Amplitude'


sig = x[0:tLen]
# define a linear space from 0 to 1/2 Fs for x-axis:
xaxis = np.linspace(0, tLen, tLen)


odmkplt.odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)



# define a sub-range for wave plot visibility
t2Len = len(fHilbert)

fnum = 3
pltTitle = 'Hilbert Filter Coefficients ('+str(t2Len)+' taps)'
pltXlabel = 'fHilbert'
pltYlabel = 'Amplitude'


sig = fHilbert
# define a linear space from 0 to 1/2 Fs for x-axis:
xaxis = np.linspace(0, t2Len, t2Len)


odmkplt.odmkPlot1D(fnum, sig, xaxis, pltTitle, pltXlabel, pltYlabel)


# // *---------------------------------------------------------------------* //
# // *---Multi Plot - source signal array vs. FFT MAG out array---*
# // *---------------------------------------------------------------------* //

xArray = np.array([x, x_shift])

fnum = 4
pltTitle = 'SSB Hilbert Transform frequency shift'
pltXlabel = 'x vs. x_shift ('+str(tLen)+' samples)'
pltYlabel = 'Amplitude'

# define a linear space from 0 to 1/2 Fs for x-axis:
xaxis = np.linspace(0, tLen, tLen)

#odmkplt.odmkMultiPlot1D(fnum, xArray[0:tLen,:], xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='gnuplot')
odmkplt.odmkMultiPlot1D(fnum, xArray[0:tLen,:], xaxis, pltTitle, pltXlabel, pltYlabel, colorMp='cool')


plt.show()