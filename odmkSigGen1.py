# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkSigGen1.py))::__
#
# Python Signal Generator
# optional CSV output
# optional .wav output
# optonal Spectral Plotting
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

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::test::---*')
print('// *--------------------------------------------------------------* //')


## Sampling Freq
#Fs = 44100.0
## Number of sample points
#N = 1024
## sample period
#T = 1.0 / Fs
#
#freqs = [5000.0, 7000.0]
#
## create composite sin source
#x = np.linspace(0.0, N*T, N)
#y = np.sin(freqs[0] * 2.0*np.pi*x) + 0.5*np.sin(freqs[1] * 2.0*np.pi*x)


# // *---------------------------------------------------------------------* //

#print('\n')
#print('// *--------------------------------------------------------------* //')
#print('// *---::Load CSV data into "dataout" (numpy array, floats)::---*')
#print('// *--------------------------------------------------------------* //')
#
#rootDir = u'C:/usr/eschei/odmkPython/odmk/audio/'
#outDir = rootDir+'csvgen/'
#os.makedirs(outDir, exist_ok=True)
#
#sigOut = 'sigOutTest1.csv'
#sigOutFull = outDir+sigOut
#
## writes data to .CSV file:
#
#outputFile = open(sigOutFull, 'w', newline='')
#outputWriter = csv.writer(outputFile)
#
#for i in range(N):
#    tmpRow = [y[i]]
#    outputWriter.writerow(tmpRow)
#
#
#outputFile.close()

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : object definition
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


class odmkSigGen1:
    ''' odmk signal generator version 1
        outputs a fs periodic signal normalized to +/- 1
        usage: mySigGen = odmkSigGen1(sigLength, fs, sigType=sin, cmplx=False)
        sigLength => total number of samples
        fs => signal sample rate
        sigType =>
        <<sin, cos, tri, sqr, saw-up, saw-dn, exp-up, exp-dn, log-up, log-dn>>
        cmplx => complex signal - real + imag, 90 deg phase shift
    '''

    def __init__(self, sigLength, fs=44100, rootDir='None'):

        # *---set primary parameters from inputs---*

        if rootDir != 'None':
            self.rootDir = rootDir
            os.makedirs(rootDir, exist_ok=True)
        else:
            self.rootDir = os.path.dirname(os.path.abspath(__file__))
        self.sigGenOutDir = os.path.join(self.rootDir, "sigGenOutDir/")
        os.makedirs(self.sigGenOutDir, exist_ok=True)
        # myfile_path = os.path.join(self.rootDir, "myfile.txt")

        self.sigLength = sigLength
        self.fs = fs
        print('\nAn odmkSigGen1 object has been instanced with:')
        print('sigLength = '+str(self.sigLength)+'; fs = '+str(fs))

    # // *-----------------------------------------------------------------* //
    # // *---gen simple periodic waveforms (sigLength # samples)
    # // *-----------------------------------------------------------------* //

    def monosin(self, freq, fs='None'):
        ''' generates mono sin wav with frequency = "freq"
            optional fs parameter: default = obj global fs '''

        if fs != 'None':
            Fs = fs
        else:
            Fs = self.fs

        # sample period
        T = 1.0 / Fs

        # create composite sin source
        x = np.linspace(0.0, self.sigLength*T, self.sigLength)
        monosin = np.sin(freq * 2.0*np.pi * x)

        return monosin

    def monocos(self, freq, fs='None'):
        ''' generates mono cos wav with frequency = "freq"
            optional fs parameter: default = obj global fs '''

        if fs != 'None':
            Fs = fs
        else:
            Fs = self.fs

        # sample period
        T = 1.0 / Fs

        # create composite sin source
        x = np.linspace(0.0, self.sigLength*T, self.sigLength)
        monocos = np.cos(freq * 2.0*np.pi * x)

        return monocos

    # // *-----------------------------------------------------------------* //
    # // *---CSV write simple periodic waveforms (sigLength # samples)
    # // *-----------------------------------------------------------------* //

    def sig2csv(self, sigIn, outNm, outDir='None'):
        ''' generates mono sin wav with frequency = "freq"
            signal output name = outNm (expects string)
            optional fs parameter: default = obj global fs '''

        if outDir != 'None':
            try:
                if not(isinstance(outDir, str)):
                    csvOutDir = outDir
                    os.makedirs(csvOutDir, exist_ok=True)
            except NameError:
                print('Error: outNm must be a string')
        else:
            csvOutDir = self.sigGenOutDir

        try:
            if not(isinstance(outNm, str)):
                sigOutFull = csvOutDir+outNm
        except NameError:
            print('Error: outNm must be a string')

        # writes data to .CSV file:
        outputFile = open(sigOutFull, 'w', newline='')
        outputWriter = csv.writer(outputFile)

        for i in range(len(sigIn)):
            tmpRow = [sigIn[i]]
            outputWriter.writerow(tmpRow)

        outputFile.close()
