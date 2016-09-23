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


# temp python debugger - use >>>pdb.set_trace() to set break
import pdb

# // *---------------------------------------------------------------------* //
# clear_all()

# // *---------------------------------------------------------------------* //

print('// //////////////////////////////////////////////////////////////// //')
print('// *--------------------------------------------------------------* //')
print('// *---::ODMK Signal Generator 1::---*')
print('// *--------------------------------------------------------------* //')
print('// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ //')


# // *---------------------------------------------------------------------* //


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
        usage:
        >>tbSigGen = sigGen.odmkSigGen1(numSamples, fs)
    '''

    def __init__(self, sigLength, fs=48000, rootDir='None'):

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
    # // *---gen simple periodic sin waveforms (sigLength # samples)
    # // *-----------------------------------------------------------------* //

    def monosin(self, freq, fs='None'):
        ''' generates mono sin wav with frequency = "freq"
            optional fs parameter: default = obj global fs
            usage:
            >>tbSigGen = sigGen.odmkSigGen1(numSamples, fs)
            >>testFreq = 5000.0
            >>sin5K = tbSigGen.monosin(testFreq2) '''

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
            optional fs parameter: default = obj global fs
            usage:
            >>tbSigGen = sigGen.odmkSigGen1(numSamples, fs)
            >>testFreq = 5000.0
            >>sin5K = tbSigGen.monocos(testFreq2) '''

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
    # // *---gen simple periodic triangle waveforms (sigLength # samples)
    # // *-----------------------------------------------------------------* //

    def monotri(self, freq, fs='None'):
        ''' generates mono triangle wav with frequency = "freq"
            optional fs parameter: default = obj global fs
            usage:
            >>tbSigGen = sigGen.odmkSigGen1(numSamples, fs)
            >>testFreq = 5000.0
            >>tri5K = tbSigGen.monotri(testFreq) '''
            
            
        if fs != 'None':
            Fs = fs
        else:
            Fs = self.fs
        
        T = 1 / Fs
        
        Tfreq = 1 / freq
        TQtrFreq = (1 / (freq * 4))
        
        cycleSamples = Tfreq / T
        qtrCycleSamples = TQtrFreq / T

        # pdb.set_trace()
        # create a triangle signal
        currentAmp = 0
        monotri = np.array([])
        for i in range(self.sigLength):
            j = (i+1) % cycleSamples
            if (j < qtrCycleSamples) or (j >= (qtrCycleSamples * 3)):
                currentAmp += 1
                monotri = np.append(monotri, currentAmp)
            elif (j >= qtrCycleSamples) or (j < (qtrCycleSamples * 3)):
                currentAmp -= 1
                monotri = np.append(monotri, currentAmp)
        monotri = monotri * (.999 / max(monotri))
        return monotri


    # // *-----------------------------------------------------------------* //
    # // *---TXT write simple periodic waveforms (sigLength # samples)
    # // *-----------------------------------------------------------------* //

    def sig2txt(self, sigIn, nChan, outNm, outDir='None'):
        ''' writes data to TXT file
            signal output name = outNm (expects string) '''

        if outDir != 'None':
            try:
                if isinstance(outDir, str):
                    txtOutDir = outDir
                    os.makedirs(txtOutDir, exist_ok=True)
            except NameError:
                print('Error: outNm must be a string')
        else:
            txtOutDir = self.sigGenOutDir

        try:
            if isinstance(outNm, str):
                sigOutFull = txtOutDir+outNm
        except NameError:
            print('Error: outNm must be a string')

        # writes data to .TXT file:
        outputFile = open(sigOutFull, 'w', newline='')

        if nChan == 0:
            print('ERROR: Number of Channels must be >= 1')
        elif nChan == 1:
            for i in range(len(sigIn)):
                outputFile.write(str(sigIn[i]) + '\n')
        else:
            for i in range(len(sigIn[0])):
                lineTmp = ""
                for j in range(len(sigIn) - 1):
                    strTmp = str(sigIn[j, i]) + str('    ')
                    lineTmp = lineTmp + strTmp
                lineTmp = lineTmp + str(sigIn[len(sigIn) - 1, i]) + '\n'
                outputFile.write(lineTmp)            

        outputFile.close()
        

    # // *-----------------------------------------------------------------* //
    # // *---CSV write simple periodic waveforms (sigLength # samples)
    # // *-----------------------------------------------------------------* //

    def sig2csv(self, sigIn, outNm, outDir='None'):
        ''' writes data to CSV file
            signal output name = outNm (expects string) '''

        if outDir != 'None':
            try:
                if isinstance(outDir, str):
                    csvOutDir = outDir
                    os.makedirs(csvOutDir, exist_ok=True)
            except NameError:
                print('Error: outNm must be a string')
        else:
            csvOutDir = self.sigGenOutDir

        try:
            if isinstance(outNm, str):
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
        
#case 3 then  // store 1 cycle of saw-up
#        table1=[(-0.999:.999/(depth/2-1):.999)'];
#case 4 then  // store 1 cycle of saw-down
#        table1=[(0.999:-.999/(depth/2-1):-.999)'];
#case 5 then  // store 1 cycle of chebychev  
#    for r=1:length(tl)
#        table1(r)=cos(13*acos(tl(r)));
#    end
#case 6 then  // store 1 cycle of pulse1  
#    for s=1:length(t2)
#        table1(s)=sin(5*%pi*(t2(s))); 
#    end 
#case 7 then  // store 1 cycle of pulse2  
#    for s=1:length(t2)
#        table1(s)=sin(9*%pi*(t2(s))); 
#    end 
#case 8 then  // store 1 cycle of pulse3 
#    for s=1:length(t2)
#        table1(s)=sin(23*%pi*(t2(s))); 
#    end
#case 9 then  // store 1 cycle of pulse4  
#    for t=1:length(t3)
#        table1(t)=cos(5*%pi*(t3(t))); 
#    end             
#else
#    for t=1:length(t1)
#        table1(t)=sin(2*%pi*(tl(t))); 
#    end


