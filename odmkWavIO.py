# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

# __::((odmkWavIO.py))::__

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

import os
import wave
import numpy as np

rootDir = 'C:/odmkDev/odmkCode/odmkPython/'
audioScrDir = 'C:/odmkDev/odmkCode/odmkPython/audio/wavsrc/'
audioOutDir = 'C:/odmkDev/odmkCode/odmkPython/audio/wavout/'

# temp python debugger - use >>>pdb.set_trace() to set break
import pdb


# // *---------------------------------------------------------------------* //

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : object definition
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 

class odmkWavIO:
    ''' odmk .wav read and write class
        usage: mySigGen = odmkSigGen1(sigLength, fs, sigType=sin, cmplx=False)
        sigLength => total number of samples
        fs => signal sample rate
        sigType =>
        cmplx => complex signal - real + imag, 90 deg phase shift
        usage:
        >>tbSigGen = sigGen.odmkSigGen1(numSamples, fs)
    '''

    def __init__(self):

        print('\nA odmkWaveIO object has been created')

    # // *-----------------------------------------------------------------* //
    # // *---gen simple periodic sin waveforms (sigLength # samples)
    # // *-----------------------------------------------------------------* //

    def wavRead(self, wavIn, wavInDir, wavLength='None'):
        ''' read a 24 bit wav file into wavOutDir
            wavOut = output file name (txt)
            wavOutDir = output directory (created if doesn't exist)
            fs = sample rate
            wavLength = optional truncated length of output .wav in seconds
            usage:
            >>tbWavIO = sigGen.odmkSigGen1(numSamples, fs) '''

        # Checks: wavOut = 24 bits, wavOut = stereo, wavOut = 44.1K, 48K, 96K

        wavInFull = wavInDir+wavIn

        fwav = wave.open(wavInFull, 'r')
  
        fSampleRate = fwav.getframerate()
        fSampleWidth = fwav.getsampwidth()
        fSampleBits = 8*fSampleWidth
        fChannels = fwav.getnchannels()
        fNumSamples = fwav.getnframes()
        # fparams = fwav.getparams()
        
        wavIn_param = {'fSampleRate': fSampleRate, 'fSampleBits': fSampleBits, 'fChannels': fChannels, 'fNumSamples': fNumSamples}


        wavIn_bytes = fwav.readframes(fNumSamples)

        fwav.close()

        wavIn_bytesL = bytearray([])
        wavIn_bytesR = bytearray([])
        wavIn_smplL = 0
        wavIn_smplR = 0
        wavIn_L = np.array([])
        wavIn_R = np.array([])
        wavIn_stereo = np.array([])

        # de-interleave bytes to L/R channels - convert to float - scale to +/-1.0 
        for ii in range(len(wavIn_bytes)):
            if (ii % 6) < 3:
                wavIn_bytesL.append(wavIn_bytes[ii])
            else:
                wavIn_bytesR.append(wavIn_bytes[ii])
            if (((ii+1) % 6) == 0):
                wavIn_smplL = int.from_bytes(wavIn_bytesL, byteorder='little', signed=True)
                wavIn_smplR = int.from_bytes(wavIn_bytesR, byteorder='little', signed=True)
                wavIn_smplL = float(wavIn_smplL / 2**23)
                wavIn_smplR = float(wavIn_smplR / 2**23)
                wavIn_L = np.append(wavIn_L, wavIn_smplL)
                wavIn_R = np.append(wavIn_R, wavIn_smplR)       
                wavIn_bytesL = bytearray([])
                wavIn_bytesR = bytearray([])    
        wavIn_stereo = np.array([wavIn_L, wavIn_R])

        return wavIn_stereo, wavIn_param

# // *---------------------------------------------------------------------* //        

    def wavWrite(self, wavOut, wavOutNm, wavOutDir, fs, channels='None', wavLength='None'):
        ''' write a 24 bit wav file into wavOutDir
            wavOut = stereo or mono 24-bit sample numpy array
            wavOutNm = output file name
            wavOutDir = output directory (created if doesn't exist)
            fs = sample rate
            wavLength = optional truncated length of output .wav in seconds
            usage:
            >>tbWavIO = sigGen.odmkSigGen1(numSamples, fs) '''

        # Checks: wavOut = 24 bits, wavOut = stereo, wavOut = 44.1K, 48K, 96K

        os.makedirs(wavOutDir, exist_ok=True)
        wavOutFull = wavOutDir+wavOutNm

        fSampleRate = fs
        if channels!='None':
            if (channels < 1 or channels > 2):
                print('ERROR: channels must be either 1 or 2')
                return
            else:
                fChannels = channels
        else:
            fChannels = 2
        if wavLength != 'None':
            if wavLength > wavOut[0]:
                print('ERROR: wavLength must be less than length of wavIn')
                return
            else:
                fNumSamples = wavLength
        else:
            fNumSamples = len(wavOut[0])

        #pdb.set_trace()               

        wavOut_L = wavOut[0]
        wavOut_R = wavOut[1]

        wavOut_bytesL = bytearray(3)
        wavOut_bytesR = bytearray(3)
        wavOut = bytearray([])
        wavOut_smplL = 0
        wavOut_smplR = 0
        
        # interleave bytes to L/R channels - convert to float - scale to +/-1.0 
        for jj in range(len(wavOut_L)):
            wavOut_smplL = int(wavOut_L[jj] * 2**23)
            wavOut_smplR = int(wavOut_R[jj] * 2**23)
            wavOut_bytesL = wavOut_smplL.to_bytes(3, byteorder='little', signed=True)
            for kk in range(3):    
                wavOut.append(wavOut_bytesL[kk])
            wavOut_bytesR = wavOut_smplR.to_bytes(3, byteorder='little', signed=True)   
            for kk in range(3):    
                wavOut.append(wavOut_bytesR[kk])          
            wavOut_bytesL = bytearray([])
            wavOut_bytesR = bytearray([])
        
        fwav = wave.open(wavOutFull, 'w')
        
        fwav.setframerate(fSampleRate)    # set the frame rate
        fwav.setsampwidth(3)            # the sample width (bytes per sample)
        fwav.setnchannels(fChannels)    # set the number of channels
        fwav.setnframes(fNumSamples)    # set the number of frames
        
        # fwav.setparams(fparams)    # set all parameters at once

        fwav.writeframes(wavOut)    # write audio frames and patch up the file header
        
        fwav.close()    # patch up the file header and close the output file

        return

# // *---------------------------------------------------------------------* //

# print('\n')
# print('// *--------------------------------------------------------------* //')
# print('// *---::done::---*')
# print('// *--------------------------------------------------------------* //')

# // *---------------------------------------------------------------------* //