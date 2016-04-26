# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

# __::((native_wav_test.py))::__

# Python .wav file read / write

# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header end-------------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************

import wave
import scipy as sp
import matplotlib.pyplot as plt

from odmkClear import *

# temp python debugger - use >>>pdb.set_trace() to set break
#import pdb

# // *---------------------------------------------------------------------* //

clear_all()

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

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::Clone .wav file & write"::---*')
print('// *--------------------------------------------------------------* //')


#wavfile_out = 'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\wavsrc\\werk\\wavclone000.wav'
wavfile_out = 'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\wavsrc\\werk\\wavclone000.wav'

fclone = wave.open(wavfile_out, 'w')


fclone.setframerate(fSampleRate)    # set the frame rate
fclone.setsampwidth(fSampleWidth)    # the sample width
fclone.setnchannels(fChannels)    # set the number of channels
fclone.setnframes(fNumSamples)    # set the number of frames

# fclone.setparams(fparams)    # set all parameters at once


fclone.writeframes(audio_data)    # write audio frames and patch up the file header

fclone.close()    # patch up the file header and close the output file

# // *---------------------------------------------------------------------* //

wavClone1_in = 'C:\\usr\\eschei\\odmkPython\\odmk\\audio\\wavsrc\\werk\\wavclone000.wav'

fclone1 = wave.open(wavClone1_in, 'r')

fc1SampleRate = fclone1.getframerate()
fc1SampleWidth = fclone1.getsampwidth()
fc1SampleBits = 8*fSampleWidth
fc1Channels = fclone1.getnchannels()
fc1NumSamples = fclone1.getnframes()
fc1params = fclone1.getparams()

print('\nClone Sample Rate = '+str(fc1SampleRate))
print('Clone Sample Bits = '+str(fc1SampleBits))
print('Clone num Channels = '+str(fc1Channels))
print('Clone num Frames = '+str(fc1NumSamples))
print('Clone Parameters = '+str(fc1params))

clone1_data = fclone1.readframes(fc1NumSamples)

fclone1.close()

# // *---------------------------------------------------------------------* //

print('\n')
print('// *--------------------------------------------------------------* //')
print('// *---::done::---*')
print('// *--------------------------------------------------------------* //')

# // *---------------------------------------------------------------------* //

























