# -*- coding: utf-8 -*-
# *****************************************************************************
# /////////////////////////////////////////////////////////////////////////////
# header begin-----------------------------------------------------------------
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *****************************************************************************
#
# __::((odmkWavGen1.py))::__
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
print('// *---::ODMK Waveform Generator 1::---*')
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

# /////////////////////////////////////////////////////////////////////////////
# #############################################################################
# begin : object definition
# #############################################################################
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 

class odmkWavGen1:
    ''' odmk wav signal generator version 1
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

    def __init__(self, sigLength, fs=48000, outDir='None'):

        # *---set primary parameters from inputs---*

        if outDir != 'None':
            self.odmkWavGen1OutDir = outDir
            os.makedirs(outDir, exist_ok=True)
        else:
            self.rootDir = os.path.dirname(os.path.abspath(__file__))
            self.odmkWavGen1OutDir = os.path.join(self.rootDir, "odmkWavGen1OutDir/")
            os.makedirs(self.odmkWavGen1OutDir, exist_ok=True)
        # myfile_path = os.path.join(self.rootDir, "myfile.txt")

        self.sigLength = sigLength
        self.fs = fs
        
        self.phaseAcc = 0        
        
        print('\nAn odmkSigGen1 object has been instanced with:')
        print('sigLength = '+str(self.sigLength)+'; fs = '+str(fs))

    # #########################################################################
    # end : object definition
    # #########################################################################

    # #########################################################################
    # begin : member method definition
    # #########################################################################


    # //////////////////////////////////////////////////////////////
    # begin : DDS functions
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    
    # //////////////////////////////////////////////////////////////
    # function: tablegen
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    
    def tablegen(self, shape, depth):        
        ''' create look-up table entries for different waveforms

        <<sin, cos, tri, saw-up, saw-dn, exp-up, exp-dn, log-up, log-dn>>
        '''

        table1 = np.zeros(depth)
        
        if shape == 1:    # store 1 cycle of sin
            t1 = np.linspace(0.0, 1.0, depth+1)
            t1 = t1[0:depth]          # crops off the last point for proper cyclic table
            for q in range(depth):
                table1[q] = np.sin(2 * np.pi*(t1[q]))
                
        if shape == 2:    # store 1 cycle of cos
            t1 = np.linspace(0.0, 1.0, depth+1)
            t1 = t1[0:depth]          # crops off the last point for proper cyclic table
            for q in range(depth):
                table1[q] = np.cos(2 * np.pi*(t1[q]))                

        elif shape == 3:    # store 1 cycle of tri (4 cases to handle arbitrary depth)
            if depth % 4 == 0:
                # first quarter cycle + 1
                table1[0:round(depth/4)+1] = np.linspace(0,1,round(depth/4)+1)
                # 2nd & 3rd quarter cycles +1 (overwrite last value of previous w/1)
                table1[round(depth/4):3*(round(depth/4))+1] = np.linspace(1,-1,2*(round(depth/4))+1)
                triQtrTmp = np.linspace(-1,0,round(depth/4)+1)
                table1[3*(round(depth/4)):depth] = triQtrTmp[0:len(triQtrTmp)-1]
            elif depth % 4 == 1:
                table1[0:round(depth/4)+1] = np.linspace(0,1,round(depth/4)+1)
                table1[round(depth/4):3*(round(depth/4))+1] = np.linspace(1,-1,2*(round(depth/4))+1)
                triQtrTmp = np.linspace(-1,0,round(depth/4)+2)
                table1[3*(round(depth/4)):depth] = triQtrTmp[0:len(triQtrTmp)-1]
            elif depth % 4 == 2:
                table1[0:round(depth/4)] = np.linspace(0,1,round(depth/4))
                table1[round(depth/4)-1:3*(round(depth/4))-1] = np.linspace(1,-1,2*(round(depth/4)))
                triQtrTmp = np.linspace(-1,0,round(depth/4)+1)
                table1[3*(round(depth/4))-2:depth] = triQtrTmp[0:len(triQtrTmp)-1]
            elif depth % 4 == 3:
                table1[0:round(depth/4)+1] = np.linspace(0,1,round(depth/4)+1)
                table1[round(depth/4):3*(round(depth/4))+1] = np.linspace(1,-1,2*(round(depth/4))+1)
                triQtrTmp = np.linspace(-1,0,round(depth/4))
                table1[3*(round(depth/4)):depth] = triQtrTmp[0:len(triQtrTmp)-1]

        elif shape == 4:    # store 1 cycle of saw-up
            table1 = np.linspace(-1,1,depth)

        elif shape == 5:    # store 1 cycle of saw-down
            table1 = np.linspace(1,-1,depth)

        elif shape == 6:    # store 1 cycle of chebychev
            t1 = np.linspace(0.0, 1.0, depth+1)
            t1 = t1[0:depth]          # crops off the last point for proper cyclic table
            for r in range(t1):
                table1[r] = np.cos(13 * np.acos(t1[r]))

        elif shape == 7:    # store 1 cycle of pulse1
            t2 = np.linspace(1,0,depth)**3
            for s in range(t2):
                table1[s] = np.sin(5 * np.pi*(t2[s]))

        elif shape == 8:    # store 1 cycle of pulse2
            t2 = np.linspace(1,0,depth)**3
            for s in range(t2):
                table1[s] = np.sin(9 * np.pi*(t2[s]))

        elif shape == 9:    # store 1 cycle of pulse3
            t2 = np.linspace(1,0,depth)**3        
            for s in range(t2):
                table1[s] = np.sin(23 * np.pi*(t2[s]))

        elif shape == 10:    # store 1 cycle of pulse4
            # create a pseudo-symmetrical exponetial pulse
            # crops off the last point for proper cyclic table
            t3_1 = np.linspace(0,1,int(np.floor(depth/2)+1))
            t3_2 = np.linspace(1,0,int(np.ceil(depth/2)+1))
            t3 = np.concatenate(( t3_1[0:len(t3_1)-1], t3_2[0:len(t3_2)-1] ))**3        
            for t in range(t3):
                table1[t] = np.cos(5 * np.pi*(t3[t]))

        else:    # default
            t1 = np.linspace(0.0, 1.0, depth+1)
            t1 = t1[0:depth]          # crops off the last point for proper cyclic table
            for q in range(depth):
                table1[q] = np.sin(2 * np.pi*(t1[q]))
        
        return table1



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


    def multiSin(self, freqArray, fs='None'):
        ''' generates an array of sin waves with frequency = "freq"
            optional fs parameter: default = obj global fs
            usage:
            >>tbSigGen = sigGen.odmkSigGen1(numSamples, fs)
            >>testFreqArray = 5000.0
            >>sin5K = tbSigGen.monosin(testFreq2) '''

        if len(freqArray) <= 1:
            print('ERROR (multiSin): freq must be a list of frequencies')
            return 1

        if fs != 'None':
            Fs = fs
        else:
            Fs = self.fs

        # sample period
        T = 1.0 / Fs

        # create composite sin source
        x = np.linspace(0.0, self.sigLength*T, self.sigLength)
        multiSin = np.zeros([self.sigLength])
        for i in range(len(freqArray)):
            for j in range(self.sigLength):
                curr = multiSin[j]
                multiSin[j] = curr + np.sin(freqArray[i] * 2.0*np.pi * x[j])
        multiSin = (0.999 / max(multiSin)) * multiSin
        return multiSin



    # wavetable oscillator function

    def odmkWTOsc1(self, numSamples, shape, freqCtrl, phaseCtrl, quant='None'):
        
        ''' *--------------------------------------------------------*   
        # odmkWTOsc1: single channel wavetable oscillator
        #
        # shape:
        # <<sin, cos, tri, saw-up, saw-dn, exp-up, exp-dn, log-up, log-dn>>
        #
        # The output frequency can be fixed or variable
        # When freqCtrl is a single scalar value, the output freq. is fixed
        # When freqCtrl is an array of length=numSamples, the output freq. varies each step
        #
        # The output phase can be fixed or variable
        # When phaseCtrl is a single scalar value, the output phase is fixed
        # When phaseCtrl is an array of length=numSamples, the output phase varies each step 
        #
        # If quant != None, quantize output to integer range +/- quant
        #
        # *--------------------------------------------------------* // '''
    
        tableDepth = 4096    

        tb = self.tablegen(shape,tableDepth);

        Fs = self.fs

        # sample period
        #T = 1.0 / Fs

        accWidth = 48
        qntWidth = np.ceil(np.log2(tableDepth))
        lsbWidth = accWidth - qntWidth    # expect 36
        lsbWidthScale = 2**lsbWidth
        lsbWidthUnScale = 2**-lsbWidth

        if not( isinstance(freqCtrl,int) or isinstance(freqCtrl,float) or isinstance(freqCtrl,list) or isinstance(freqCtrl, np.ndarray)):
            print('ERROR (odmkWTOsc1): freqCtrl must be a single freq value, a list, or a numpy array of frequency values')
            return 1            
        elif ( (isinstance(freqCtrl,list) or isinstance(freqCtrl, np.ndarray)) and len(freqCtrl) < numSamples ):
            print('ERROR (odmkWTOsc1): freqCtrl array must be at least numSamples long')
            return 1
        elif (isinstance(freqCtrl,int) or isinstance(freqCtrl,float)):    # fixed freq
            # scale freq to match quantizer
            freqCtrlScale = lsbWidthScale * freqCtrl
            skipInc = ((2**qntWidth) * freqCtrlScale) / Fs
        elif (isinstance(freqCtrl,list) or isinstance(freqCtrl, np.ndarray)):     # variable freq
            freqCtrlScale = lsbWidthScale * freqCtrl[0]           
            skipInc = ((2**qntWidth) * freqCtrlScale) / Fs

        if not( isinstance(phaseCtrl,int) or isinstance(phaseCtrl,float) or isinstance(phaseCtrl,list) or isinstance(phaseCtrl, np.ndarray)):        
            print('ERROR (odmkWTOsc1): phaseCtrl must be a single phase value, or an array of phase values')
            return 1
        elif ( (isinstance(phaseCtrl,list) or isinstance(phaseCtrl, np.ndarray)) and len(phaseCtrl) < numSamples ):
            print('ERROR (odmkWTOsc1): phaseCtrl array must be at least numSamples long')
            return 1
        # initialize variables so that the OSC starts at the beginning of the table            
        elif (isinstance(phaseCtrl,int) or isinstance(phaseCtrl,float)):    # fixed phase            
            # scale phase offset
            # converts radians into a scaled phase offset to be added to the output of the acc
            # dependent on the table depth
            phaseOffset = round( ((phaseCtrl * tableDepth) / (2 * np.pi)) % tableDepth )
            accAddr = phaseOffset * lsbWidthScale    # init to 'zero plus phase shift'
        elif (isinstance(phaseCtrl,list) or isinstance(phaseCtrl, np.ndarray)):
            phaseOffset = round( ((phaseCtrl[0] * tableDepth) / (2 * np.pi)) % tableDepth )
            accAddr = phaseOffset * lsbWidthScale    # init to 'zero plus phase shift'

        # ***initialize***

        # init osc output vectors
        odmkOsc = np.zeros([numSamples])
        odmkOsc90 = np.zeros([numSamples])
        odmkSqrPulse = np.zeros([numSamples])
        if quant!='None': 
            odmkOscQuant = np.zeros([numSamples])
            

        # used to add a 90 deg offset for complex sinusoid generation
        offset90 = tableDepth / 4

        accAddr90 = (offset90 + phaseOffset) * lsbWidthScale    # init to 'zero plus offset90 plus phase shift'
        qntAddr = int(phaseOffset)
        qntAddr90 = int(qntAddr + offset90)
        yLow = tb[qntAddr]
        yLow90 = tb[qntAddr90]

        # **main loop**
        # ___::((Interpolated WaveTable))::___
        # generates main osc waveform out, 90deg shifted out (sin/cos), square pulse out
    
        for i in range(numSamples):    # osc main loop
            
            if (isinstance(freqCtrl,list) or isinstance(freqCtrl, np.ndarray)):
                freqCtrlScale = lsbWidthScale * freqCtrl[i]    # scalar or array

            accAddrP1 = accAddr + lsbWidthScale
            accAddr90P1 = accAddr90 + lsbWidthScale
    
            # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            # generate oscillator output
            # ///////////////////////////////////////////////////        
            if (qntAddr == tableDepth-1):     # avoid addr overflow for linear interpolation addressing - should investigate
                yHigh = tb[qntAddr]
                yHigh90 = tb[qntAddr90 + 1]
                odmkOsc[i] = tb[tableDepth-1]
                # linear interpolation
                odmkOsc90[i] = yLow90 + (yHigh90-yLow90) * ((accAddr90P1-(qntAddr90 * lsbWidthScale)) * lsbWidthUnScale)
                # temp
                #yHigh_tap[i] = yHigh
                #yHigh90_tap[i] = yHigh90
            elif (qntAddr90 == tableDepth-1):
                yHigh = tb[qntAddr + 1]
                yHigh90 = tb[qntAddr90]
                # linear interpolation
                odmkOsc[i] = yLow + (yHigh - yLow) * ((accAddrP1 - (qntAddr * lsbWidthScale)) * lsbWidthUnScale)
                odmkOsc90[i] = tb[tableDepth-1]
            else:
                yHigh = tb[qntAddr + 1]
                yHigh90 = tb[qntAddr90 + 1]
                odmkOsc[i] = yLow + (yHigh-yLow) * ((accAddrP1 - (qntAddr * lsbWidthScale)) * lsbWidthUnScale)
                odmkOsc90[i] = yLow90 + (yHigh90 - yLow90) * ((accAddr90P1 - (qntAddr90 * lsbWidthScale)) * lsbWidthUnScale)

            if quant!='None':
                if isinstance(quant, int):
                    odmkOscQuant[i] = int(round(quant*odmkOsc[i]))
                else:
                    print('quantization value must be an integer')
                    return 1


            # generate square pulse output
            if odmkOsc[i] >= 0:
                odmkSqrPulse[i] = 1
            else:
                odmkSqrPulse[i] = 0
  
            # phase accumulator
            #accAddr = (accAddr + skipInc[i]) % 2**accWidth
            #accAddr90 = (accAddr90 + skipInc[i]) % 2**accWidth
            accAddr = (accAddr + skipInc) % 2**accWidth
            accAddr90 = (accAddr90 + skipInc) % 2**accWidth            
    
            # quantize
            qntAddr = int(np.floor(accAddr * lsbWidthUnScale))
            qntAddr90 = int(np.floor(accAddr90 * lsbWidthUnScale))
    
            yLow = tb[qntAddr]
            yLow90 = tb[qntAddr90]
            # temp
            #yLow_tap[i] = yLow
            #yLow90_tap[i] = yLow90
 
        if quant!='None':
            return odmkOscQuant
        else:
            return odmkOsc, odmkOsc90, odmkSqrPulse


    # #########################################################################
    # begin : waveform generators
    # #########################################################################

    phaseAcc = 0

    def pulseWidthMod(self, phaseInc, pulseWidth):
        ''' generates pulse width modulated square wav
            usage:
            >>phaseInc = output base freq - Fo / Fs
            >>pulseWidth = (% of cycle)/100 -> [0 - 1] '''

        # create pulse width modulated square wave
        if self.phaseAcc >= 1.0:
            self.phaseAcc -= 1
        if self.phaseAcc > pulseWidth:
            pwm = 0
        else:
            pwm = 1
            
        self.phaseAcc += phaseInc

        return pwm



    # #########################################################################
    # begin : file output
    # #########################################################################

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
