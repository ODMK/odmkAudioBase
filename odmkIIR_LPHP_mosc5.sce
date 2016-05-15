///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : header
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////
//
//
//___::((odmkIIR_LPHP_mosc4))::___

//__<<name=> "odmkIIR_LPHP_mosc4" (.sce)>>__

//___::((JIROBATA Programming Industries))::___
//___::((ODMK:odorousbeast:BarutanBreaks:djoto:2014:2015:2016))::___
//___::((created by eschei))___



//___::((IIR 2nd Order, Low-pass/High-pass filter Scilab script))::___


//include components: 
//___::((odmk_osc3))::___
//___::((odmk_multiOsc4))::___
//___::((odmk_IIR_2ndorder))::___


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//This script implements 2nd Order IIR Filtering examples
//the filters are parameterizable
//The script generates osc waveforms & tests filtering
//performs FFT for analysis & plotting
//Plots transfer functions, mag & phase plots


//Real time controls: - cutoff, type: LP/HP, 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : header
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


xdel(winsid()) //-> closes all open graphs (= matlab.close('all'))
clear all;
clc;


exec('/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/osc/odmk_osc5.sce');


///////////////////////////////////////////////////////////////////////////////////////////////
//begin : stack adjustment
///////////////////////////////////////////////////////////////////////////////////////////////

//try for ~10MBytes
//stacksize input is doubles
//
max_bytes = 25*10e6;
//max_bits=16;
max_bits=24;
bytespersample=ceil(max_bits/8);
max_data_bytes=max_bytes-(12+24);    //total size - header,format,etc. approx
max_stack=max_data_bytes/8;

//if size > max then
//    error('wav file too large');
//else
stacksize(max_stack)
//    stacksize('max')
//    stacksize('min')
//    sz=stacksize()
//end

///////////////////////////////////////////////////////////////////////////////////////////////
//end : stack adjustment
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : function definitions
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : IIR LPHP filter functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//A 2nd order IIR filter described by the transfer function:

//temp
/////////////////////////////////////
//tunable 2nd order allpass filter
/////////////////////////////////////
//A(z) = [-c + d*(1-c)*(z^-1) + (z^-2)]/[1 + d*(1-c)*(z^-1) - c*(z^-2)]

//parameter c - adjusts the bandwidth
//c = [tan(%pi*fb/fs) - 1]/[tan(%pi*fb/fs) + 1]

//parameter d - adjusts cutoff frequency
//d = -cos*(2*%pi*fc/fs)

//which leads to the following difference equations
//y(n) = -c*x(n) + d*(1-c)*x(n-1) + x*(n-2) -d*(1-c)*y(n-1) + c*y(n-2)
//temp 


/////////////////////////////////////
//Generic 2nd order filter
/////////////////////////////////////
//A(z) = [b0 + b1*(z^-1) + b2*(z^-2)]/[1 + a1*(z^-1) + a2*(z^-2)]

//which leads to the following difference equations
//y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) - a1*y(n-1) - a2*y(n-2)


function [b,a] = odmk_IIR_2ndorder(fc,fs,ftype)
//    //Derive coefficients for a shelving filter with a given amplitude 
//    //and cutoff frequency. All coefficients are calculated as 
//    //described in Zolzer’s DAFX book (p. 50 -55). 
//    //
//    //Usage: [b,a] = odmk_IIR_2ndorder(fc,fs,ftype);
//    //
//    //        b: filter zeros (numerator coefficients)
//    //        a: filter poles (denominator coefficients)
//
//    //        FC is the center frequency
//    //        Fs is the sampling rate
//    //        Choices are: ’LP’ or ’HP’
//
//
//    //Error Check 
    if((strcmp(ftype,'LP')) & (strcmp(ftype,'HP')))
        error(['Unsupported Filter Type: ' ftype]);
    end

    K = tan((%pi*fc)/fs);
    K2 = K^2;
    r2 = sqrt(2);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    Lowpass
    //%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (~(strcmp(ftype,'LP')))        
        b0 = K2/(1 + r2*K + K2);
        b1 = (2*K2)/(1 + r2*K + K2);
        b2 = K2/(1 + r2*K + K2);
        a0 = 1;
        a1 = (2*(K2 - 1))/(1 + r2*K + K2);
        a2 = (1 - r2*K + K2)/(1 + r2*K + K2);

//        b0 = 0.1013120;
//        b1 = 0.1244110;
//        b2 = 0.1013120;
//        a0 = 1;
//        a1 = -1.125814;
//        a2 = 0.4812868;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    Highpass
    //%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (~(strcmp(ftype,'HP')))
        b0 = 1/(1 + r2*K + K2);
        b1 = -2/(1 + r2*K + K2);
        b2 = 1/(1 + r2*K + K2);
        a0 = 1;
        a1 = (2*(K2 - 1))/(1 + r2*K + K2);
        a2 = (1 - r2*K + K2)/(1 + r2*K + K2);
    end
          
    //%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    ALL-PASS
    //%%%%%%%%%%%%%%%%%%%%%%%%%% 
//    else
//        b0 = V0;
//        b1 = 0;
//        b2 = 0;
//        a1 = 0;
//        a2 = 0;
//    end      
    
    //return values    
    a=[a0,a1,a2];    //den
    b=[b0,b1,b2];    //num

endfunction

    
////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : IIR LPHP filter functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : function definitions
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : main script
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : odmk fir 3 main script
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//sampling frequency
fs = 44100;
Nfft = 16384; //define length of FFT for analysis

////////////////////////////////////////////////////////////////
//begin: generate test signal using odmk_osc1
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
//main: osc1
//oscillator outputs sin and cos waveforms - uses 1/4 wave symmetry
//and linear interpolation


//DDS implementation outputs sin, cos. 
//fo = output frequency
//fs = sampling frequency (default=>44100), defined in parent
//td = table depth

//fo = fs/(2^N)*(delta(acc))
//delta(acc) = floor(fo*((2^N)/fs)+0.5)


//output frequency: fo = fs/(2^N)*(phacc))
//create an array of frequency objects - fo_n = [fo_1,fo_2,...,fo_n] - assume Hertz
//fObj = [77,254,440,528,666];
//fObj = [1666,3666,5666,7666,9666];
foArray = [666,2333,3000,7770,11000,16666];
phArray = [0, 0, 0, 0, 0, 0];
//phArray = [0, %pi/3, 2*%pi/3, %pi, 4*%pi/3, 5*%pi/3];

nChan = length(foArray);

//number of steps through the algorithm -> determines length of output sample
nItera = 16384;    //2048 4096 8192 16384 32768 65536 131072 262144

//Nfft = 16384; //define length of FFT for analysis - defined in parent
 
td = 4096;    //table depth

//define the waveform shape stored in the table by selecting from list below:
//1=sin
//2=tri
//3=saw-up
//4=saw-down
//5=chebychev
//6=pulse1
//7=pulse2
//8=pulse3
//9=pulse4
shape = 1;    //waveform shape  
 


//initialize oscillator wavetable output
odmkOsc = zeros(nItera,nChan);
//initialize oscillator square pulse output
odmkSqrPulse = zeros(nItera,nChan);

//initialize Sum of nChan waveforms
oscSum = zeros(nItera,1);
oscSumObj = zeros(nChan,1);
//initialize Mix of 2 waveforms
oscMix1 = zeros(nItera,1);

//init analysis of osc output
if nItera > Nfft then
    odmkOscFFT = zeros(Nfft,nChan);
    odmkSqrPulseFFT = zeros(Nfft,nChan);
    odmkOscMag = zeros(Nfft,nChan);
    odmkSqrPulseMag = zeros(Nfft,nChan);
else     
    odmkOsc1FFT = zeros(nItera,nChan);
    odmkSqrPulseFFT = zeros(nItera,nChan);
    odmkOscMag = zeros(nItera,nChan);
    odmkSqrPulseMag = zeros(nItera,nChan);
end

//init filter outputs
oscFiltx = zeros(nItera,1);
zfx = zeros(nItera,1);

//init wave file outputs 
odmkOscWav = zeros(nItera,2*nChan);
odmkSqrPulse = zeros(nItera,2*nChan);

odmkOsc1_out = zeros(nItera,2);
odmkOsc2_out = zeros(nItera,2);

oscSumWav = zeros(nItera,2);
oscMix1Wav = zeros(nItera,2);




// main osc function call
[odmkOsc, odmkOsc90, odmkSqrPulse] = odmk_multiOsc(nItera, fs, td, shape, foArray, phArray)




//Sum and mix ouput waveforms
for m=1:nItera
    //accumulate row of osc value
    oscSumObj=0;
    for n=1:nChan
        oscSumObj = oscSumObj+odmkOsc(m,n);
        //oscSumSx(m) = oscSumSx(m)+oscSumObj(n);
    end
    oscSumSx(m) = oscSumObj;
    oscMix1(m) = odmkOsc(m,1)*odmkOsc(m,2);
end
oscSum = (1/max(abs(oscSumSx))).*oscSumSx

  

//analysis
//check spectrum
//if nItera > Nfft then
//    for k=1:nChan
//        odmkOscFFT(:,k) = fft(odmkOsc(1:Nfft,k));
//        odmkOscMag(:,k) = abs(odmkOscFFT(1:Nfft,k));
//        odmkOscPhase(:,k) = atan(imag(odmkOscFFT(1:Nfft,k)),real(odmkOscFFT(1:Nfft,k)));
//    end    
//else
//    for k=1:nChan    
//        odmkOscFFT(:,k) = fft(odmkOsc(:,k));
//        odmkOscMag(:,k) = abs(odmkOscFFT(:,k));
//        odmkOscPhase(:,k) = atan(imag(odmkOscFFT(:,k)),real(odmkOscFFT(:,k)));
//    end
//end


////////////////////////////////////////////////////////////////
//end: generate test signal using odmk_osc1
////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////
//begin: filter .wav audio data using ::((odmk_IIR_2ndorder))::
////////////////////////////////////////////////////////////////

//implements shelving LP or HP filtering


    //Derive coefficients for a shelving filter with a given amplitude 
    //and cutoff frequency. All coefficients are calculated as 
    //described in Zolzer’s DAFX book (p. 50 -55). 
    //
    //Usage: [B,A] = shelving(G, Fc, Fs, Q, type);
    //
    //        

    //        G is the logrithmic gain (in dB
    //        FC is the center frequency
    //        Fs is the sampling rate
    //        Q adjusts the slope be replacing the sqrt(2) term
    //        type is a character string defining filter type
    //        Choices are: ’Base_Shelf’ or ’Treble_Shelf’


//initializw coefficients & filter type
//fs = 44100;    //defined above
fc = 0.15*fs;    //cutoff frequency 
ftype = "LP";    //"LP","HP"

//Filter func call - static fc
[iirZero,iirPole] = odmk_IIR_2ndorder(fc,fs,ftype)


//Filter func call - variable fc
xfxzro =((0:0.5/(nItera):0.5)).*fs;
xfx = xfxzro(2:length(xfxzro));

iirZerox = zeros(nItera,3);
iirPolex = zeros(nItera,3);

for j=1:length(xfx)
    [iirZerox(j,:),iirPolex(j,:)] = odmk_IIR_2ndorder(xfx(j),fs,ftype)
end


////generate filter freq response - iir example use
//hz=iir(2,'lp','cheb1',[.1 0],[.050 .05]);
//[hzm,fr]=frmag(hz,256);
//plot(fr,hzm)
////hz=iir(3,'bp','ellip',[.15 .25],[.08 .03]);
////[hzm,fr]=frmag(hz,256);
////plot(fr,hzm,'r')
//
//hzm_dB = 20*log10(hzm);


//[y,zf] = filter(num,den,x [,zi])
[oscFilt,zf] = filter(iirZero,iirPole,oscSum)

//[y,zf] = filter(num,den,x [,zi])
//[oscFiltx,zfx] = filter(iirZerox,iirPolex,oscSum)


//hz=iir(2,'lp','ellip',[.1 .2],[.08 .03]);
//[hzm,fr]=frmag(hz,512);
//scf(1)
//plot(fr,hzm,'r')


//temp - check 2nd order A(z) frequency response
[hzmIir,frIir]=frmag(iirZero,iirPole,512);

scf(2)
plot(frIir,hzmIir,'red')

//test filter with impulse response


////check coefficients using impulse response
//figure(16)
//plot(lpFIR1coef,'cyan')
//plot(odmk_fir1_out(1:nTap),'red')
//a = gca();
//xlabel('Example 8, odmk FIR coefficients vs FIR impulse response output - M=63');
//xgrid(2)

if nItera >= Nfft then
        oscSumFFT = fft(oscSum(1:Nfft));
        oscSumMag = abs(oscSumFFT);
        oscSumMagdB = 20*log10(oscSumMag);
        oscSumPhase = atan(imag(oscSumFFT(1:Nfft)),real(oscSumFFT(1:Nfft)));
        
        oscMix1FFT = fft(oscMix1(1:Nfft));
        oscMix1Mag = abs(oscMix1FFT);
        oscMix1MagdB = 20*log10(oscMix1Mag);
        oscMix1Phase = atan(imag(oscMix1FFT(1:Nfft)),real(oscMix1FFT(1:Nfft)));
        
        //processed audio using odmk filter
        oscFiltFFT = fft(oscFilt(1:Nfft));
        oscFiltMag = abs(oscFiltFFT);
        oscFiltMagdB = 20*log10(oscFiltMag);
        oscFiltPhase = atan(imag(oscFiltFFT(1:Nfft)),real(oscFiltFFT(1:Nfft)));
        
        freqAxis=fs*(0:(Nfft-1))/Nfft; //associated frequency vector   
else   
        oscSumFFT = fft(oscSum);
        oscSumMag = abs(oscSumFFT);
        oscSumMagdB = 20*log10(oscSumMag);
        oscSumPhase = atan(imag(oscSumFFT),real(oscSumFFT));
        
        oscMix1FFT = fft(oscMix1);
        oscMix1Mag = abs(oscMix1FFT);
        oscMix1MagdB = 20*log10(oscMix1Mag);
        oscMix1Phase = atan(imag(oscMix1FFT),real(oscMix1FFT));
        
        //processed audio using odmk filter
        oscFiltFFT = fft(oscFilt);
        oscFiltMag = abs(oscFiltFFT);
        oscFiltMagdB = 20*log10(oscFiltMag);
        oscFiltPhase = atan(imag(oscFiltFFT),real(oscFiltFFT));
        
        freqAxis=fs*(0:(length(oscMix1FFT)-1))/length(oscMix1FFT); //associated frequency vector
end


//'error' between reference filter and odmk filter
//lp8_diff=odmk_fir1_out-lp8_ref;
//wavAfltL_diff=wavAfltL-wavAfltL_ref;
//wavAfltR_diff=wavAfltR-wavAfltR_ref;


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//end: filter .wav audio data using (odmk_shelving_2ndorder)
////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : main script
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : plotting
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : osc test signal filtering plots
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////plot coefficients (from scilab 'wfir' function)
//figure(1)
//plot(lpFIR1coef)
//a = gca();
//xlabel('Example 8, odmk FIR coefficients, M=63');
//xgrid(2)
//
//
////plot filter frequency response (from 'frmag' scilab function)
//figure(2)
//plot(fr99,lpFIR1mag_dB)
//a = gca();
//xlabel('Normalized Digital Frequency fr ');
//ylabel('Magnitude');
//title('Ex8, Freq Response of FIR LPF using wfir, M = 63, cutoff @ 0.25')
//xgrid(2)
//

//plot fft mag of input .wav audio vs. filtered audio
scf(3)
    plot(freqAxis,oscSumMag,'blue')
    //a=get("current_axes");//get the handle of the newly created axes
    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
xlabel("$oscSumMag - FFT\ Output\ Magnitude : 0\le x\le fs\ (Sample Rate)$","fontsize",4,"color","red");
ylabel("FFT Out","fontsize",3,"color","cyan");
title("$oscSum\ FFT\ Mag$","color","black","fontsize",3); 
a=gca(); // Handle on axes entity
//a.grid=[2,2];

//plot fft mag of input .wav audio vs. filtered audio
scf(4)
    plot(freqAxis,oscFiltMag,'cyan')
    //a=get("current_axes");//get the handle of the newly created axes
    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
xlabel("$oscSumMag - FFT\ Output\ Magnitude : 0\le x\le fs\ (Sample Rate)$","fontsize",4,"color","red");
ylabel("FFT Out","fontsize",3,"color","cyan");
title("$oscFilt\ FFT\ Mag\ (blue)$","color","black","fontsize",3); 
a=gca(); // Handle on axes entity
//a.grid=[2,2];

//plot fft mag of input .wav audio vs. filtered audio
scf(5)
    plot(freqAxis,oscSumMag,'blue')
    plot(freqAxis,oscFiltMag,'cyan')
    //a=get("current_axes");//get the handle of the newly created axes
    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
xlabel("$oscSumMag - FFT\ Output\ Magnitude : 0\le x\le fs\ (Sample Rate)$","fontsize",4,"color","red");
ylabel("FFT Out","fontsize",3,"color","cyan");
title("$oscSum\ FFT\ Mag\ (blue)\ vs.\ oscFilt\ FFT\ Mag\ (cyan)$","color","black","fontsize",3); 
a=gca(); // Handle on axes entity
//a.grid=[2,2];

//
////plot reference filter output vs. odmk filter output
//scf(4)
//    plot(freqAxis,wavAfltL_refMag,'blue')
//    plot(freqAxis,wavAfltR_refMag,'cyan')
//    plot(freqAxis,wavAfltLMag,'red')
//    plot(freqAxis,wavAfltRMag,'magenta')
//    //a=get("current_axes");//get the handle of the newly created axes
//    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
//xlabel("$FFT Output Magnitude : 0\le x\le fs (Sample Rate)$","fontsize",4,"color","red");
//ylabel("FFT Out","fontsize",3,"color","cyan");
//title("$chAData FFT Mag (blue/cyan) vs. wavflt FFT Mag (red/magenta)$","color","black","fontsize",3); 
//a=gca(); // Handle on axes entity
//a.grid=[2,2];
//
//
////plot the difference 'error' between reference filter and odmk filter
//figure(5)
//plot(wavAfltL_diff,'red')
//plot(wavAfltL_diff,'magenta')
//a = gca();
//xlabel('odmk_fir1 filtered output vs. scilab convol results difference');
//
//
//////plot the difference 'error' between reference filter and odmk filter
////figure(6)
////plot(wavflt_diff)
////a = gca();
////xlabel('Example 8, odmk_fir1 filtered output vs. scilab convol results difference');


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : osc test signal filtering plots
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//scf(18)
//    plot(freqAxis,oscSumMag)
//    a=gca(); // Handle on axes entity
//    a.x_location = "origin";
//    //a=get("current_axes");//get the handle of the newly created axes
//    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
//xlabel("$oscSum - FFT Output Magnitude : 0\le x\le fs (Sample Rate)$","fontsize",4,"color","red");
//ylabel("FFT Out","fontsize",3,"color","cyan");
//title("oscSum FFT Out Magnitude","color","black","fontsize",3); 
////a=gca(); // Handle on axes entity
////a.x_location = "origin";
//a.grid=[2,2];

//scf(5)
//    plot(freqAxis,oscSumPhase)
//    a=gca(); // Handle on axes entity
//    a.x_location = "origin";
//    //a=get("current_axes");//get the handle of the newly created axes
//    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
//xlabel("$oscSum - FFT Output Phase : 0\le x\le fs (Sample Rate)$","fontsize",4,"color","red");
//ylabel("FFT Out","fontsize",3,"color","cyan");
//title("FFT Out","color","red","fontsize",3); 
////a=gca(); // Handle on axes entity
////a.x_location = "origin";
//a.grid=[2,2];
//
//scf(6)
//    plot(freqAxis,oscMix1Mag)
//    a=gca(); // Handle on axes entity
//    a.x_location = "origin";
//    //a=get("current_axes");//get the handle of the newly created axes
//    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
//xlabel("$oscMix1 - FFT Output Magnitude : 0\le x\le fs (Sample Rate)$","fontsize",4,"color","red");
//ylabel("FFT Out","fontsize",3,"color","cyan");
//title("FFT Out","color","red","fontsize",3); 
////a=gca(); // Handle on axes entity
////a.x_location = "origin";
//a.grid=[2,2];
//
//scf(7)
//    plot(freqAxis,oscMix1Phase)
//    a=gca(); // Handle on axes entity
//    a.x_location = "origin";
//    //a=get("current_axes");//get the handle of the newly created axes
//    //a.data_bounds=[0,-ceil(max(real(odmkOscFFT(:,k))));8192,ceil(max(real(odmkOscFFT(:,k))))];
//xlabel("$oscMix1 - FFT Output Phase : 0\le x\le fs (Sample Rate)$","fontsize",4,"color","red");
//ylabel("FFT Out","fontsize",3,"color","cyan");
//title("FFT Out","color","red","fontsize",3); 
////a=gca(); // Handle on axes entity
////a.x_location = "origin";
//a.grid=[2,2];



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : plotting
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : auxiliary
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : auxiliary
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////