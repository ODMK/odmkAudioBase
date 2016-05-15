// /////////////////////////////////////////////////////////////////////////////////////////////
// #############################################################################################
// begin : header
// #############################################################################################
// /////////////////////////////////////////////////////////////////////////////////////////////
//
// __<<name=> "odmk_STFT_iSTFT_phrot" (.sce)>>__

// ___::((JIROBATA Programming Industries))::___
// ___::((ODMK:odorousbeast:BarutanBreaks:djoto:2014:2015:2016))::___
// ___::((created by eschei))___

// impllementation of a STFT - short-time fourier transform <-> iSTFT re-synthesis
// rotate bins between each frame


// include components:

// ___::((STFT))::___
// ___::((iSTFT))::___
// ___::((odmkWaterfall))::___
//

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This script implements short-time fourier transform (STFT) and
// inverse short-time fourier transform (iSTFT) using overlap-add re-synthesis
// an input waveform is broken into overlapped time-slices, which are then windowed and
// fourier transformed using FFT. The results are plotted using odmkWaterfall for 3D
// visualizations. Re-synthesis is performed by using overlap-add method and iFFT.

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// /////////////////////////////////////////////////////////////////////////////////////////////
// #############################################################################################
// end : header
// #############################################################################################
// /////////////////////////////////////////////////////////////////////////////////////////////


xdel(winsid()) //-> closes all open graphs (= matlab.close('all'))
clear;
clc;


exec('/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/osc/odmk_osc5.sce');


// /////////////////////////////////////////////////////////////////////////////////////////////
// begin : stack adjustment
// /////////////////////////////////////////////////////////////////////////////////////////////

// //try for ~10MBytes
// //stacksize input is doubles
// //
// max_bytes = 20*10e7;
// max_bits=24;
// //bits=24;
// bytespersample=ceil(max_bits/8);
// max_data_bytes=max_bytes-(12+24);    //total size - header,format,etc. approx
// max_stack=max_data_bytes/8;
//
// //if size > max then
// //    error('wav file too large');
// //else
// stacksize(max_stack)
// //    stacksize('max')
// //    stacksize('min')
// //    sz=stacksize()
// //end

// /////////////////////////////////////////////////////////////////////////////////////////////
// end : stack adjustment
// /////////////////////////////////////////////////////////////////////////////////////////////


// /////////////////////////////////////////////////////////////////////////////////////////////
// #############################################################################################
// begin : function definitions
// #############################################################################################
// /////////////////////////////////////////////////////////////////////////////////////////////


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// begin : cyclicZn functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cZn] = cyclicZn(n)
    // generate cyclic group Zn member elements
    
    cZn = zeros(n, 1);
    for k = 1:n
        //z(k) = %e^(((k-1)*2*%pi*%i)/n)                          //Define cyclic group Zn points
        cZn(k) = cos(((k - 1) * 2*%pi)/n) + %i*sin(((k - 1) * 2*%pi)/n)   //Euler's identity    
    end
    
endfunction

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// end : cyclicZn functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rotFrame] = rotateFrame(inFrame, rotVal)

    rotFrame = zeros(length(inFrame), 1)
    for j = 1:length(inFrame)
        if j <= rotVal then
            //rotFrame(length(inFrame) - rotVal-1 + j) = inFrame(j);
            rotFrame(length(inFrame) - rotVal + j) = inFrame(j);
        else
            rotFrame(j - rotVal) = inFrame(j);
        end
    end

endfunction

//j=1
//rotframe(97) = inframe(1)
//
//j=32
//rotframe(128) = inframe(32)
//
//j=33
//rotframe(1) = inframe(33)
//
//j=128
//rotframe(96) = inframe(128)


// //////////////////////////////////////////////////////////////
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// begin : phase vocoder functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// //////////////////////////////////////////////////////////////
// function: phaseVoc
// phase vocoder top-level function - Uses pvMod, STFT, iSTFT

function wavDataOut = phaseVoc(wavDataIn, Ra, Rs, nfft, aWin, sWin, zp, fplots)
    //    wavDataOut = phaseVoc(wavDataIn, Ra, Rs, nfft, aWin, sWin, zp, fplots)
    //    wavDataIn is an input waveform. nfft is the FFT size.  
    //    phaseVoc calls STFT and iSTFT and optional plotting
    //    Ra => analysis hop factor, spaces time instances t(a,u) = u*Ra, where u is 1,2,3,4,5,.... 
    //    Rs => synthesis hop factor, spaces time instances t(s,u) = u*Rs, where u is 1,2,3,4,5,....
    //    nfft => length of fft
    //    aWin => analysis window [0=square, integer /=0 = cosine window, length of aWin, or win array]     
    //    sWin => synthesis window

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //anaysis stage
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // overlapped time-slices of wavDataIn are windowed, (zero-padded), then FFT is performed 
    // we need 25% window overlap for smooth reconstruction(??) -> hop = nfft/4;

    mprintf('\nBegin phase Vocoder (phaseVoc)\n');

    wavInLength = length(wavDataIn);    //total length of input audio sample - analysis 
    mprintf('\nphaseVoc.wavInLength = length of wavIn = %i\n',wavInLength);

    // Calculate the basic STFT, magnitude scaled
    wavX = STFT(wavDataIn', nfft, aWin, Ra, zp, fs);    //check if iSTFT needs to change with zp?
    //wavX = STFT(wavDataIn', nfft, 0, Ra, fs);
    
    [nSTFTRow, nSTFTCol] = size(wavX);
    //nRow = number of frequency bins = nfft/2 + 1
    //nCol = number of audio time-slices (windowed and fft)
    
    mprintf('\nphaseVoc.STFT.nSTFTRow = number of rows of STFT output matrix (nBins) = %i\n',nSTFTRow);
    mprintf('phaseVoc.STFT.nSTFTCol = number of columns of STFT output matrix (nSlice)= %i\n',nSTFTCol);
    

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //synthesis
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // Invert to a waveform
    //wavDataOut = iSTFT(pvocX, nfft, nfft, Rs)';
    wavDataOut = iSTFT(wavX, nfft, sWin, Rs, zp)';

    //crop wavDataOut to match length of wavDataIn
    wavDataOut = wavDataOut(1:length(wavDataIn));

    wavOutLength = length(wavDataOut);    //total length of output audio sample - Synthesis
    mprintf('phaseVoc.wavOutLength = number of columns of iSTFT output array = %i\n',wavOutLength);    

    //STFT plotting
    if fplots==1 then
        if zp ==0 then
            nfftx = nfft;
        else
            nfftx = nfft*2;
        end
        odmkWaterfall(wavX, nfftx)
    end    

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // check results
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    scf(40)    
        plot(wavDataIn,'blue')
        plot(wavDataOut,'red')
        title("STFT - iSTFT compare - wavDataIn - blue; wavDataOut - red","color","black","fontsize",3);                

    wavError = wavDataIn-wavDataOut;
    scf(41)  
        plot(wavError)
        title("STFT - iSTFT error - wavDataIn-wavDataOut","color","black","fontsize",3);


endfunction


function stftMatrixOut = STFT(wavIn, nfft, aWin, Ra, zp, fs)
    //    stftMatrixOut = STFT(wavDataIn, nfft, aWin, Ra, zp, fs)  Short-time Fourier transform.
    //    Returns frames of the short-time Fourier transform of wavIn
    //    wavIn currently a real-only input waveform (ex. mono audio data)
    //    rows of the result are the indexed frequency bin data outputs
    //    column of the result is one nfft-point fft (default 256); each
    //    successive frames are offset by Ra samples until entire wavIn is processed  
    //    frames are hann-windowed at nfft points when aWin is integer/=0, or rectangular if W=0
    //    or with user defined window if aWin is a vector.
    //    currently uses built-in scilab fft

    mprintf('\n// *------------------------------------------------------------------------* //\n')
    mprintf('\nBegin Short Time Fourier Transform (STFT)\n')


    //if nargin < 2;  nfft = 256; end
    //if nargin < 3;  aWin = nfft; end
    //if nargin < 4;  Ra = 0; end
    //if nargin < 5;  fs = 8000; end

    // if wavIn is a column, convert to row
    if size(wavIn,1) > 1
        wavIn = wavIn';
    end


    //calculate window
    //either a specific window array is passed to the STFT function,
    //or a window of length aWin is generated
    //currently expecting aWin = nfft, or zero (square case)
    if length(aWin) == 1    
        if aWin == 0
            // special case: rectangular window
            aWindow = ones(1,nfft);
        else
            halfLength = floor(aWin/2);   // midpoint of win
            //generate 1/2 of a cos window function
            halfwin = 0.5*(1 + cos(%pi*(0:halfLength)/halfLength));
            aWindow = zeros(1,aWin);
            //fill in top-half of windoe array
            aWindow((halfLength+1):(halfLength+halfLength)) = halfwin(1:halfLength);
            //fill in bottom-half of windoe array
            aWindow((halfLength+1):-1:2) = halfwin(1:halfLength);
        end
    else
        aWindow = aWin;    //in case of explicit window array passed to funcfrun
    end

    wavInLength = length(wavIn);    //total length of audio - analysis
    mprintf('\nLength of wavin = %i\n', wavInLength)
    
    // ::Calculate wavIn Zero Padding::
    // zero pad wavIn to take care of last frame
    // adds enough zeros so that wav length is evenly divided by overlapped windows
    // no trailing samples left unprocessed
    // in case when wavIn is already an even muliple of overlapped windows,
    // nothing should be added, wavInZpLength=wevInLength

    Ra_hop = fix(Ra * nfft);    // num samples each hop
    mprintf('Ra_hop = %i\n', Ra_hop)
    extn1 = 3 * Ra_hop;
    wavInLessExtn1 = wavInLength-extn1;    // final overlap
    

    wavInZpad = modulo(Ra_hop - modulo(wavInLessExtn1, Ra_hop), Ra_hop);
    mprintf('temp wavInZpad = %i\n', wavInZpad)
    
    wavInZp = {wavIn,zeros(1, wavInZpad)};
    wavInZpLength = wavInLength + wavInZpad;
    mprintf('phaseVoc.STFT.wavInZpLength = length of internal wav (zero-padded wav)= %i\n', wavInZpLength);
    
    odmkSTFTWinPlot(aWindow, Ra, wavInZpLength)


//// TEMP C++
//    for(posin=posout=0; posin < input_size; posin+=hopsize) {
//        mod = posin%fftsize;
//        //c++ window and rotate frame
//        for(i=0; i < fftsize; i++)
//            if(posin+i < input_size)
//                sigframe[(i+mod)%fftsize] = input[posin+i]*window[i];
//            else
//                sigframe[(i+mod)%fftsize] = 0;
//// TEMP end

    scf(500)
    plot(wavInZp(1:nfft));
    title("$wavInZP$","color","black","fontsize",3);

    //perform windowing and fft on slices of audio data
    nCol = 1;

   //check if zp is enabled
    if zp == 0 then
        //no frame zero padding
        // pre-allocate output array - exploit symmetry - only save lower bins + 1 (center bin)
        stftMatrix = zeros((1 + nfft/2),1 + fix((wavInZpLength-nfft)/(Ra * nfft)));
        //b runs length of sample minus nfft
        for b = 0:(Ra * nfft):(wavInZpLength - nfft)
            //create temp matrix u = window dotMPY with rotated audio time-slice             
            mod = modulo(b, nfft)
            u_norot = aWindow.*wavInZp((b + 1):(b + nfft));    // before rot mod
            //u(modulo(b+mod,nfft)) = aWindow.*wavInZP((b+1):(b+nfft));
            u = rotateFrame(u_norot, mod);
            // ***TEMP visual inspection***
//            if b == 0*Ra*nfft then
//                mprintf('\nWindowed sig with no rotation (length = %i)', length(u_norot))
//                //disp(u_norot')
//                scf(555)
//                plot(u_norot);
//                title("$u_norot check$","color","black","fontsize",3);
//                mprintf('\nWindowed sign with rotation (length = %i)', length(u))
//                disp(u)
//                scf(556)
//                plot(u);
//                title("$u check$","color","black","fontsize",3);                
//            end
            // ***TEMPEND***
            X = fft(u_norot');
            //truncate top half of freq bins
            //stftMatrix(:,nCol) = X(1:(1+nfft/2))';
            stftMatrix(:, nCol) = X(1:(1 + nfft/2));
            nCol = nCol + 1;
        end;
    else
        //frame zero padding
        // pre-allocate output array - exploit symmetry - only save lower bins + 1 (center bin)
        stftMatrix = zeros((1+nfft),1+fix((wavInZpLength-nfft)/(Ra*nfft)));
        //b runs length of sample minus nfft 
        for b = 0:(Ra*nfft):(wavInZpLength-nfft)        
            //create temp matrix u = window dotMPY with current audio time-slice with zero padding
            mod = modulo(b,nfft)
            u_norot = aWindow.*wavInZp((b+1):(b+nfft));
            //u = aWindow.*wavInZP((b+mod+1):(b+mod+nfft));    // ***check this***
            u = rotateFrame(u_norot, mod);
            X = fft([u,zeros(1,nfft)]);
            //truncate top half of freq bins (length of X actually 2*nfft)
            stftMatrix(:,nCol) = X(1:(1+nfft))';
            nCol = nCol+1;
        end;
    end;

    //scf(501)
    //plot(stftMatrix(:,1));
    //title("$stftMatrix(:,1)$","color","black","fontsize",3);

    //odmkSTFTFramePlot(stftMatrix, 1, "stftMatrix frame 1", 501)
    odmkSTFTFramePlot(stftMatrix, 2, 501)


//    else
//        //frame zero padding
//        // pre-allocate output array - exploit symmetry - only save lower bins + 1 (center bin)
//        stftMatrix = zeros((1+nfft),1+fix((wavInZpLength-nfft)/(Ra*nfft)));
//        //b runs length of sample minus nfft 
//        for b = 0:(Ra*nfft):(wavInZpLength-nfft)        
//            //create temp matrix u = window dotMPY with current audio time-slice with zero padding           
//            u = aWindow.*wavInZP((b+1):(b+nfft));
//            X = fft([u,zeros(1,nfft)]);
//            //truncate top half of freq bins (length of X actually 2*nfft)
//            stftMatrix(:,nCol) = X(1:(1+nfft))';
//            nCol = nCol+1;
//        end;
//    end;


    stftMatrixOut = stftMatrix;
    
    [nSTFTRow, nSTFTCol] = size(stftMatrixOut);
    //nRow = number of frequency bins = nfft/2 + 1
    //nCol = number of audio time-slices (windowed and fft)

    mprintf('\nCreated 2D STFT Output Array <<stftMatrix>>, dimentions: %i rows, %i cols', nSTFTRow, nSTFTCol)

    mprintf('\nEnd Short Time Fourier Transform (STFT)\n')
    mprintf('\n// *------------------------------------------------------------------------* //\n')


    //number of frequency bins = nfft/2 = rows of stft output matrix
    //[stftBin, stftSlice] = size(stftMatrix);
    
    //mprintf('STFT.stftBin = number of rows of STFT stftMatrix matrix = %i\n',stftBin);
    //mprintf('STFT.stftSlice = number of cols of STFT stftMatrix matrix = %i\n',stftSlice);

    //Matlab: plot a spectrogram
    //    tt = [0:size(stftMatrix,2)]*Ra/fs;
    //    ff = [0:size(stftMatrix,1)]*fs/nfft;
    //    imagesc(tt,ff,20*log10(abs(stftMatrix)));
    //    axis('xy');
    //    xlabel('time / sec');
    //    ylabel('freq / Hz')


endfunction


function wavOut = iSTFT(istftMatrixIn, nfft, sWin, Rs, zp)
    //    wavOut = iSTFT(istftMatrixIn, nfft, nWin, Rs)  Inverse short-time Fourier transform.
    //    Performs overlap-add resynthesis from the short-time Fourier transform 
    //    rows of istftMatrixIn are the indexed frequency bin data from STFT
    //    columns of istftMatrixIn are one nfft-point fft frame from STFT    
    //    frames are hann-windowed at nfft points when nWin is integer/=0, or rectangular if W=0
    //    or with user defined window if nWin is a vector.
    //    ??sort out over-lap add windows??
    //    currently uses built-in scilab ifft

    mprintf('\n// *------------------------------------------------------------------------* //\n');
    mprintf('Begin Inverse Short Time Fourier Transform (iSTFT)\n');

    //if nargin < 2; nfft = 2*(size(d,1)-1); end
    //if nargin < 3; w = 0; end
    //if nargin < 4; Rs = 0; end  // will become winlen/2 later

    [niRow, niCol] = size(istftMatrixIn);
    //iSTFTRow = number of frequency bins = nfft/2 + 1
    //iSTFTCol = number of audio time-slices (windowed and fft)
    
    
//    if niRow ~= ((nfft/2)+1 | (nfft+1))
//        error('number of rows should be nfft/2+1 for non-zp, nfft+1 for zp')
//    end

    //calculate window
    //either a specific window array is passed to the STFT function,
    //or a window of length nWin is generated
    //currently expecting nWin = nfft, or zero (square case)
    if length(sWin) == 1    
        if sWin == 0
            // special case: rectangular window
            synWin = ones(1,nfft);
        else
            halfLength = floor(sWin/2);   // midpoint of win
            //generate 1/2 of a cos window function
            halfwin = 0.5*(1 + cos(%pi*(0:halfLength)/halfLength));
            synWin = zeros(1,sWin);
            //fill in top-half of windoe array
            synWin((halfLength+1):(halfLength+halfLength)) = halfwin(1:halfLength);
            //fill in bottom-half of windoe array
            synWin((halfLength+1):-1:2) = halfwin(1:halfLength);
        end
    else
        synWin = sWin;    //in case of explicit window array passed to funcfrun
    end

    sLength = nfft + (niCol-1)*(Rs*nfft);
    mprintf('\nphaseVoc.STFT.sLength = length of internal wavOut (zero-padded)= %i\n',sLength);
        
    wavOut = zeros(1,sLength);
    Rs_hop = fix(Rs * nfft);    // num samples each hop
    //b runs length of sample minus nfft
    for b = 0:niCol-1
        v = istftMatrixIn(:,b+1)';
        //create full frame from conjugate of lower bins (fft symmetry)
        if zp == 0 then
            v = [v, conj(v([((nfft/2)):-1:2]))];
            // = real(ifft(v));
            F = real(ifft(v));
        else
            v = [v, conj(v([((nfft)):-1:2]))];
            Fzp = real(ifft(v));
            F_rot = Fzp(1:nfft);
        end
        // unrotate frames
        mod = modulo(b*Rs_hop, nfft)
        //F = rotateFrame(F_rot, mod);
        // ***TEMP visual inspection***
//        if b == 0 then
//            mprintf('\niSTFT sig with rotation (length = %i)', length(F_rot))
//            mprintf('direct input to ifft')
//            disp(v')
//            scf(557)
//            plot(F_rot);
//            title("$u norot check$","color","black","fontsize",3);
//            mprintf('\niSTFT sign un-rotation (length = %i)', length(F))
//            //disp(F)
//            scf(558)
//            plot(F);
//            title("$u check$","color","black","fontsize",3);                
//        end
        // ***TEMPEND***
        //F = F'
        wavOut((b*(Rs*nfft)+1):(b*(Rs*nfft)+nfft)) = wavOut((b*(Rs*nfft)+1):(b*(Rs*nfft)+nfft))+F.*synWin;
    end;
       
    
    //sLength = length(wavOut);    //total length of audio - synthesis
    odmkiSTFTWinPlot(synWin,Rs,sLength)
    
endfunction


// //////////////////////////////////////////////////////////////             
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// end : phase vocoder functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



// //////////////////////////////////////////////////////////////             
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// begin : STFT frame plot function
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function odmkSTFTFramePlot(stftMatrixIn, frame, idx)
    // plots one frame of a STFT matrix
    // stftMatrixIn is a 3D stack of FFT frames (STFT output)
    // frame is the frame # to plot
    // idx is the index used by the scilab plotting function to make unique plot window

    scf(idx)
    plot(stftMatrix(:,frame));
    title("$stftMatrix single frame:$","color","black","fontsize",3);

endfunction

// //////////////////////////////////////////////////////////////             
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// end : STFT window plot functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



// //////////////////////////////////////////////////////////////             
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// begin : STFT window plot function
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function odmkSTFTWinPlot(win,Ra,wavInZpLength)
    
    //mprintf('phaseVoc.wavInZpLength = length of wavIn = %i\n',wavInZpLength);

    //plot window information
    tw = (1:1:length(win));
    scf(2)
    plot(tw,win,'red')
    title("$STFT Analysis Window$","color","black","fontsize",3);
    xlabel("$Time axis : 0\le x\le window length (Nfft)$","fontsize",3,"color","black");
    ylabel("$Window Amplitude$","fontsize",3,"color","black");

 
    //plot a window for each audio time-slice
    tt = (1:1:wavInZpLength);
    scf(3)    
    winPlotSum = zeros(1,wavInZpLength);
    for i = 0:(Ra*nfft):(wavInZpLength-nfft)      
        winPlotArray = zeros(1,wavInZpLength);
        winPlotArray((i+1):(i+nfft)) = win;
        winPlotSum = winPlotSum+winPlotArray;
        plot2d(tt,[winPlotArray],style=[5]); 
    end;
    plot2d(tt,winPlotSum,style=[1]);
    title("$STFT overlapped Analysis Windows (red) ; overlapped sum (black)$","color","black","fontsize",3);
    xlabel("$Time axis : 0\le x\le wav length$","fontsize",3,"color","black");
    ylabel("$Window Amplitude$","fontsize",3,"color","black"); 

endfunction

// //////////////////////////////////////////////////////////////             
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// end : STFT window plot functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// //////////////////////////////////////////////////////////////             
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// begin : iSTFT window plot function
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function odmkiSTFTWinPlot(win,Rs,sLength)

    //plot window information
    tw = (1:1:length(win));
    scf(30)
    plot(tw,win,'red')
        title("$iSTFT Synthesis Window$","color","black","fontsize",3);
    xlabel("$Time axis : 0\le x\le window length (Nfft)$","fontsize",3,"color","black");
    ylabel("$Window Amplitude$","fontsize",3,"color","black");

 
    //plot a window for each audio time-slice
    tt = (1:1:sLength);
    scf(31)    
    winPlotSum = zeros(1,sLength);
    for i = 0:(Ra*nfft):(sLength-nfft)      
        winPlotArray = zeros(1,sLength);
        winPlotArray((i+1):(i+nfft)) = win;
        winPlotSum = winPlotSum+winPlotArray;
        plot2d(tt,[winPlotArray],style=[5]); 
    end;
    plot2d(tt,winPlotSum,style=[1]);
    title("$iSTFT overlapped Synthesis Windows (red) ; overlapped sum (black)$","color","black","fontsize",3);
    xlabel("$Time axis : 0\le x\le wav length$","fontsize",3,"color","black");
    ylabel("$Window Amplitude$","fontsize",3,"color","black"); 

    scf(32)
    plot2d(tt,winPlotSum,style=[1]);
//    title("$iSTFT overlapped windows sum - time domain$","color","blue","fontsize",3);
//    xlabel("$Time axis : 0\le x\le wav length$","fontsize",3,"color","black");
//    ylabel("$Window Amplitude$","fontsize",3,"color","black");

endfunction

////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : iSTFT window plot functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : waterfall plot function
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function odmkWaterfall(zArrayIn,nfft,fs)
    //odmk waterfall FFT 3D plot
    //
    //zArrayIn -> expect rows = nBin; columns = # of time slices
    //currently symmetric matrix, nBin=nColumns = Nfft/2+1??
    
    //nBin = length(zArrayIn(1,:));
    nBin = nfft/2+1;
    //nSlice = length(zArrayIn(:,1));
    nSlice = length(zArrayIn(1,:));


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //begin: data formatting
    ///////////////////////////////////////////////////////////////////////////////////////////////
   
    //remove zeros from array and calculate magnitude and phase
    for e=1:nSlice
        for f=1:nBin
            if (zArrayIn(f,e) == 0) then
                zArrayInMAG(f,e) = 0.00001;
                zArrayInPHASE(f,e) = 0.00001;
            else    
                zArrayInMAG(f,e) = abs(zArrayIn(f,e));
                zArrayInPHASE(f,e) = atan(imag(zArrayIn(f,e)),real(zArrayIn(f,e)));
            end
        end
    end
    //disp(zArrayInMAG)
    
    //create a symmetrical arrays from lower bins & conj of lower bins (fft symm property)
    for k=1:nSlice
        for h=1:nfft/2-1
            zArrayConjMAG(h,k) = zArrayInMAG(nBin-h,k);
            zArrayConjPHASE(h,k) = zArrayInPHASE(nBin-h,k);
        end
    end    
    zArraySymmMAG = [zArrayInMAG(2:nBin,:);zArrayConjMAG];
    zArrayLbinMAG = [zArrayInMAG(2:nBin,:)];
    zArraySymmPHASE = [zArrayInPHASE(2:nBin,:);zArrayConjPHASE];
    
    
    //create linear spaces for 3D plotting functions
    //define x axis as frequency - number of freq. bins nBin
    //define y axis as time-slice - mumber of time-slices
    x = linspace(1,nBin,nBin);        //x = 1 x nBin array
    y = linspace(1,nSlice,nSlice);    //y = 1 x nSlice
    
    
    freqAxis = fs.*(1:length(zArrayInMAG(:,1)))/length(zArrayInMAG(:,1));
    //create log frequency axis, scaled by 1/100 (0 - 22050 Hz / 100)
    freqAxis_log = ((fs/2)/(max(log(x)))*log(x))./100;
    
    //xx = array of vectors of linear increasing integers (FFT bin indexii)
    //yy = array of vectors of constant y offsets (spacing for waterfall plots)
    xx = zeros(length(x),length(y));    //xx = 2D array nBins x nSlices
    yy = zeros(length(x),length(y));    //yy = 2D array nBins x nSlices
    for m=1:nSlice
        xx(:,m) = x(1:nBin)';
        xxlog(:,m) = freqAxis_log';
        yy(:,m) = (20*m).*ones(nBin,1);
    end


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //end: data formatting
    ///////////////////////////////////////////////////////////////////////////////////////////////


    //plot few slices
   
    scf(4)
    for i=1:5
        plot2d(freqAxis,[20*log10(zArrayInMAG(:,i))],style=[i]); 
    end
    title("$zArrayInMAG STFT Magnitude, first 5 frames$","color","black","fontsize",3);
    xlabel("$Linear Frequency : 0\le x\le fs (Sample Rate)$","fontsize",3,"color","black");
    ylabel("$Log-Scaled Magnitude$","fontsize",3,"color","black"); 
    a=gca(); // Handle on axes entity
    a.grid=[2,2];


    scf(5)
    for i=1:5
        plot2d(freqAxis_log,[20*log10(zArrayInMAG(:,i))],style=[i]); 
    end
    title("$zArrayInMAG STFT Magnitude, first 5 frames$","color","black","fontsize",3);
    xlabel("$Linear Frequency : 0\le x\le fs (Sample Rate)$","fontsize",3,"color","black");
    ylabel("$Log-Scaled Magnitude$","fontsize",3,"color","black"); 
    a=gca(); // Handle on axes entity
    a.grid=[2,2];
    

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //plot a 3D suface using plot3d
    ///////////////////////////////////////////////////////////////////////////////////////////////


    //use plot3d to plot surface spectrum (Magnitude)
    scf(55)
    plot3d(freqAxis_log,y,zArrayInMAG,-23,69, flag=[2 6 4]);
    f=gcf();           //get figure handle
    f.color_map = rainbowcolormap(32);    //set color table
    title 'waterfall STFT plot (1 channel)' 
    xlabel 'x-axis = frequency bins scaled by 1/100 (rows of sfftArray)'
    ylabel 'y-axis = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'



    //use plot3d to plot surface spectrum (log-scale Magnitude)
    scf(56)
    plot3d(freqAxis_log,y,20*log10(zArrayInMAG),-23,69, flag=[2 6 4]);
    f=gcf();           //get figure handle
    f.color_map = rainbowcolormap(32);    //set color table
    title 'waterfall STFT surface plot (1 channel)' 
    xlabel 'x-axis = frequency bins scaled by 1/100 (rows of sfftArray)'
    ylabel 'y-axis = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //linear waterfall plots
    ///////////////////////////////////////////////////////////////////////////////////////////////


    //use param3d1 to plot waterfall spectrum
    
    //create 3D coordinate vectors:
    //x = frequency axis, 0 - Nfft/2
    //y = constant offset spacing for each slice FFT magnitude array 
    //create array of color id values (what is scilab max?)    
    scf(7)    
    colors = 20*ones(length(xx(1,:)),1);
    //param3d1(xx,yy,list(zArrayInMAGxx,colors),-46,86,"X@Y@Z",[2,3])
    param3d1(xxlog,yy,list(zArrayInMAG,colors),-46,86,"X@Y@Z",[2,3])   
    title 'waterfall STFT array plot - linear Mag (1 channel)' 
    xlabel 'x-axis = frequency bins (rows of sfftArray)'
    ylabel 'y-axis = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'
     
    scf(8)    
    //param3d1(x,y,list(z,colors),[theta,alpha,leg,flag,ebox])
    colors = 19*ones(length(xx(1,:)),1);
    param3d1(xxlog,yy,list(20*log10(zArrayInMAG),colors),-46,86,"X@Y@Z",[2,3])   
    title 'waterfall STFT array plot - 20log10*Mag (1 channel)' 
    xlabel 'x-axis = frequency bins (rows of sfftArray)'
    ylabel 'y-axis = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'
 

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //radial waterfall plots
    /////////////////////////////////////////////////////////////////////////////////////////////// 
 
 
    //use param3d1 to plot radial waterfall spectrum (# bins = NFFT/2 - lower bins) 
    //generate radial "x" axi (Nfft/2 length spokes in the x,y plane)
    zn1 = cyclicZn(nSlice);

    //define color map - temp - increment/modulo colors    
    colors2=4+modulo(linspace(1,nSlice,nSlice),30);
    
    //create an array of baseVector - radial unit vectors:
    for b=1:nSlice
        //baseVector(:,b) = [0+%i*0;zn(b)];
        for c=1:nBin
            baseXY(c,b) = c*[zn1(b)];
        end
    end


    for j=1:nSlice
        //colors(j) = j;
        for k=1:nBin
            //create arrays of scaled x, constant y 
            xZn(j,k) = real(baseXY(k,j));
            yZn(j,k) = imag(baseXY(k,j));
        end
        //param3d(xx,yy,zArrayInx,[70,30,"X@Y@Z",[2,3]])        
    end    
    
    scf(69)    
    plot(real(baseXY),imag(baseXY))
    param3d1(xZn',yZn',list(zArrayInMAG,colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - linear Mag - nff2/2+1 lower bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'


    scf(70)    
    plot(real(baseXY),imag(baseXY))
    param3d1(xZn',yZn',list(20*log10(zArrayInMAG),colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - 20log10*Mag - nff2/2+1 lower bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'

    
    //create an array of baseVector - radial unit vectors:
    for b=1:nSlice
        //baseVector(:,b) = [0+%i*0;zn(b)];
        for c=1:length(zArraySymmMAG(:,1))
            base2XY(c,b) = c*[zn1(b)];
        end
    end
    //disp(base2XY)
    
    for j=1:nSlice
        //colors(j) = j;
        for k=1:length(zArraySymmMAG(:,1))
            //create arrays of scaled x, constant y 
            xxZn(k,j) = real(base2XY(k,j));
            yyZn(k,j) = imag(base2XY(k,j));
        end       
    end
    
    scf(71)    
    plot(real(base2XY),imag(base2XY))
    param3d1(xxZn,yyZn,list(zArraySymmMAG,colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - linear Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'
    

    scf(72)    
    plot(real(base2XY),imag(base2XY))
    param3d1(xxZn,yyZn,list(20*log10(zArraySymmMAG),colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - 20log10*Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'
    

    scf(73)    
    plot(real(base2XY),imag(base2XY))
    param3d1(xxZn,yyZn,list(fftshift(zArraySymmMAG),colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - 20log10*Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'    


    scf(74)    
    plot(real(base2XY),imag(base2XY))
    param3d1(xxZn,yyZn,list(20*log10(fftshift(zArraySymmMAG)),colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - 20log10*Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'


///////////////////////////////////////////////////////////////////////////////////////////////
//polar plots
///////////////////////////////////////////////////////////////////////////////////////////////
 
//    //polar plotting
//
//    x = 0:2*%pi/(length(zArraySymmMAG(:,1))-1):2*%pi;
//
//    scf(20)
//    //polarplot(x,odmkOscMAGwin(:,5),style=[23])
//    //polarplot(x,zArrayInMAG(2:1+nfft/2,1),style=[23])
//    polarplot(x,zArraySymmMAG(:,1),style=[23])
//    legend('y = odmkOscPHASEwin(:,3)',4)
//
//
//    scf(21)
//    //polarplot(x,odmkOscPHASEwin(:,5),style=[23])
//    //polarplot(x,zArrayInPHASE(2:1+nfft/2,1),style=[23])
//    polarplot(x,zArraySymmPHASE(:,1),style=[23])
//    legend('y = odmkOscPHASEwin(:,3)',4)
//
//    
//    //polar plot of composite frames 1-16 (arbitrary)
//    scf(22)
//        for i=1:15
//            //polarplot(x,zArrayInMAG(2:1+nfft/2,i),style=[i]);
//            polarplot(x,zArraySymmMAG(:,i),style=[i]);
//        end
//    legend('y = odmkOscPHASEwin(:,3)',4)
//
//
//    scf(23)
//        for i=1:15
//            //polarplot(x,zArrayInPHASE(2:1+nfft/2,i),style=[i]);
//            polarplot(x,zArraySymmPHASE(:,i),style=[i]); 
//        end
//    legend('y = odmkOscPHASEwin(:,3)',4)

    
    
endfunction


////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : waterfall plot function
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

mprintf('\n// *------------------------------------------------------------------------* //\n');
mprintf('      Begin Main')
mprintf('\n// *------------------------------------------------------------------------* //\n');

////////////////////////////////////////////////////////////////
//begin: osc
////////////////////////////////////////////////////////////////

//oscillator outputs selection of waveforms
//uses wavetable with linear interpolation


//DDS implementation outputs sin, cos. 
//fo = output frequency
//fs = sampling frequency
//td = table depth

//fo = fs/(2^N)*(delta(acc))
//delta(acc) = floor(fo*((2^N)/fs)+0.5)

fs = 44100.0;
mprintf('\nSampling frequency (fs) = %f\n',fs)



//number of steps through the algorithm -> determines length of output sample
nItera = 2048;    //2048 4096 8192 16384 32768 65536 131072 262144
mprintf('\nMnumber of main program iterations (nItera) = %i\n\n', nItera)


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

mprintf('\n// *------------------------------------------------------------------------* //\n'); 


////////////////////////////////////////////////////////////////
// begin: odmk_osc5stereo
////////////////////////////////////////////////////////////////


//--<<stereo>>--//

//output frequency: fo = fs/(2^N)*(phacc))
//foArray = 528;    //440 528 Hz

//foFixStereo = [528, 528];
foFixStereo = [7777, 7777];

//continually variable Frequency array:
//foArrayStereo = zeros(numSamples,2);
//foArrayStereo(:,1) = linspace(930,130,numSamples)'
//foArrayStereo(:,2) = linspace(130,930,numSamples)'


//define phase value (0-2pi radians)
//phase = [0,%pi/4,%pi/2,3*%pi/4,%pi];

phFixStereo = [0, 0];
//phFixStereo = [%pi/4, 3*%pi/4];

//continually variable Phase array:
//phArrayStereo = zeros(numSamples,2);
//foArrayStereo(:,1) = linspace(930,130,numSamples)'
//foArrayStereo(:,2) = linspace(130,930,numSamples)'

mprintf('\ninput source frequency - foFixStereo:\n')
disp(foFixStereo)    //displays hztext

mprintf('\ninput source phase offset - phFixStereo:\n')
disp(phFixStereo)    //displays hztext

nChan = 1;


//init wave file outputs
odmkOscStereo = zeros(nItera,2);
odmkOsc90Stereo = zeros(nItera,2);
odmkSqrPulseStereo = zeros(nItera,2);

[odmkOscStereo,odmkOsc90Stereo,odmkSqrPulseStereo] = odmk_osc5Stereo(nItera, fs, td, shape, foFixStereo, phFixStereo) 


//create mono signal
odmkOscStereo_L = odmkOscStereo(:,1);
//odmkOsc90Stereo_L = odmkOsc90Stereo(:,1);
//odmkSqrPulseStereo_L = odmkSqrPulseStereo(:,1);

pvInputSignal = odmkOscStereo_L
 
 
//////////////////////////////////////////////////////////////////
//// begin: odmk_multiOsc
////////////////////////////////////////////////////////////////// 
//
//output frequency: fo = fs/(2^N)*(phacc))
//create an array of frequency objects - fo_n = [fo_1,fo_2,...,fo_n] - assume Hertz
//fObj = [77,254,440,528,666];
//fObj = [1666,3666,5666,7666,9666];

//foArray = [666,2048,2333,3000,4096,7770,8192,11000,16666];
//phArray = [0, 0, 0, 0, 0, 0, 0, 0, 0];

//foArray = [7770];
//phArray = [0];
//
////phArray = [0, %pi/3, 2*%pi/3, %pi, 4*%pi/3, 5*%pi/3];
//
//nChan = length(foArray);
//
//mprintf('\nArray of frequencies for additive synthesis input source:\n')
//disp(foArray)    //displays hztext
//
//mprintf('\nArray of phase offsets for additive synthesis input source:\n')
//disp(phArray)    //displays hztext


////initialize variables,vectors 
//tb = zeros(td,1);    //table created by calling function tablegen
//
//
////initialize oscillator wavetable output
//odmkOsc = zeros(nItera,nChan);
////odmkosc2 = zeros(nItera,1);
////initialize oscillator square pulse output
//odmkSqrPulse = zeros(nItera,nChan);
////odmksqrpulse2 = zeros(nItera,1);
//
////initialize Sum of nObj waveforms
//odmkSum = zeros(nItera,1);
//odmkSumObj = zeros(nChan,1);
////initialize Mix of 2 waveforms
//odmkMix1 = zeros(nItera,1);
//
////init analysis of osc output
////if nItera > Nfft then
////    odmkOscFFT = zeros(Nfft,nObj);
////    odmkSqrPulseFFT = zeros(Nfft,nObj);
////    odmkOscMAG = zeros(Nfft,nObj);
////    odmkSqrPulseMag = zeros(Nfft,nObj);
////else     
////    odmkOsc1FFT = zeros(nItera,nObj);
////    odmkSqrPulseFFT = zeros(nItera,nObj);
////    odmkOscMAG = zeros(nItera,nObj);
////    odmkSqrPulseMag = zeros(nItera,nObj);
////end
//
////init wave file outputs 
//odmkOscWav = zeros(nItera,2*nChan);
//odmkSqrPulse = zeros(nItera,2*nChan);
//
//odmkOsc1_out = zeros(nItera,2);
//odmkOsc2_out = zeros(nItera,2);
//
//odmkSumWav = zeros(nItera,2);
//odmkMix1Wav = zeros(nItera,2);
//
////test taps
//acc_addr_tap = zeros(nItera,1);
//qnt_addr_tap = zeros(nItera,1);
//
//
//
//// main osc function call
//[odmkOsc, odmkOsc90, odmkSqrPulse] = odmk_multiOsc(nItera, fs, td, shape, foArray, phArray)
//
//mprintf('\nGenerating input signal (odmk_multiOsc func call)\n')
//
////Sum and mix ouput waveforms
//for m=1:nItera
//    //frequency sum
//    //accumulate row of osc values
//    odmkSumObj=0;
//    for n=1:nChan
//        odmkSumObj = odmkSumObj+odmkOsc(m,n);
//        //odmkSumSx(m) = odmkSumSx(m)+odmkSumObj(n);
//    end
//    odmkSum(m) = odmkSumObj;
//    
//    //frequency mix - only mix first two signals for 'sum/diff'
//    odmkMixIn1(m) = 0.5*(odmkOsc(m,1)+odmkOsc(m,2));    //??0.5 used to match FFT Mag outputs??
//    odmkMix1(m) = odmkOsc(m,1)*odmkOsc(m,2);
//end
////random waveforms:
//pvInputSignal = (1/max(abs(odmkSum))).*odmkSum    //normalize odmkSum
//odmkOsc1 = odmkOsc(:,1);


mprintf('\nPV input signal name: pvInputSignal\n')
mprintf('\nPV input signal length: %i', length(pvInputSignal))
mprintf('\nPV input signal nChan (# oscs): %i\n', nChan)

scf(111)    
plot(pvInputSignal,'blue')
title("$PV Input Signal (n-channel mixed sine)$","color","black","fontsize",3);

mprintf('\n// *------------------------------------------------------------------------* //\n');

////////////////////////////////////////////////////////////////
// end: osc
////////////////////////////////////////////////////////////////


pvNFFT = 128;
aWin = pvNFFT;
sWin = 0;           //frames are hann-windowed at nfft points when nWin is integer/=0, or rectangular if W=0  
Ra = .25;        //analysis hop factor, %of time-slice
Rs = .25;        //synthesis hop factor, %of time-slice
zp = 0;            //zero-padding (nfft -> 2*nfft, increases freq resolution)
specPlots = 0;     //3D waterfall plots of spectral data
noSpecPlots = 0;

//1024 samples is about 60 ms at 16kHz, a good window
//wavVocoded1=phaseVoc(odmkOsc1,Ra_hop,Rs_hop,pvNFFT,aWin,sWin,zp,specPlots);
//wavVocodedR=phaseVoc(wavDataR,Ra_hop,Rs_hop,pvNFFT,pvNFFT,noSpecPlots);
// Compare original and resynthesis
//

wavVocoded1=phaseVoc(pvInputSignal,Ra,Rs,pvNFFT,aWin,sWin,zp,specPlots);


mprintf('PV FFT length: %i\n', pvNFFT)
mprintf('PV Analysis Hop Length (Ra_hop) = %f\n', fix(Ra*pvNFFT))
mprintf('PV Synthesis Hop Length (Rs_hop) = %f\n', fix(Rs*pvNFFT))
mprintf('PV Zero padding (0=off, 1=on) = %i\n', zp)

mprintf('\n// *------------------------------------------------------------------------* //\n');

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

    

///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : plotting
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////