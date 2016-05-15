///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : header
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////
//
//__<<name=> "phaseVocoder_transcendental" (.sce)>>__

//___::((ODMK:RHOTAH:ODERUS MONK:2012))::___


//include components:

//___::((phaseVocoder2))::___
//___::((24 but .wav file read & write))::___
//

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//This script frequency shifting using a phase vocoder, overlap FFT algorithm

//Real time controls:
//step size (frequency/velocity); 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : header
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


xdel(winsid()) //-> closes all open graphs (= matlab.close('all'))
clear;
clc;


///////////////////////////////////////////////////////////////////////////////////////////////
//begin : stack adjustment
///////////////////////////////////////////////////////////////////////////////////////////////

//try for ~10MBytes
//stacksize input is doubles
//
max_bytes = 20*10e7;
max_bits=24;
//bits=24;
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


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : cyclicZn functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//generate cyclic group Zn member elements
function [cZn] = cyclicZn(n)
    
    cZn = zeros(n,1);
    for k=1:n
        //z(k) = %e^(((k-1)*2*%pi*%i)/n)                          //Define cyclic group Zn points
        cZn(k) = cos(((k-1)*2*%pi)/n)+%i*sin(((k-1)*2*%pi)/n)   //Euler's identity    
    end
endfunction

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : cyclicZn functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : phase vocoder functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

////////////////////////////////////////////////////////////////
//function: phaseVoc
//phase vocoder top-level function - Uses pvMod, STFT, iSTFT

function wavDataOut = phaseVoc(wavDataIn, Ra, Rs, nfft, aWin, sWin, zp, fs, pkDet, maxPeaks, fplots)
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

    wavInLength = length(wavDataIn);    //total length of input audio sample - analysis 
    mprintf('phaseVoc.wavInLength = length of wavIn = %i\n',wavInLength);

    // Calculate the basic STFT, magnitude scaled
    wavX = STFT(wavDataIn', nfft, aWin, Ra, zp, fs);    //check if iSTFT needs to change with zp?
    //wavX = STFT(wavDataIn', nfft, 0, Ra, zp, fs);
    
    //***temp***
//    //plot out first 8 results of first 4 frames of wavX
//    for ii=1:4
//        for jj=1:8
//            wavXtemp(jj,ii) = wavX(jj,ii);
//        end
//    end
//    disp(wavXtemp)
    
    [nSTFTRow, nSTFTCol] = size(wavX);
    //nRow = number of frequency bins = nfft/2 + 1
    //nCol = number of audio time-slices (windowed and fft)
    
    mprintf('phaseVoc.STFT.nSTFTRow = number of rows of STFT output matrix (nBins) = %i\n',nSTFTRow);
    mprintf('phaseVoc.STFT.nSTFTCol = number of columns of STFT output matrix (nSlice)= %i\n',nSTFTCol);
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //time-scale modification
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //modification factor
    Rmod = Ra/Rs;

    //Calculate the new timebase samples
    t = 0:Rmod:(nSTFTCol-2);
    // Have to stay two columns off end because (a) counting from zero, and 
    // (b) need col n AND col n+1 to interpolate

    // Generate the new spectrogram
    pvocX = pvMod(wavX, t, Rmod);
    
    //***temp***
//    //plot out first 8 results of first 4 frames of pvocX
//    for ii=1:4
//        for jj=1:8
//            pvocXtemp(jj,ii) = pvocX(jj,ii);
//        end
//    end
//    disp(pvocXtemp)
    
    
    if pkDet==1 then
        [peakLocArray,peakArray] = pvPeaks(wavX, nfft, maxPeaks)
    end
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //synthesis
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // Invert to a waveform
    wavDataOut = iSTFT(pvocX, nfft, sWin, Rs, zp)';
    //wavDataOut = iSTFT(wavX, nfft, sWin, Rs, zp)';
    
    //crop wavDataOut to match length of wavDataIn
    //wavDataOut = wavDataOut(1:length(wavDataIn));
    
    wavOutLength = length(wavDataOut);    //total length of output audio sample - Synthesis
    mprintf('phaseVoc.wavOutLength = number of columns of iSTFT output array = %i\n',wavOutLength);    
    
    //STFT plotting
    if fplots==1 then
        if zp ==0 then
            nfftx = nfft;
        else
            nfftx = nfft*2;
        end
        odmkWaterfall(wavX, nfftx, fs)
        
//        if pkDet==1 then
//        //call pvPeaks plotting function
//        odmkpvPeaksPlot(wavX, peakLocArray, peakArray, nfft, fs)          
//           
//        end   
            
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //check results
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        scf(22)    
            plot(wavDataIn,'blue')
            plot(wavDataOut,'red')
            title("STFT - iSTFT compare - wavDataIn - blue; wavDataOut - red","color","black","fontsize",3);                
    
        //wavError = wavDataIn-wavDataOut;
        //scf(23)  
            //plot(wavError)
            //title("STFT - iSTFT error - wavDataIn-wavDataOut","color","black","fontsize",3);
    
    end

    if pkDet==1 then
        //call pvPeaks plotting function
        odmkpvPeaksPlot(wavX, peakLocArray, peakArray, nfft, fs)          
    end 

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
    //zero pad wavIn to take care of last frame
    //adds enough zeros so that wav length is evenly divided by overlapped windows
    //no trailing samples left unprocessed
    //in case when wavIn is already an even muliple of overlapped windows,
    //nothing should be added, aLength=wevInLength 
    aLength = wavInLength + modulo(fix(wavInLength/((1-Ra)*nfft))*((1-Ra)*nfft)+((1-Ra)*nfft),wavInLength);
    mprintf('phaseVoc.STFT.aLength = length of internal wav (zero-padded wav)= %i\n',aLength);
    
    wavInZP = {wavIn,zeros(1,(aLength-wavInLength))};
    
    //odmkSTFTWinPlot(aWindow,Ra,aLength)
    odmkSTFTWinPlot(aWindow,Ra,length(aWindow)*5)

    //perform windowing and fft on slices of audio data
    nCol = 1;
    
   //check if zp is enabled
    if zp == 0 then
        //no frame zero padding
        // pre-allocate output array - exploit symmetry - only save lower bins + 1 (center bin)
        stftMatrix = zeros((1+nfft/2),1+fix((aLength-nfft)/(Ra*nfft)));
        //b runs length of sample minus nfft 
        for b = 0:(Ra*nfft):(aLength-nfft)        
            //create temp matrix u = window dotMPY with current audio time-slice
            u = aWindow.*wavInZP((b+1):(b+nfft));
            X = fft(u);
            //truncate top half of freq bins
            stftMatrix(:,nCol) = X(1:(1+nfft/2))';
            nCol = nCol+1;
        end;
    else
        //frame zero padding
        // pre-allocate output array - exploit symmetry - only save lower bins + 1 (center bin)
        stftMatrix = zeros((1+nfft),1+fix((aLength-nfft)/(Ra*nfft)));
        //b runs length of sample minus nfft 
        for b = 0:(Ra*nfft):(aLength-nfft)        
            //create temp matrix u = window dotMPY with current audio time-slice
            //zero padding           
            u = aWindow.*wavInZP((b+1):(b+nfft));
            X = fft([u,zeros(1,nfft)]);
            //truncate top half of freq bins (length of X actually 2*nfft)
            stftMatrix(:,nCol) = X(1:(1+nfft))';
            nCol = nCol+1;
        end;
    end;

    stftMatrixOut = stftMatrix;
   
    
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


//    //plot window information - moved to separate function
//    ttw = (1:1:length(win));
//    scf(1)
//    plot(ttw,win,'red')
//    title("window function","color","black","fontsize",3);
//    xlabel("Time axis : 0\le x\le window length (Nfft)$","fontsize",4,"color","black");
//    ylabel("Window Magnitude","fontsize",3,"color","black");
    
        
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
            sWindow = ones(1,nfft);
        else
            halfLength = floor(sWin/2);   // midpoint of win
            //generate 1/2 of a cos window function
            halfwin = 0.5*(1 + cos(%pi*(0:halfLength)/halfLength));
            sWindow = zeros(1,sWin);
            //fill in top-half of windoe array
            sWindow((halfLength+1):(halfLength+halfLength)) = halfwin(1:halfLength);
            //fill in bottom-half of windoe array
            sWindow((halfLength+1):-1:2) = halfwin(1:halfLength);
        end
    else
        sWindow = sWin;    //in case of explicit window array passed to funcfrun
    end

    sLength = nfft + (niCol-1)*(Rs*nfft);
    mprintf('phaseVoc.STFT.sLength = length of internal wavOut (zero-padded)= %i\n',sLength);
        
    wavOut = zeros(1,sLength);
//    //b runs length of sample minus nfft
//    for b = 0:Rs:(Rs*(cols-1))
//        disp(b)
//        v = istftMatrixIn(:,1+b/Rs)';
//        //create full frame from conjugate of lower bins (fft symmetry)
//        v = [v, conj(v([((nfft/2)):-1:2]))];
//        F = real(ifft(v));
//        //overlap-add windowed ifft frames
//        wavOut((b+1):(b+nfft)) = wavOut((b+1):(b+nfft))+F.*sWindow;
//    end;
    
        //b runs length of sample minus nfft
    for b = 0:niCol-1
        v = istftMatrixIn(:,b+1)';
        //create full frame from conjugate of lower bins (fft symmetry)
        if zp == 0 then
            v = [v, conj(v([((nfft/2)):-1:2]))];
            F = real(ifft(v));
        else
            v = [v, conj(v([((nfft)):-1:2]))];
            Fzp = real(ifft(v));
            F = Fzp(1:nfft);
        end
        //overlap-add windowed ifft frames
        wavOut((b*(Rs*nfft)+1):(b*(Rs*nfft)+nfft)) = wavOut((b*(Rs*nfft)+1):(b*(Rs*nfft)+nfft))+F.*sWindow;
    end;
    
    //normalize gain (**window overlap > 50% causes linear increase in gain)
    istftGain = 1/(0.5/Rs);
    wavOut = istftGain*wavOut;
    
    
    //sLength = length(wavOut);    //total length of audio - synthesis
    odmkiSTFTWinPlot(sWindow,Rs,length(sWindow)*5)
    
endfunction



function pvArray = pvMod(stftArrayIn, t, nfft, Ra)
    //     pvArray = pvMod(stftArrayIn, t, nfft, Ra)  time-scale modifications: time-stretching / time-compression  
    //     Interpolate an STFT array according to the 'phase vocoder' standard method using phase unwrapping
    //     stftArrayIn is a 2D STFT array, with nfft/2+1 rows (dft freq bins), and n colums (time slices, or frames).
    //     t is a vector of (real) time-samples, which specifies a path through 
    //     the time-base defined by the columns of stftArrayIn.  For each value of t, 
    //     the spectral magnitudes in the columns of stftArrayIn are interpolated, and 
    //     the phase difference between the successive columns of stftArrayIn is 
    //     calculated; a new column is created in the output array pvArray that 
    //     preserves this per-step phase advance in each bin.
    //     Ra is the STFT hop size, defaults to Nfft/2
    //     Ra is needed to calculate the 'null' phase advance expected in each bin.
    //     Note: t is defined relative to a zero origin, so 0.1 is 90% of 
    //     the first column of stftArrayIn, plus 10% of the second.

    //if nargin < 3
        //Ra = 0;
    //end

    [nBin,nSlice] = size(stftArrayIn);

    nfft2 = nfft/2;

    if Ra == 0
        // default value
        Ra = nfft2;
    end

    // Empty output array
    pvArray = zeros(nBin, length(t));

    // Expected phase advance in each bin
    delPhi = zeros(1,nfft2+1);
    delPhi(2:(nfft2+1)) = (2*%pi*Ra)./(nfft./(1:(nfft2)));

    // Phase accumulator
    // Preset to phase of first frame for perfect reconstruction
    // in case of 1:1 time scaling
    //ph = angle(stftArrayIn(:,1));
    ph = atan(imag(stftArrayIn(:,1)),real(stftArrayIn(:,1)));

    // Append a 'safety' column on to the end of stftArrayIn to avoid problems 
    // taking *exactly* the last frame (i.e. 1*stftArrayIn(:,cols)+0*stftArrayIn(:,cols+1))
    stftArrayIn = [stftArrayIn,zeros(nBin,1)];

    nCol = 1;
    for tt = t
        // Grab the two columns of b
        stftArray_cols = stftArrayIn(:,floor(tt)+[1 2]);
        tf = tt - floor(tt);
        stftArray_mag = (1-tf)*abs(stftArray_cols(:,1)) + tf*(abs(stftArray_cols(:,2)));
        // calculate phase advance
        //dp = angle(stftArray_cols(:,2)) - angle(stftArray_cols(:,1)) - delPhi';
        dp = atan(imag(stftArray_cols(:,2)),real(stftArray_cols(:,2))) - atan(imag(stftArray_cols(:,1)),real(stftArray_cols(:,1))) - delPhi';
        // Reduce to -pi:pi range
        dp = dp - 2 * %pi * round(dp/(2*%pi));
        // Save the column
        pvArray(:,nCol) = stftArray_mag .* exp(%i*ph);
        // Cumulate phase, ready for next frame
        ph = ph + delPhi' + dp;
        nCol = nCol+1;
    end
endfunction


function [peakLocArray,peakArray] = pvPeaks(stftArrayIn, nfft, maxPeaks)
    //     pvArray = pvPeaks(stftArrayIn, nfft, maxPeaks)  peak detection & define regions
    //     stftArrayIn is a 2D STFT array of DFT Magnitudes with rows = freq. bins, cols = overlapped DFT output frames
    //     pvPeaks outputs a 2D array of peak locations, and regions bounaries
    //     the peak locations are dft bin values corresponding to the dft Magnitudes of the detected peaks
    //     the region boundaries are the dft bin values corresponding to the mid-points between detected peaks
    //     the peakArray is initialized to values evenly distributed across the freq. spectrum then updated as peaks are found
    [nBin,nSlice] = size(stftArrayIn);
    nfft2 = nfft/2;
    //peakLocStep = fix(nfft2/maxPeaks);
    
    //initialize output array
    //maxPeaks*2-1 = alternating peak locations, + region boundaries
    peakLocArray = zeros(maxPeaks*2-1,nSlice);
    peakArray = zeros(maxPeaks,nSlice);
    
    //initilize sub array with evenly distributed peak locations
//    peakSubArrayInit(1) = peakLocStep;
//    for i = 1:maxPeaks-1
//        peakSubArrayInit(i+1) = peakSubArrayInit(i) + peakLocStep;
//    end

    //initilize sub array - one frame to temporarily store peaks
    peakSubArrayInit = zeros(maxPeaks,1);
    
    for j = 1:nSlice        //loop through time-slices
        //binCnt = 0;
        numPeaks = 0;
        //the output of the fft is scaled by Nfft/2
        minPeak = nfft2;        //set to max Magnitude so that first deteted peak will become initial minimum
        minPeakLoc = 0;    //location in subarray of minimum peak
        replacedLoc = 0;
        
        prevMag2 = nfft2;        //set to max Magnitude to prevent detection of DC **Change to direct calc
        prevMag1 = nfft2;        //set to max Magnitude to prevent detection of DC
        currentMag = 0;
        nextMag1 = 0;
        nextMag2 = 0;
        
        //initilize sub array with evenly distributed peak locations
        peakLocSubArray = zeros(maxPeaks,1);   //one column for temp storage of peak locations
        peakSubArray = zeros(maxPeaks,1);   //one column for temp storage of peak magnitudes
        regionSubArray = zeros(maxPeaks-1,1);
        
        for k = 1:nBin        //loop through DFT freq. bins
            SubArraylength = length(peakSubArray);
            peakFound = 0;

            prevMag2 = prevMag1;
            prevMag1 = currentMag;
            currentMag = nextMag1;
            nextMag1 = nextMag2;
            nextMag2 = abs(stftArrayIn(k,j));
            
            //disp(prevMag2)
            //disp(prevMag1)
            //disp(currentMag)
            //disp(nextMag1)
            //disp(stftArrayIn(k,j))
            
            //detect if peak
            //we want to build an array of the largest detected peaks, length of array is maxPeaks 
            //if a new peak is detected that is greater than the current minimum peak, then discard old minimum peak
            if ((prevMag2 < currentMag) & (prevMag1 < currentMag) & (nextMag1 < currentMag) & (nextMag2 < currentMag)) then
                numPeaks = numPeaks +1;
                if (numPeaks <= maxPeaks) then
                    //peak detected - write k-th bin#(-2 delay compensation) into slot of peakLocSubArray
                    peakLocSubArray(numPeaks) = k-2;
                    //peak detected - write k-th bin Mag into slot of peakSubArray    
                    peakSubArray(numPeaks) = currentMag;    
                    //keep track of minimum peak value & location
                    if (currentMag < minPeak) then
                        minPeakLoc = numPeaks;
                        minPeak = currentMag;
                    end
                //***change to replace based on dynamic range of current bin vs. surrounding bins
                //update peakSubArray  & minPeakLoc with current bin location
                elseif ((numPeaks > maxPeaks) & (currentMag > minPeak)) then    
                //subArray full, so start replacement routine
                    minPeak = currentMag;                 //temporarily initialize minPeak to currentMag
                    replacedLoc = minPeakLoc;             //store temporary minpeak location
                    minPeakLoc = SubArraylength;    //temporarily set minPeakLoc to last position of array
                    for p = 1:(SubArraylength-1)
                        //shift array locations to accomodate a new detected peak
                        if (p >= replacedLoc) then
                            peakLocSubArray(p) = peakLocSubArray(p+1);
                            peakSubArray(p) = peakSubArray(p+1);    
                        end
                        if (peakSubArray(p) < minPeak) then
                            minPeakLoc = p;
                            minpeak = peakSubArray(p);
                        end
                    end
                    //add new detected peak to end of array
                    peakLocSubArray(SubArraylength) = k-2;
                    peakSubArray(SubArraylength) = currentMag;
                end
                //this point, should have a linearly ascending list of peak locations, and a final min peak location
            end    //end detect peak
        end    //end loop through DFT freq. bins
        //***temp***
        if j==5 then
            disp(peakSubArray,'5th slice peakSubArray values')
            disp(peakLocSubArray,'5th slice peakLocSubArray values')
        end
        //***end temp***
        
        //mprintf('phaseVoc.pvPeaks.numPeaks = number of peaks detected = %i\n',numPeaks);
        //calculate region boundaries between peak locations (**regionSubArray 1 less than peakLocSubArray)
        //peakLocArray contains alternating - peak - region boundary - peak ...
        for q = 1:length(regionSubArray)
            regionSubArray(q) = (peakLocSubArray(q+1)-peakLocSubArray(q))/2;
            peakLocArray((2*q),j) = regionSubArray(q); 
        end
        for r = 1:length(peakLocSubArray)
            peakLocArray((2*r-1),j) = peakLocSubArray(r);
        end
        peakArray(:,j) = peakSubArray;                
    end

    
    [npeakArrayRow, npeakArrayCol] = size(peakArray);
    //nRow = number of frequency bins = nfft/2 + 1
    //nCol = number of audio time-slices (windowed and fft)
    
    mprintf('phaseVoc.pvPeaks.npeakArrayRow = number of rows of peakArray matrix (nBins = maxPeaks) = %i\n',npeakArrayRow);
    mprintf('phaseVoc.pvPeaks.npeakArrayCol = number of columns of peakArray matrix (nSlice)= %i\n',npeakArrayCol);
    
    [npeakLocArrayRow, npeakLocArrayCol] = size(peakLocArray);
    //nRow = number of frequency bins = nfft/2 + 1
    //nCol = number of audio time-slices (windowed and fft)
    
    mprintf('phaseVoc.pvPeaks.npeakLocArrayRow = number of rows of peakLocArray matrix (nBins = 2 x maxPeaks -1) = %i\n',npeakLocArrayRow);
    mprintf('phaseVoc.pvPeaks.npeakLocArrayCol = number of columns of peakLocArray matrix (nSlice)= %i\n',npeakLocArrayCol);
    
    
    //transcendental - fixed initial peak array replace algo..            
//                   if (replacedLoc == 0) then
//                        if ((k - ((k/nfft2)*peakLocStep)) < ((k+1) - (((k+1)/nfft2)*peakLocStep))) then
//                            peakSubArray((k - ((k/nfft2)*peakLocStep))) = k;
//                        else
//                            peakSubArray(((k+1) - (((k+1)/nfft2)*peakLocStep))) = k;   
//                        end
//                    else
//                        peakSubArray(((k+1) - (((k+1)/nfft2)*peakLocStep))) = k;    
//                    end
     
endfunction


////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : phase vocoder functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : STFT window plot function
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function odmkSTFTWinPlot(win,Ra,aLength)
    
    //mprintf('phaseVoc.aLength = length of wavIn = %i\n',aLength);

    //plot window information
    tw = (1:1:length(win));
scf(1)
    plot(tw,win,'red')
    title("$STFT Analysis Window$","color","black","fontsize",3);
    xlabel("$Time axis : 0\le x\le window length (Nfft)$","fontsize",3,"color","black");
    ylabel("$Window Amplitude$","fontsize",3,"color","black");

 
    //plot a window for each audio time-slice
    tt = (1:1:aLength);
scf(2)    
    winPlotSum = zeros(1,aLength);
    for i = 0:(Ra*nfft):(aLength-nfft)      
        winPlotArray = zeros(1,aLength);
        winPlotArray((i+1):(i+nfft)) = win;
        winPlotSum = winPlotSum+winPlotArray;
        plot2d(tt,[winPlotArray],style=[5]); 
    end;
    plot2d(tt,winPlotSum,style=[1]);
    title("$STFT overlapped Analysis Windows (red) ; overlapped sum (black)$","color","black","fontsize",3);
    xlabel("$Time axis : 0\le x\le wav length$","fontsize",3,"color","black");
    ylabel("$Window Amplitude$","fontsize",3,"color","black"); 

//scf(3)
//    plot2d(tt,winPlotSum,style=[1]);
//    title("$STFT overlapped windows sum - time domain$","color","blue","fontsize",3);
//    xlabel("$Time axis : 0\le x\le wav length$","fontsize",3,"color","black");
//    ylabel("$Window Amplitude$","fontsize",3,"color","black");

endfunction

////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : STFT window plot functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : iSTFT window plot function
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function odmkiSTFTWinPlot(win,Rs,sLength)

    //plot window information
    tw = (1:1:length(win));
scf(19)
    plot(tw,win,'red')
        title("$iSTFT Synthesis Window$","color","black","fontsize",3);
    xlabel("$Time axis : 0\le x\le window length (Nfft)$","fontsize",3,"color","black");
    ylabel("$Window Amplitude$","fontsize",3,"color","black");

 
    //plot a window for each audio time-slice
    tt = (1:1:sLength);
scf(20)    
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

scf(21)
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
//begin : odmkpvPeaksPlot plot function
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function odmkpvPeaksPlot(zArrayIn,peakLocArray,peakArray,nfft,fs)
    
    
    //nBin = length(zArrayIn(1,:));
    nBin = nfft/2+1;
    //nSlice = length(zArrayIn(:,1));
    nSlice = length(zArrayIn(1,:));
    
    numPeaks = length(peakArray(:,1));    //same as maxPeaks
    
    
    //***temp***
    
    //plot one frame of peakArray
    tt = [1:1:length(peakArray(:,5))];
    
    scf(23)
    axfig23 = gca();
    plot(axfig23,[tt; tt],[zeros(1:length(peakArray(:,5))); peakArray(:,5)'],'color','red')
    title("$5th slice peakArray$","color","black","fontsize",3);    
    disp(peakArray(:,5),'5th slice peakArray values')
    a=get("current_axes");//get the handle of the newly created axes
    a.data_bounds=[0,0;10,70];
    
    ///plot one frame of peakLocArray
    scf(24)
    plot(peakLocArray(:,5));
    title("$5th slice peakLocArray$","color","black","fontsize",3);
    disp(peakLocArray(:,5),'5th slice peakLocArray values');

    
    //***end temp***
    for e=1:nSlice
        pk1Array = zeros(nBin,nSlice);
        j=1;
        LocIndex = 0;
        for g=1:nBin
            LocIndex = LocIndex+1;
            if (LocIndex == peakLocArray(j)) then
                pk1Array(e,LocIndex) = peakArray(e,j);
                j=j+1;           
            end
        end        
    end            
        
    //plot one frame of pk1Array
    scf(25)
    plot(pk1Array(:,5));
    title("$5th slice pk1Array$","color","black","fontsize",3);    
    disp(pk1Array(:,5),'5th slice pk1Array values')
    
//    ///plot one frame of peakLocArray
//    scf(24)
//    plot(peakLocArray(:,5));
//    title("$5th slice peakLocArray","color","black","fontsize",3);
//    disp(peakLocArray(:,5),'5th slice peakLocArray values');    
    
    //create a matrix for plotting overlayed peaks
    peakOverlayArray = zeros(nBin,nSlice);
    
    for c=1:nSlice
        peakLocInc = 1;
        for d=1:nBin
            if (peakLocArray(peakLocInc) == d) then
                peakOverlayArray(d,c) = peakArray(peakLocInc);
                peakLocInc = peakLocInc + 1;
            else
                peakOverlayArray(d,c) = 0;    
            end
        end
    end
    
    //create a matrix for plotting overlayed regions
    peakOverlayArray = zeros(nBin,nSlice);
    
    for c=1:nSlice
        peakLocInc = 1;
        for d=1:nBin
            if (peakLocArray(peakLocInc) == d) then
                peakOverlayArray(d,c) = peakArray(peakLocInc);
                peakLocInc = peakLocInc + 1;
            else
                peakOverlayArray(d,c) = 0;    
            end
        end
    end
    
    
    
    
    //nBin = length(zArrayIn(1,:));
    nBin = nfft/2+1;
    //nSlice = length(zArrayIn(:,1));
    nSlice = length(zArrayIn(1,:));


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //begin: data formatting
    ///////////////////////////////////////////////////////////////////////////////////////////////

    
    zArrayInMAG = abs(zArrayIn);
    zArrayInPHASE = atan(imag(zArrayIn),real(zArrayIn));
    
    //disp(zArrayInMAG)
    
    //create a symmetrical arrays from lower bins & conj of lower bins (fft symm property)
    for k=1:nSlice
        for h=1:nfft/2-1
            zArrayConjMAG(h,k) = zArrayInMAG(nBin-h,k);
            zArrayConjPHASE(h,k) = zArrayInPHASE(nBin-h,k);
        end
    end    
    zArraySymmMAG = [zArrayInMAG(2:nBin,:);zArrayConjMAG];
    zArraySymmPHASE = [zArrayInPHASE(2:nBin,:);zArrayConjPHASE];
    
    
    //create linear spaces for 3D plotting functions
    //define x axis as frequency - number of freq. bins nBin
    //define y axis as time-slice - mumber of time-slices
    x = linspace(1,nBin,nBin);        //x = 1 x nBin array
    y = linspace(1,nSlice,nSlice);    //y = 1 x nSlice
    

    freqAxis_log = ((fs/2)/(max(log(x)))*log(x));
    
    
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
    //linear waterfall plots
    ///////////////////////////////////////////////////////////////////////////////////////////////


    //use param3d1 to plot waterfall spectrum
    
    //create 3D coordinate vectors:
    //x = frequency axis, 0 - Nfft/2
    //y = constant offset spacing for each slice FFT magnitude array 
    //create array of color id values (what is scilab max?)    
//scf(7)    
//    colors = 6*ones(length(xx(1,:)),1);
//    //param3d1(xx,yy,list(zArrayInMAGxx,colors),-46,86,"X@Y@Z",[2,3])
//    param3d1(xxlog,yy,list(zArrayInMAG,colors),-46,86,"X@Y@Z",[2,3])   
//    title 'waterfall STFT array plot - linear Mag (1 channel)' 
//    xlabel 'x-axis = frequency bins (rows of sfftArray)'
//    ylabel 'y-axis = audio slices (columns of sfftArray)'
//    zlabel 'z-axis = frequency bin magnitude'
//     
//scf(8)    
//    //param3d1(x,y,list(z,colors),[theta,alpha,leg,flag,ebox])
//    colors = 7*ones(length(xx(1,:)),1);
//    param3d1(xxlog,yy,list(20*log10(zArrayInMAG),colors),-46,86,"X@Y@Z",[2,3])   
//    title 'waterfall STFT array plot - 20log10*Mag (1 channel)' 
//    xlabel 'x-axis = frequency bins (rows of sfftArray)'
//    ylabel 'y-axis = audio slices (columns of sfftArray)'
//    zlabel 'z-axis = frequency bin magnitude'
    
    
endfunction    
    
    
    

////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : odmkpvPeaksPlot plot functions
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
    
//old - creates a symmetrical nxn 2D array  
//    //define x axis as frequency - number of freq. bins nBin
//    //define y axis as time-slice - mumber of time-slices
//    if nSlice<=nBin then
//        //initialie truncated Mag array
//        zArrayInMAGxx = zeros(nSlice,nSlice);
//        for m=1:nSlice
//            xx(:,m) = x(1:nSlice)';
//            yy(:,m) = (4*m).*ones(nSlice,1);
//            zArrayInMAGxx(m,:) = zArrayInMAG(m,1:nSlice);
//        end
//    else
//        //initialie truncated Mag array
//        zArrayInMAGxx = zeros(nBin,nBin);
//        for n=1:nBin
//            xx(:,n) = y(1:nBin)';
//            yy(:,n) = (4*n).*ones(nBin,1);;
//            dlength = nBin;
//            zArrayInMAGxx(:,n) = zArrayInMAG(1:nBin,n);
//            //zArrayInPHASExx = zArrayInPHASE(1:nBin);
//        end
//    end


///////////////////////////////////////////////////////////////////////////////////////////////
//end: data formatting
///////////////////////////////////////////////////////////////////////////////////////////////

//    //create frequency axis vectors
//    freqAxis=fs*(0:(Nfft-1))/Nfft;
//
//    freqAxisx=fs*(1:(Nfft))/Nfft;
//    freqAxis_log = log(freqAxisx);
// 
//    freqAxis_zp=fs*(0:(2*Nfft-1))/(2*Nfft);
//
//    freqAxis_zpx=fs*(1:(2*Nfft))/Nfft;
//    freqAxis_zplog = log(freqAxis_zpx);
//
//    //freqAxis_zp=fs*(0:(2*Nfft-1))/Nfft;
//
//    if (Nfft > 100) then
//        shortFreqAxis = [1:1:100];
//    else
//        shortFreqAxis = [1:1:Nfft];
//    end


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


    //use plot3d to plot surface spectrum (Magnitude)
scf(5)
    plot3d(freqAxis_log,y,zArrayInMAG,-23,69, flag=[2 6 4]);
    f=gcf();           //get figure handle
    f.color_map = rainbowcolormap(32);    //set color table
    title 'waterfall STFT plot (1 channel)' 
    xlabel 'x-axis = frequency bins scaled by 1/100 (rows of sfftArray)'
    ylabel 'y-axis = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'


///////////////////////////////////////////////////////////////////////////////////////////////
//plot a 3D suface using plot3d
///////////////////////////////////////////////////////////////////////////////////////////////


    //use plot3d to plot surface spectrum (log-scale Magnitude)
scf(6)
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
    
    //create an array of baseVector - radial unit vectors:
    for b=1:nSlice
        //baseVector(:,b) = [0+%i*0;zn(b)];
        for c=1:nBin
            baseXY(c,b) = c*[zn1(b)];
        end
    end
scf(9)    
    plot(real(baseXY),imag(baseXY))
       
    colors2 = 21*ones(nSlice,1);
//    for j=1:nSlice
//        //colors(j) = j;
//        for k=1:nBin
//            //create arrays of scaled x, constant y 
//            xZn(j,k) = real(baseXY(k,j));
//            yZn(j,k) = imag(baseXY(k,j));
//        end
//        //param3d(xx,yy,zArrayInx,[70,30,"X@Y@Z",[2,3]])        
//    end
    
    for j=1:nSlice
        //colors(j) = j;
        for k=1:nBin
            //create arrays of scaled x, constant y 
            xZn(j,k) = real(baseXY(k,j));
            yZn(j,k) = imag(baseXY(k,j));
        end
        //param3d(xx,yy,zArrayInx,[70,30,"X@Y@Z",[2,3]])        
    end
    
    param3d1(xZn',yZn',list(zArrayInMAG,colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - linear Mag - nff2/2+1 lower bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'

scf(10)    
    plot(real(baseXY),imag(baseXY))
       
    colors2 = 22*ones(nSlice,1);
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
scf(11)    
    plot(real(base2XY),imag(base2XY))
       
    colors2 = 21*ones(nSlice,1);
    for j=1:nSlice
        //colors(j) = j;
        for k=1:length(zArraySymmMAG(:,1))
            //create arrays of scaled x, constant y 
            xxZn(k,j) = real(base2XY(k,j));
            yyZn(k,j) = imag(base2XY(k,j));
        end
        //param3d(xx,yy,zArrayInx,[70,30,"X@Y@Z",[2,3]])        
    end
    param3d1(xxZn,yyZn,list(zArraySymmMAG,colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - linear Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'
    

scf(12)    
    plot(real(base2XY),imag(base2XY))
       
    colors2 = 9*ones(nSlice,1);
    param3d1(xxZn,yyZn,list(20*log10(zArraySymmMAG),colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - 20log10*Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'
    

scf(13)    
    plot(real(base2XY),imag(base2XY))
       
    colors2 = 19*ones(nSlice,1);
    param3d1(xxZn,yyZn,list(fftshift(zArraySymmMAG),colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - 20log10*Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'    


scf(14)    
    plot(real(base2XY),imag(base2XY))
       
    colors2 = 18*ones(nSlice,1);
    param3d1(xxZn,yyZn,list(20*log10(fftshift(zArraySymmMAG)),colors2),-46,86,"X@Y@Z",[2,3])
    title 'radial waterfall STFT plot - 20log10*Mag - symmetrical bins (1 channel)' 
    xlabel 'radii = 0-nfft/2+1 frequency bins (rows of sfftArray)'
    ylabel 'angle = audio slices (columns of sfftArray)'
    zlabel 'z-axis = frequency bin magnitude'

///////////////////////////////////////////////////////////////////////////////////////////////
//polar plots
/////////////////////////////////////////////////////////////////////////////////////////////// 
    //polar plotting

    x = 0:2*%pi/(length(zArraySymmMAG(:,1))-1):2*%pi;

    scf(15)
    //polarplot(x,odmkOscMAGwin(:,5),style=[23])
    //polarplot(x,zArrayInMAG(2:1+nfft/2,1),style=[23])
    polarplot(x,zArraySymmMAG(:,1),style=[23])
    legend('y = odmkOscPHASEwin(:,3)',4)


    scf(16)
    //polarplot(x,odmkOscPHASEwin(:,5),style=[23])
    //polarplot(x,zArrayInPHASE(2:1+nfft/2,1),style=[23])
    polarplot(x,zArraySymmPHASE(:,1),style=[23])
    legend('y = odmkOscPHASEwin(:,3)',4)

    
    //polar plot of composite frames 1-16 (arbitrary)
    scf(17)
        for i=1:15
            //polarplot(x,zArrayInMAG(2:1+nfft/2,i),style=[i]);
            polarplot(x,zArraySymmMAG(:,i),style=[i]);
        end
    legend('y = odmkOscPHASEwin(:,3)',4)


    scf(18)
        for i=1:15
            //polarplot(x,zArrayInPHASE(2:1+nfft/2,i),style=[i]);
            polarplot(x,zArraySymmPHASE(:,i),style=[i]); 
        end
    legend('y = odmkOscPHASEwin(:,3)',4)



    //spherical coordinates - r = radius, theta = incination/polar angle(z,r plane), phi = azimuth(x,y plane)
    //spheric -> cartesian:
    //    x = r*sin(theta)*cos(phi)
    //    y = r*sin(theta)*sin(phi)
    //    z = r*cos(theta)
    //
    //cartesian -> spheric
    //    r = sqrt((x^2)+(y^2)+(z^2))
    //    theta = arccos(z/r)
    //    phi = arctan(y/x)



    //for j=1:n
    //zRo(j) = abs(zn(j));
    //zTheta(j) = atan(imag(zn(j)),real(zn(j)));
 
 
 
 
    

        //create 3D coordinate vectors:
    //x = frequency axis, 0 - Nfft/2
    //y = constant offset spacing for each slice FFT magnitude array 
    //create array of color id values (what is scilab max?)    
//scf(7)    
    //param3d1(x,y,list(z,colors),[theta,alpha,leg,flag,ebox])
    //colors = 16*ones(nSlice,1);
    //for j=1:nSlice
        //colors(j) = j;
        //xx(:,j) = x';
        //yy(:,j) = (4*j).*ones(nSlice,1);
        //param3d(xx,yy,zArrayInx,[70,30,"X@Y@Z",[2,3]])        
    //end
    //param3d1(xx,yy,list(zArrayIn,colors),-46,86,"X@Y@Z",[2,3])   
    //a=get("current_axes");//get the handle of the newly created axes
    //a.rotation_angles=[65,75]; 
    //a.data_bounds=[-1,-1,-1;1,1,2]; //boundaries given by data_bounds
    //a.thickness = 2;
    
    
endfunction


////////////////////////////////////////////////////////////////             
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : waterfall plot function
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : Wave file functions : xwavread, xwavewrite
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

////////////////////////////////////////////////////////////////
//function: xwavread
//Read Microsoft .wav sound file

function [y,Fs,bits]=xwavread(wavfile,ext)
    y=[];Fs=[];bits=[];
    //Read Microsoft .wav sound file.
    //y=wavread(wavfile) reads a .wav file specified by the string wavfile,
    //returning the sampled data in y. The .wav extension is appended
    //if no extension is given.  
    //Amplitude values are in the range [-1,+1].
    //[y,fs,bits]=wavread(wavfile) returns the sample rate (fs) in Hertz
    //and the number of bits per sample (bits) used to encode the
    //data in the file.
    //[...]=wavread(wavfile,n) returns only the first n samples from each
    //      channel in the file.
    //[...]=wavread(wavfile,[n1 n2]) returns only samples n1 through n2 from
    //     each channel in the file.
    //siz=wavread(wavfile,'size') returns the size of the audio data contained
    //    in the file in place of the actual audio data, returning the
    //     vector siz=[samples channels].

    //
    // xwavread modifications from wavread:
    // supports 24 and 32bits.
    // License GNU GPL
  

    // Append .wav extension if necessary
    if strindex(wavfile,'.')==[] then    //searches string input (wav file name) for '.'
        wavfile = wavfile+'.wav';
    end
    
    // Open the file (little-endian=1, probably don't need)
    // fid = positive integer for identification
    // 'rb' = open a binary file for reading
    [fid,err] = mopen(wavfile,'rb',1);    
    if err<0 then  
        error('Cannot open '+wavfile),
    end

    // Handle ext optional argument
    if argn(2)<2 then
        ext = [];
    elseif type(ext)==10 then    //10 indicates a string (scilab type code)
        ext=convstr(ext)
        //if ext<>'size' then    //<> is different than
        if ext~='size' then
            error('second argument must be ""size"", an integer or a vector of 2 integers")
        end
    //elseif type(ext)==1 then    //1 indicates a double (scilab type code)
    elseif type(ext)==8 then    //8 indicates an integer (scilab type code)
        exts = size(ext,'*');    //returns the product of the dimensions
        if exts>2 then
            error('Index range must be specified as a scalar or 2-element vector.');
        elseif exts==1 then
            if ext==0 then
                ext = 'size';  // synonym for size (confusing, is this necessary?)
            else
                ext = [1,ext];
            end     
        end
    else
        error('second argument must be ""size"", an integer or a vector of 2 integers")
    end
  


    Data=[];
    //The mget function reads data in the input specified by the stream parameter fd (here fid)
    //and returns a vector of floating point data
    ID=stripblanks(ascii(mget(4,'c',fid)));    //4 bytes chunk ID
    //check to make sure the ID matches with riff, (microsoft wave format ID)
    mprintf('chunk ID: %s\n',ID);
    if (convstr(ID)~='riff') then
        error('.wav file does not contain the RIFF identifier')
    end
    Size=mget(1,'ui',fid);    //?4 bytes Chunk Size in bytes
    mprintf('chunk size: %i\n',Size);
    rifftype=mget(4,'c',fid);
    dtype = convstr(ascii(rifftype)');
    mprintf('dtype: %s\n\n',convstr(ID));
    if dtype~='wave' then
        error('.wav file does not contain the wave identifier')
    end
  
    // Find optional chunks
    found_fmt = 0; found_data = 0;
    while ~found_data then
        [ID,Size]=find_cktype(fid);
        select ID
            case 'fact' then
                total_bytes=Size;
                orig_pos = mtell(fid);
                nbytes = 4;
                // # of required bytes in <fact-ck> header
                if total_bytes < nbytes then
                    error('Error reading .wav file');
                end
                factdata=mget(1,'ui',fid); // Samples per second
                rbytes = total_bytes-(mtell(fid)-orig_pos);    //mtell shows current position in fid
                if rbytes then
                    if mseek(rbytes,fid,'cur')==-1 then
                        error('Error reading <fact-ck> chunk.');
                    end
                end
            case 'fmt' then
                found_fmt = 1;
                [wFormatTag,nChannels,nSamplesPerSec,nAvgBytesPerSec,...
                nBlockAlign,nBitsPerSample,cbSize]=read_wavefmt(fid,Size);
                mprintf('Loaded WAV File: %s\n',wavfile);
                mprintf(' wFormatTag = %i\n',wFormatTag);
                mprintf(' Sample Rate = %f (Hz)\n',nSamplesPerSec);
                mprintf(' Ave Bytes per sec = %i\n',nAvgBytesPerSec);
                mprintf(' Bits Per Sample = %i\n',nBitsPerSample);
                mprintf(' Number of Channels = %i\n',nChannels);
                mprintf(' Block Alignment = %i\n\n',nBlockAlign);
            case 'data' then
                found_data = 1;
                if ~found_fmt then
                    error('Invalid .wav file: found data before format information.');
                end
                if (ext=='size')|(~(ext==[]))&and(ext==0) then    //tricky, enter 0 <=> 'size'
                    // Caller just wants data size:
                    samples=read_wavedat(fid, Size ,wFormatTag, nChannels, nBitsPerSample,-1)
                    mclose(fid);
                    y = [samples,nChannels];
                else
                    //y is [#data samples][#channels] vector of two integers
                    y=read_wavedat(fid, Size, wFormatTag, nChannels, nBitsPerSample, ext)
                    mclose(fid);
                    szy=size(y);
                    total_samples=szy(1)*szy(2);
                    mprintf('Total samples = %i\n\n',total_samples);
                end
                //mprintf('Size (# of bytes) = %i\n',Size)
                //if mseek(Size,fid,'cur')==-1 then
//                if mseek(Size,fid,'cur')<0 then
//                    error('Incorrect chunk size information in RIFF file.');
//                end
            end
        end
    Fs = nSamplesPerSec;
    if wFormatTag==1 then
        bits = nBitsPerSample;
    elseif wFormatTag==3 then
        bits = nBitsPerSample;
    else
        bits = [];
    end
endfunction

////////////////////////////////////////////////////////////////
//function: find_cktype
//find Microsoft .wav sound file ID type

function [ID,Size]=find_cktype(fid)
    ID=stripblanks(ascii(mget(4,'c',fid)));
    mprintf('check for ID type: %s\n',ID)
    Size=mget(1,'ui',fid);
    mprintf('check for Size (# of bytes) = %ui\n\n',Size)
endfunction

////////////////////////////////////////////////////////////////
//function: read_wavefmt
//Read Microsoft .wav sound file format

function [wFormatTag,nChannels,nSamplesPerSec,nAvgBytesPerSec,nBlockAlign,nBitsPerSample,cbSize]=read_wavefmt(fid,total_bytes)
    orig_pos = mtell(fid); 
    nbytes = 14; // # of required bytes in  header

    if total_bytes<nbytes then
        error('Error reading .wav file');
    end
  
    // Read wav data:
    wFormatTag=mget(1,'us',fid);  // Data encoding format : actually 3 is used??
    nChannels=mget(1,'us',fid); // Number of channels
    nSamplesPerSec=mget(1,'ui',fid); // Samples per second
    nAvgBytesPerSec=mget(1,'ui',fid); // Avg transfer rate
    nBlockAlign=mget(1,'us',fid); // Block alignment
  
    //if (wFormatTag==1) | (wFormatTag==3) then
    if wFormatTag==1 then
        [cbSize,nBitsPerSample]=read_fmt_pcm(fid,total_bytes)
        bits = nBitsPerSample;
//    if wFormatTag==1 then
//        [cbSize,nBitsPerSample]=read_fmt_pcm(fid,total_bytes)
//        bits = nBitsPerSample;    
//    elseif wFormatTag==3 then
//        [cbSize,nBitsPerSample]=read_fmt_pcm(fid,total_bytes)
//        bits = nBitsPerSample;         
    else
        mprintf(' wFormatTag = %us\n',wFormatTag);
        error('Invalid wav format');
    end

    rbytes = total_bytes-(mtell(fid)-orig_pos);
    if rbytes then
        if mseek(rbytes,fid,'cur')==-1 then
            error('Error reading wav file');
        end
    end
  
endfunction

////////////////////////////////////////////////////////////////
//function: read_fmt_pcm
//find pcm format

function [cbSize,nBitsPerSample]=read_fmt_pcm(fid,total_bytes)
    nbytes = 14; cbSize=[]; nBitsPerSample=[]; 
    // # of bytes already read 
    if total_bytes < nbytes+2 then
        error('Error reading wav file');
    end
    nBitsPerSample = mget(1, 'us', fid);
    nbytes = nbytes+2;
    if total_bytes > nbytes then
        if total_bytes>=nbytes+2 then
            cbSize = mget(1,'us',fid);
            nbytes = nbytes+2;
        end
//    if total_bytes>nbytes then
//        if mseek(total_bytes-nbytes,fid,'cur')==-1 then
//            error('Error reading .wav file');
//        end
//    end
    end
endfunction


////////////////////////////////////////////////////////////////
//function: read_wavedat
//read .wav file data

function Data=read_wavedat(fid, Size, wFormatTag, nChannels, nBitsPerSample, ext)
    fmt_msg = [];
    select wFormatTag
        case 1 then
            // PCM Format: limited to 16bit 2 channels
            Data=read_dat_pcm(fid, Size, nChannels, nBitsPerSample, ext)
        case 2 then
            fmt_msg = 'Microsoft ADPCM';
        case 3 then
            // #define WAVE_FORMAT_IEEE_FLOAT 3, nBitsPerSample 32i? (float)
            Data=read_dat_pcm(fid, Size, nChannels, nBitsPerSample, ext)
        case 6 then
            fmt_msg = 'CCITT a-law';
        case 7 then
            fmt_msg = 'CCITT mu-law';
        case 17 then
            fmt_msg = 'IMA ADPCM';
        case 34 then
            fmt_msg = 'DSP Group TrueSpeech TM';
        case 49 then
            fmt_msg = 'GSM 6.10';
        case 50 then
            fmt_msg = 'MSN Audio';
        case 257 then
            fmt_msg = 'IBM Mu-law';
        case 258 then
            fmt_msg = 'IBM A-law';
        case 259 then
            fmt_msg = 'IBM AVC Adaptive Differential';
        else
            fmt_msg = 'Format #'+string(wFormatTag);
        end
    if ~(fmt_msg==[]) then
        error('Data compression format '+fmt_msg+' is not supported');
    end
endfunction


////////////////////////////////////////////////////////////////
//function: read_dat_pcm
//read .wav file pcm data

function Data=read_dat_pcm(fid, total_bytes, nChannels, nBitsPerSample, ext)
// Determine # bytes/sample - format requires rounding
//  to next integer number of bytes:
    BytesPerSample = ceil(nBitsPerSample/8);
    select BytesPerSample
        case 1 then  // unsigned 8-bit
            dtype = 'uc';
        case 2 then // signed 16-bit
            dtype = 's';
        case 3 then // signed 24-bit    WAVE_FORMAT_EXTENSIBLE
            dtype = 'c';
        case 4 then // signed 32-bit (long)    WAVE_FORMAT_EXTENSIBLE
            dtype = 'l';
        else
            error('Cannot read .wav file  with this number of bits per sample');
        end//  select BytesPerSample
        
    // # bytes in this chunk
    total_samples = total_bytes/BytesPerSample;
    SamplesPerChannel = total_samples/nChannels;
    
    if (~(ext==[]))&(ext==-1) then
        // Just return the samples per channel, and seek past data:
        Data = SamplesPerChannel;
//        if mseek(total_bytes,fid,'cur')==-1 then
//            error( 'Error reading .wav file');
//        end
        return
    elseif ext==[] then
        ext = [1,SamplesPerChannel];
    elseif prod(size(ext))~=2 then
       error('Sample limit vector must have 2 entries');
       return
    elseif ext(1)<1|ext(2)>SamplesPerChannel then
       error('Sample limits out of range.');
    elseif ext(1)>ext(2) then
       error('Invalid sample limits (use ascending order)');
    elseif ext(1)>1 then
        // Skip if specified:
        if mseek(BytesPerSample*(ext(1)-1)*nChannels,fid,'cur')==-1 then
        error('Error reading .wav file');
        return
        end
    end
  
    // Read data:
    nSPCext = ext(2)-ext(1)+1;
    // # samples per channel in extraction range
    extSamples = nChannels*nSPCext;
    //24-bits files need special treatment
    if (BytesPerSample)==3 then
        Data_tmp=[];
        Data_tmp=(mget(3*nChannels*nSPCext,dtype,fid));
        oct1=uint8(Data_tmp(1:3:$-1));
        oct2=uint8(Data_tmp(2:3:$));
        oct3=Data_tmp(3:3:$);
        Data_tmp2=(double(oct1)*(2^0)+double(oct2)*(2^8)+double(oct3)*(2^16));
        //y=matrix(v,[sizes]) = v=a vector or a matrix; sizes= vector of integers
        Data=matrix(Data_tmp2,[nChannels,nSPCext]);
    else
        Data=matrix(mget(nChannels*nSPCext,dtype,fid),[nChannels,nSPCext])
    end;//if 24bits

    // Skip trailing samples:
    //  if mseek(BytesPerSample*(SamplesPerChannel-ext(2))*nChannels,fid,'cur')==-1 then
    //    error('Error reading .wav file');
    //    return
    //  end
  
    // Determine if a pad-byte is appended and skip if present:
    junk = Size;
    if junk-fix(junk./2).*2 then
        mseek(1,fid,'cur');
    end
    Data=Data';
    // Normalize data range in [-1 1] (min will hit -1)
    select BytesPerSample
        case 1 then
            Data=(Data-128)/128;
        case 2 then 
            Data=Data/32768;
        case 3 then 
            Data=Data/(2^23);
        case 4 then
            Data=Data/(2^31);
        end;//normalisation in range [-1 +1]
        
endfunction



////////////////////////////////////////////////////////////////
//function: xwavwrite
//write .wav file

function []=xwavwrite(y,Fs,nbits,wavefile)
    [nargout,nargin] = argn(0)
  //WAVWRITE Write Microsoft WAVE ("""".wav"""") sound file.
  //   WAVWRITE(Y,WAVEFILE) writes a WAVE file specified by the
  //   string WAVEFILE.  The data should be arranged with one channel
  //   per column.  Amplitude values outside the range [-1,+1] are
  //   clipped prior to writing.
  // 
  //   WAVWRITE(Y,FS,WAVEFILE) specifies the sample rate FS of the data
  //   in Hertz.
  // 
  //   WAVWRITE(Y,FS,NBITS,WAVEFILE) forces an NBITS-bit file format to
  //   be written, where NBITS<=16.
  // 
  //   Supports multi-channel 8-24-bit WAVE data.


  
  // Get user default preferences:
    Fs_pref = 44100.000;
    nbits_pref = 24;

    // Parse inputs:
    if nargin<2|nargin>4 then
        error('Incorrect number of input arguments.');
    elseif nargin<3 then
         wavefile = Fs;
         Fs = Fs_pref;
         nbits = nbits_pref;
    elseif nargin<4 then
        wavefile = nbits;
        nbits = nbits_pref;
    end
  
    // Open file for output:
    if ~(type(wavefile)==10) then
        error('wavefile must be a string.');
    end
    if mtlb_findstr(wavefile,'.')==[] then
        wavefile = wavefile+'.wav';
    end
    
    [fid,%v] = mopen(wavefile,'wb',1)
    if %v<0 then fid = -1;end
    // Little-endian
    if fid==(-1) then
        error('Can''t open WAVE file for output.');
    end
  
    // If input is a vector, force it to be a column:

    if size(y,2) > 2 then
        error('Data array must have 1- or 2-dimensions, only.');
    end
    [samples,channels] = size(y);
     if samples==1 then
        y = y(:);
        [samples,channels] = size(y);
    end
  
    // Clip data to normalized range [-1,+1]:
    i = matrix(find(abs(y)>1),1,-1);
    if ~(i==[]) then
        y(i) = sign(y(i)) 
        warning('Data clipped during write to file.');
    end
  
    // # bytes per sample to write
    bytes_per_sample = ceil(nbits/8);
    total_samples = samples*channels;
    total_bytes = total_samples*bytes_per_sample;
  
    // Determine number of bytes in RIFF chunk
    // (not including pad bytes, if needed):
    // ----------------------------------
    //  'RIFF'           4 bytes
    //  size             4 bytes (ulong)
    //  'WAVE'           4 bytes
    //  'fmt '           4 bytes
    //  size             4 bytes (ulong)
    // <wave-format>     14 bytes
    // <format_specific> 2 bytes (PCM)
    //  'data'           4 bytes
    //  size             4 bytes (ulong)
    // <wave-data>       N bytes
    // ----------------------------------
    riff_cksize = 36+total_bytes;
    // Don't include 'RIFF' or its size field
    fmt_cksize = 16;
    // Don't include 'fmt ' or its size field
    data_cksize = total_bytes;
    // Don't include 'data' or its size field
  
    // Determine pad bytes:
    data_pad = data_cksize-fix(data_cksize./2).*2;
    riff_cksize = riff_cksize+data_pad;
    // + fmt_pad, always 0

    ck= tlist(['ck','fid','Size','ID']) 
    // Write RIFF chunk:
    ck('fid')=fid
    ck('Size')=riff_cksize
    ck('ID')='RIFF';
    write_ckinfo(ck);

    // Write WAVE:
    ck('ID')='WAVE';
    write_ckinfo(ck,1);
  
    // Write <fmt-ck>:
    ck('ID')='fmt ';
    ck('Size')=fmt_cksize
    write_ckinfo(ck);
    // Write <wave-format>:
    fmt=tlist(['fmt','wFormatTag','nChannels','nSamplesPerSec','nAvgBytesPerSec','nBlockAlign','nBitsPerSample'])
    fmt('wFormatTag')=1
    //fmt('wFormatTag')=3    //for IEEE floating point
    // Data encoding format = PCM
    fmt('nChannels')=channels
    // Number of channels
    fmt('nSamplesPerSec')=Fs
    // Samples per second
    fmt('nAvgBytesPerSec')=channels*bytes_per_sample*Fs

    // Avg transfer rate
    fmt('nBlockAlign')=channels*bytes_per_sample
    // Block alignment
    fmt('nBitsPerSample')=nbits
    // standard <PCM-format-specific> info
    
    mprintf(' \nCreating WAV File: %s\n',wavefile);
    mprintf(' wavefile format (fmt)');
    mprintf(' wFormatTag = %i\n',fmt('wFormatTag'));
    mprintf(' Sample Rate = %f (Hz)\n',fmt('nSamplesPerSec'));
    mprintf(' Ave Bytes per sec = %i\n',fmt('nAvgBytesPerSec'));
    mprintf(' Bits Per Sample = %i\n',fmt('nBitsPerSample'));
    mprintf(' Number of Channels = %i\n',fmt('nChannels'));
    mprintf(' Block Alignment = %i\n\n',fmt('nBlockAlign'));
    
    status = write_wavefmt(fid,fmt);
    //mprintf('<xwavwrite status>) <write_wavefmt> success!  = %i\n',status);
    // Write <data-ck>:
    ck('ID') = 'data';
    ck('Size')=data_cksize
    write_ckinfo(ck);
  
    // Write <wave-data>, and its pad byte if needed:
    status = write_wavedat(fid,fmt,y);
    mclose(fid);
    // Close file
    // end of wavwrite()
endfunction

function write_ckinfo(ck,sflg)
    [nargout,nargin] = argn(0)
    // WRITE_CKINFO: Writes next RIFF chunk, but not the chunk data.
    //   If optional sflg is set to nonzero, write SUBchunk info instead.
    //   Expects an open FID pointing to first byte of chunk header,
    //   and a chunk structure.
    //   ck.fid, ck.ID, ck.Size, ck.Data
    //if length(ck('ID'))<>4 then    //<> op: "is different than"
    if length(ck('ID'))~=4 then
        error('write_ckinfo: ck has a wrong length')
    end
    mput(ascii(ck('ID')),'c',ck('fid'))
    // Error condition
    if nargin==1 then
        // Write chunk size (skip if subchunk):
        mput(ck('Size'),'ui',ck('fid'))
    end
endfunction

function [status]=write_wavedat(fid,fmt,data)
    status=[];
    // WRITE_WAVEDAT: Write WAVE data chunk
    //   Assumes fid points to the wave-data chunk
    //   Requires <wave-format> structure to be passed.
  
    status = 0;
  
    //if fmt('wFormatTag')==1s then
    if fmt('wFormatTag')==1 then
        // PCM Format:
        // Determine # bytes/sample - format requires rounding
        //  to next integer number of bytes:
        BytesPerSample = ceil(fmt('nBitsPerSample')/8);
        mprintf(' Bytes Per Sample = %i\n',BytesPerSample);
        select BytesPerSample
            case 1 then
                dtype = 'uc'; // unsigned 8-bit
                // Scale data according to bits/samples: [-1,+1] -> [0,255]
                data = round((data+1)*255/2);
            case 2 then
                dtype = 's';
                // signed 16-bit
                // Scale data according to bits/samples: [-1,+1] -> [-32768,+32767]
                data = round((data+1)*65535/2)-32768;
            case 3 then
                dtype='c'    //c = char, is this correct?
                // signed 24-bit  WAVE_FORMAT_EXTENSIBLE
                // Scale data according to bits/samples: [-1,+1] -> [-8 388 608,+8 388 607]
                data = round((data+1)*(2^24-1)/2)-(2^23);
            case 4 then 
                dtype='l'
                // signed 32-bit  WAVE_FORMAT_EXTENSIBLE
                // Scale data according to bits/samples: [-1,+1] -> [-2 147 483 648,+2 147 483 647]
                data = round((data+1)*(2^32-1)/2)-(2^31);
            else
                error('Cannot write WAVE files with this bits/sample.');
            end
        //mprintf('<write_wavedat dtype>) data type = %s\n',dtype);
        // Write data, one row at a time (one sample from each channel):
        [samples,channels] = size(data);
        total_samples = samples*channels;
        mprintf(' Total samples = %i\n\n',total_samples);

        //24-bits needs special treatment
        if BytesPerSample==3 then
            oct3=(floor((data)/(2^16)));//msb
            oct2=(floor((data-(oct3*2^16))/(2^8)));
            oct1=(floor(data-(oct3*2^16)-(oct2*2^8)));//lsb
            data_line=zeros(3*total_samples,1);
            data_line(1:6:$)=(oct1(:,1));  // (k:n:$)=>start at k, skip n steps of array addresses)
            data_line(2:6:$)=(oct2(:,1));
            data_line(3:6:$)=(oct3(:,1));
            data_line(4:6:$)=(oct1(:,2));
            data_line(5:6:$)=(oct2(:,2));
            data_line(6:6:$)=(oct3(:,2));
            data_line=data_line';
            
            szdl=size(data_line(1,:));
            //mprintf('<write_wavedat data_line) size n = %i, m = %i\n',szdl(1),szdl(2));
            //if mput(data_line,"d",fid)~=0 then     
            mput(data_line,"c",fid);
            //mprintf('<write_wavedat status>) data written to file');
                //status = -1;
                //return
            //end;
        //if format is different from 24-bits
        else
            mput(matrix(data',total_samples,1),dtype,fid);
            //status = -1;
            //return
        end
    
        // Error condition?
        // Determine if a pad-byte is appended to data chunk:
        %v2_1$1 = total_samples*BytesPerSample
        if %v2_1$1-fix(%v2_1$1./2).*2 then
            mtlb_fwrite(fid,0,'uchar');
        end   
    else
        // Unknown wave-format for data.
        error('Unknown data format.');
    end
    return
endfunction


function [status]=write_wavefmt(fid,fmt)
  status = 0;
  // WRITE_WAVEFMT: Write WAVE format chunk.
  //   Assumes fid points to the wave-format subchunk.
  //   Requires chunk structure to be passed, indicating
  //   the length of the chunk.
  
  // Create <wave-format> data:

  mput(fmt('wFormatTag'),'us',fid);
  mput(fmt('nChannels'),'us',fid);
  mput(fmt('nSamplesPerSec'),'ui',fid);
  mput(fmt('nAvgBytesPerSec'),'ui',fid);
  mput(fmt('nBlockAlign'),'us',fid);
  
  // Write format-specific info:
  if (fmt('wFormatTag')==1) then
    // Write standard <PCM-format-specific> info:
    mput(fmt('nBitsPerSample'),'us',fid)
  else
     error('Unknown data format.');
  end
endfunction


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : Wave file functions : xwavread, xwavewrite
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


//read .wav audio file from directory

//wavfile_in='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavSource/odiusEMattack002.wav';
//wavfile_in='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavSource/death_by_black_hole_vox.wav';
//wavfile_in='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavSource/YT-KONGA2x02.wav';
wavfile_in='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavSource/601_90ALTOSAX_C_chop.wav';
//wavfile_in='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavSource/The_Amen_Break_odmk.wav';


[audio_in,fs,bits]=xwavread(wavfile_in);

szA=length(audio_in(:,1));
nchanA=length(audio_in(1,:));
mprintf(' <xwavread out> # of audio samples = %i x # of channels = %i\n',szA,nchanA);
mprintf(' <xwavread fs output> Sample Rate = %f (Hz)\n',fs);
mprintf(' <xwavread bits output>  bits per sample = %us\n\n',bits)

audioDataL=audio_in(:,1);
audioDataR=audio_in(:,2);


Dlength = length(audioDataL);

wavDataL = audio_in((1:Dlength),1);
wavDataR = audio_in((1:Dlength),2);


pvNFFT = 1024;
aWin = pvNFFT;
sWin = 0;       //frames are hann-windowed at nfft points when nWin is integer/=0, or rectangular if W=0
Ra_hop = .5;    //analysis hop factor, %of time-slice
Rs_hop = .5;    //synthesis hop factor, %of time-slice
zp = 0;            //zero-padding (nfft -> 2*nfft, increases freq resolution)
pkDet = 1;        //turns on peak detection
maxPeaks = 15;    //max number of peaks calculated by peak detector
specPlots = 0;   //3D waterfall plots of spectral data
noSpecPlots = 0;

//1024 samples is about 60 ms at 16kHz, a good window
wavVocodedL=phaseVoc(wavDataL,Ra_hop,Rs_hop,pvNFFT,aWin,sWin,zp,fs,pkDet,maxPeaks,specPlots);
wavVocodedR=phaseVoc(wavDataR,Ra_hop,Rs_hop,pvNFFT,aWin,sWin,zp,fs,pkDet,maxPeaks,noSpecPlots);
// Compare original and resynthesis


wavVocodedLscale = .99*(wavVocodedL./(max(wavVocodedL)));
wavVocodedRscale = .99*(wavVocodedR./(max(wavVocodedR)));


//wavVocoded(:,1) = wavVocodedL;
//wavVocoded(:,2) = wavVocodedR;

wavVocoded(:,1) = wavVocodedLscale;
wavVocoded(:,2) = wavVocodedRscale;
//fs=44100;
//bits=24;

wavVocoded_out='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavResult/wavVocoded_out';
xwavwrite(wavVocoded,fs,bits,wavVocoded_out);
mprintf(' <xwavwrite out> wrote .wav file %s\n\n',wavVocoded_out)




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



    //darkmagenta
    //darkcyan
    //steelblue
    //darkred
    //darkblue
    //darkred
    //darkslategray1-4
    //turquoise1-4
    //cadetblue1-4
    //chartreuse1
    //lightgoldenrod1
    //gold1
    //darkgoldenrod1
    //coral1
    //tomato1
    //orange
    //orangered1
    //deeppink1
    //hotpink1
    //pink1
    //palevioletred1
    //maroon1
    //violetred1
    //magenta1
    //orchid1
    //plum1
    //mediumorchid1
    //purple1
    //mediumpurple1
    //darkgrey
    

///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : plotting
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////