// /////////////////////////////////////////////////////////////////////////////////////////////////
// ################################################################################################
// begin : header
// ################################################################################################
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// ___::((odmk_osc5_tb))::___
//
// ___<<name=> "odmk_osc5.sce" >>___
//
// ___::((JIROBATA Programming Industries))::___
// ___::((ODMK:odorousbeast:BarutanBreaks:djoto:2014:2015:2016))::___
// ___created by eschei
//
//
// include components:
//
// ___::((osc5))::___
// ___::((odmk_multiOsc))::___
// ___::((table generator))::___
// ___::((Interpolated WaveTable))::___
// ___::((24 but .wav file write))::___
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This script implements wavetable complex audio oscillator
//
// odmk_osc5: single channel complex wavetable oscillator
//
// odmk_osc5Stereo: Stereo channel complex wavetable oscillator
//
// odmk_multiOsc: multi-channel complex wavetable oscillator
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// /////////////////////////////////////////////////////////////////////////////////////////////////
// ################################################################################################
//end : header
// ################################################################################################
// /////////////////////////////////////////////////////////////////////////////////////////////////


xdel(winsid()) //-> closes all open graphs (= matlab.close('all'))
clear
clc

//include .wav read/write functions
exec('/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/wav_rdwr_scilab/odmkWav_rw.sce');

///////////////////////////////////////////////////////////////////////////////////////////////
//begin : stack adjustment
///////////////////////////////////////////////////////////////////////////////////////////////

////try for ~10MBytes
////stacksize input is doubles
////
//max_bytes = 10*10e6;
//max_bits=16;
////bits=24;
//bytespersample=ceil(max_bits/8);
//max_data_bytes=max_bytes-(12+24);    //total size - header,format,etc. approx
//max_stack=max_data_bytes/8;
//
////if size > max then
////    error('wav file too large');
////else
//stacksize(max_stack)
////    stacksize('max')
////    stacksize('min')
////    sz=stacksize()
////end

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
//begin : DDS functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////////////////////////////////////////////////////////////////
//function: tablegen
//look-up table generator


function [tbout] = tablegen(shape,depth)
    
    //create look-up table entries for different waveforms
    //define an n deep table
    t1=(0:1/(depth):1);
    t1=t1(1:depth)            //crops off the last point for proper cyclic table
    t2=(1:-1/(depth-1):0)^3;
    t3=[0:1/(depth/2-1):1,1:-1/(depth/2-1):0]^3;  
    table1=zeros(length(t1),1);
    select shape
        case 1 then  // store 1 cycle of sin
            for q=1:length(t1)
                table1(q)= sin(2*%pi*(t1(q))); 
            end
        case 2 then  // store 1 cycle of tri
                table1=[(0:.999/(depth/4-1):.999)';(.999:-.999/(depth/4-1):0)';(0:-.999/(depth/4-1):-.999)';(-.999:.999/(depth/4-1):0)']; 
        case 3 then  // store 1 cycle of saw-up
                table1=[(-0.999:.999/(depth/2-1):.999)'];
        case 4 then  // store 1 cycle of saw-down
                table1=[(0.999:-.999/(depth/2-1):-.999)'];
        case 5 then  // store 1 cycle of chebychev  
            for r=1:length(tl)
                table1(r)=cos(13*acos(tl(r)));
            end
        case 6 then  // store 1 cycle of pulse1  
            for s=1:length(t2)
                table1(s)=sin(5*%pi*(t2(s))); 
            end 
        case 7 then  // store 1 cycle of pulse2  
            for s=1:length(t2)
                table1(s)=sin(9*%pi*(t2(s))); 
            end 
        case 8 then  // store 1 cycle of pulse3 
            for s=1:length(t2)
                table1(s)=sin(23*%pi*(t2(s))); 
            end
        case 9 then  // store 1 cycle of pulse4  
            for t=1:length(t3)
                table1(t)=cos(5*%pi*(t3(t))); 
            end             
        else
            for t=1:length(t1)
                table1(t)=sin(2*%pi*(tl(t))); 
            end
    end
    tbout = table1;
    
endfunction


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//function: odmk_osc5 - begin
////////////////////////////////////////////////////////////////

//wavetable oscillator function

function [odmkOsc, odmkOsc90, odmkSqrPulse] = odmk_osc5(numSamples, Fs, tableDepth, shape, foArray, phArray)
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // odmk_osc5: single channel complex wavetable oscillator
    //
    // The output frequency can be fixed or variable
    // When foArray is a single scalar value, the output freq. is fixed
    // When foArray is an array of length=numSamples, the output freq. varies each step
    //
    // The output phase can be fixed or variable
    // When phaseArray is a single scalar value, the output phase is fixed
    // When phaseArray is an array of length=numSamples, the output phase varies each step 
    //
    // Real time controls:
    // frequency ; phase
    //
    // data format: input foMatrix, phaseMatrix => (numSamples,1)
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tb = zeros(tableDepth,1);
    tb = tablegen(shape,tableDepth);
    
    //check inputs
    if (length(foArray) == 1) then
        for n=1:numSamples
            foArr(n) = foArray;
        end
    elseif (length(foArray) == numSamples) then
        foArr = foArray;
    else
        error('length of foArray must be 1 or match numSamples')
    end
    
    if (length(phArray) == 1) then
        for n=1:numSamples
            phaseArr(n) = phArray;
        end
    elseif (length(phArray) == numSamples) then
        phaseArr = phArray;
    else
        error('length of phArray must be 1 or match numSamples')
    end
    
    
    //***initialize***
    
    skipInc = zeros(numSamples,1);
    phOffset = zeros(numSamples,1);
    
    //init osc output vectors
    odmkOsc = zeros(numSamples,1);
    odmkOsc90 = zeros(numSamples,1);
    odmkSqrPulse = zeros(numSamples,1);
    
    accWidth = 48;
    qntWidth = ceil(log2(tableDepth));
    lsbWidth = accWidth-qntWidth;    //expect 36
    //mprintf('lsb_width %i\n\n',lsb_width)
    
    foArrScale = (2^lsbWidth)*foArr;
    
    
    //phase accumulator increment value (skipping value)
    //skipInc = floor(fo*(2^N/fs)+0.5);
    //skipInc = floor(foscale*(2^qntWidth)/fs+0.5);
    //skipInc = floor((foscale*(2^qntWidth))/fs);
    skipInc = ((2^qntWidth).*foArrScale)./fs; 

    //scale phase offset
    //converts radians into a scaled phase offset to be added to the output of the acc
    //dependent on the table depth
    phOffset = round(modulo(((phArray.*tableDepth)./(2*%pi)),tableDepth));
    

    //used to add a 90 deg offset for complex sinusoid generation
    offset90 = tableDepth/4; 

    //**main loop**
    //___::((Interpolated WaveTable))::___
    //generates main osc waveform out, 90deg shifted out (sin/cos), square pulse out

    //initialize variables so that the OSC starts at the beginning of the table
    accAddr = phOffset(1)*2^lsbWidth;    //init to 'zero plus phase shift'
    accAddr90 = (offset90+phOffset(1))*2^lsbWidth;    //init to 'zero plus offset90 plus phase shift'
    qntAddr = 1+phOffset(1);
    qntAddr90 = qntAddr+offset90;
    yLow = tb(qntAddr);
    yLow90 = tb(qntAddr90);

    for i=1:numSamples;    //osc main loop

        accAddrP1 = accAddr+(2^lsbWidth);
        accAddr90P1 = accAddr90+(2^lsbWidth);

        // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        // generate oscillator output
        // ///////////////////////////////////////////////////        
        if (qntAddr==tableDepth) then    //avoid addr overflow for linear interpolation addressing - should investigate
            yHigh = tb(qntAddr);
            yHigh90 = tb(qntAddr90+1);
            odmkOsc(i) = tb(tableDepth);
            //linear interpolation
            odmkOsc90(i) = yLow90 + (yHigh90-yLow90)*((accAddr90P1-(qntAddr90*2^lsbWidth))*2^-lsbWidth);
            //temp
            yHigh_tap(i) = yHigh;
            yHigh90_tap(i) = yHigh90;
        elseif (qntAddr90==tableDepth) then
            yHigh = tb(qntAddr+1);
            yHigh90 = tb(qntAddr90);
            //linear interpolation
            odmkOsc(i) = yLow + (yHigh-yLow)*((accAddrP1-(qntAddr*2^lsbWidth))*2^-lsbWidth);
            odmkOsc90(i) = tb(tableDepth);
        else
            yHigh = tb(qntAddr+1);
            yHigh90 = tb(qntAddr90+1);
            odmkOsc(i) = yLow + (yHigh-yLow)*((accAddrP1-(qntAddr*2^lsbWidth))*2^-lsbWidth);
            odmkOsc90(i) = yLow90 + (yHigh90-yLow90)*((accAddr90P1-(qntAddr90*2^lsbWidth))*2^-lsbWidth);
        end 

        //generate square pulse output
        if odmkOsc(i) >= 0 then
            odmkSqrPulse(i) = 1;
        else
            odmkSqrPulse(i) = 0;
        end

        //phase accumulator
        accAddr = modulo(accAddr+skipInc(i),2^accWidth);
        accAddr90 = modulo(accAddr90+skipInc(i),2^accWidth);

        //quantize
        //plus one required since scilab address begins at 1
        qntAddr = floor(accAddr*(2^-lsbWidth))+1;
        qntAddr90 = floor(accAddr90*(2^-lsbWidth))+1;

        yLow = tb(qntAddr);
        yLow90 = tb(qntAddr90);
        //temp
        yLow_tap(i) = yLow;
        yLow90_tap(i) = yLow90;

        //end osc
    end
    //__**end main loop**__

    //noise reduction - Taylor Series noise reduction (future)

    scf(1)
    plot(tb,'black')
    xlabel("$0\le x\le 4096$","fontsize",4,"color","black");
    ylabel("Wavetable Wave Shape","fontsize",4,"color","black");
    title("Wavetable Wave Shape","color","black","fontsize",4);


endfunction

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//function: odmk_osc5 - end
////////////////////////////////////////////////////////////////


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//function: odmk_osc5Stereo - begin
////////////////////////////////////////////////////////////////

//wavetable oscillator function

function [odmkOsc, odmkOsc90, odmkSqrPulse] = odmk_osc5Stereo(numSamples, Fs, tableDepth, shape, foArray, phArray)
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // odmk_osc5Stereo: Stereo channel complex wavetable oscillator
    //
    // The output frequency can be fixed or variable
    // When foArray is a single scalar value, the output freq. is fixed
    // When foArray is an array of length=numSamples, the output freq. varies each step
    //
    // The output phase can be fixed or variable
    // When phaseArray is a single scalar value, the output phase is fixed
    // When phaseArray is an array of length=numSamples, the output phase varies each step 
    //
    // Real time controls:
    // frequency ; phase
    //
    // data format: input foMatrix, phaseMatrix => (numSamples,1)
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    tb = zeros(tableDepth,1);
    tb = tablegen(shape,tableDepth);
    
    //check inputs
    if (length(foArray(:,1)) == 1) & (length(phArray(:,1)) == 1) then
        foArr = foArray;
        phArr = phArray;
        skipInc = [0, 0];
        phOffset = [0, 0];
    elseif (length(foArray) == numSamples) | (length(phArray) == numSamples) then
        foArr = foArray;
        phArr = phArray;
        skipInc = zeros(numSamples,2);
        phOffset = zeros(numSamples,2);
    else
        error('length of phArray must be 1 or match numSamples')
    end


    //***initialize***

    //init osc output vectors
    odmkOsc = zeros(numSamples,2);
    odmkOsc90 = zeros(numSamples,2);
    odmkSqrPulse = zeros(numSamples,2);
    
    accWidth = 48;
    qntWidth = ceil(log2(tableDepth));
    lsbWidth = accWidth-qntWidth;    //expect 36
    //mprintf('lsb_width %i\n\n',lsb_width)
    
    //scale foArray in place]
    foArr = (2^lsbWidth)*foArr;
       
    //phase accumulator increment value (skipping value)
    //skipInc = floor(fo*(2^N/fs)+0.5);
    //skipInc = floor(foscale*(2^qntWidth)/fs+0.5);
    //skipInc = floor((foscale*(2^qntWidth))/fs);
    skipInc = ((2^qntWidth).*foArr)./fs; 
        
    //scale phase offset
    //converts radians into a scaled phase offset to be added to the output of the acc
    //dependent on the table depth
    phOffset = round(modulo(((phArr.*tableDepth)./(2*%pi)),tableDepth));
        
    //used to add a 90 deg offset for complex sinusoid generation
    offset90 = tableDepth/4;     
    
    //**main loop**
    //___::((Interpolated WaveTable))::___
    //generates main osc waveform out, 90deg shifted out (sin/cos), square pulse out
    
    //initialize variables so that the OSC starts at the beginning of the table
    accAddr = phOffset(1,:)*2^lsbWidth;    //init to 'zero plus phase shift'
    accAddr90 = (offset90+phOffset(1,:))*2^lsbWidth;    //init to 'zero plus offset90 plus phase shift'
    qntAddr = 1+phOffset(1,:);
    qntAddr90 = qntAddr+offset90;
    for h = 1:2
        yLow(h) = tb(qntAddr(h));
        yLow90(h) = tb(qntAddr90(h));
    end
    yLow = yLow';
    yLow90 = yLow90';

    
    for i=1:numSamples;
        //osc main loop
    
        //switch between fixed or continuously variable frequency 
           
        if (length(skipInc(:,1))==1) then            
        //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        //generate fixed frequency oscillator output
        /////////////////////////////////////////////////////        
    
            //accAddr plus one
            accAddrP1 = accAddr+(2^lsbWidth);
            accAddr90P1 = accAddr90+(2^lsbWidth);
       
            
            //Stereo Processing - calculate new output sample
            for j = 1:2
                //avoid addr overflow for linear interpolation addressing - should investigate
                if (qntAddr(j)==tableDepth) then    
                    yHigh(j) = tb(qntAddr(j));
                    yHigh90(j) = tb(qntAddr90(j)+1);
                    odmkOsc(i,j) = tb(tableDepth);
                    //linear interpolation
                    odmkOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*((accAddr90P1(j)-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                    yHigh90_tap(i,j) = yHigh90(j);
                elseif (qntAddr90(j)==tableDepth) then
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j));
                    //linear interpolation
                    odmkOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*((accAddrP1(j)-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkOsc90(i,j) = tb(tableDepth);
                else
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j)+1);

                    odmkOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*((accAddrP1(j)-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*((accAddr90P1(j)-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                end 
            
                //generate square pulse output
                if odmkOsc(i,j) >= 0 then
                    odmkSqrPulse(i,j) = 1;
                else
                    odmkSqrPulse(i,j) = 0;
                end
                        
                //phase accumulator
                accAddr(j) = modulo(accAddr(j)+skipInc(j),2^accWidth);
                accAddr90(j) = modulo(accAddr90(j)+skipInc(j),2^accWidth);
            
                //quantize
                //plus one required since scilab address begins at 1
                qntAddr(j) = floor(accAddr(j)*(2^-lsbWidth))+1;
                qntAddr90(j) = floor(accAddr90(j)*(2^-lsbWidth))+1;
            
                yLow(j) = tb(qntAddr(j));
                yLow90(j) = tb(qntAddr90(j));
            end
            
        else
        //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        //generate variable frequency oscillator output
        /////////////////////////////////////////////////////              

            //accAddr plus one
            accAddrP1 = accAddr+(2^lsbWidth);
            accAddr90P1 = accAddr90+(2^lsbWidth);

            //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            //generate oscillator output
            /////////////////////////////////////////////////////

            //Stereo Processing - calculate new output sample
            for j = 1:2
                //avoid addr overflow for linear interpolation addressing - should investigate
                if (qntAddr(j)==tableDepth) then    
                    yHigh(j) = tb(qntAddr(j));
                    yHigh90(j) = tb(qntAddr90(j)+1);
                    odmkOsc(i,j) = tb(tableDepth);
                    //linear interpolation
                    odmkOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*(((accAddr90(i,j)+(2^lsbWidth))-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                    yHigh90_tap(i,j) = yHigh90(j);
                elseif (qntAddr90(j)==tableDepth) then
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j));
                    //linear interpolation
                    odmkOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*(((accAddr(i,j)+(2^lsbWidth))-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkOsc90(i,j) = tb(tableDepth);
                else
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j)+1);
                    odmkOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*((accAddrP1(j)-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*((accAddr90P1(j)-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                end 

                //generate square pulse output
                if odmkOsc(i,j) >= 0 then
                    odmkSqrPulse(i,j) = 1;
                else
                    odmkSqrPulse(i,j) = 0;
                end


                //phase accumulator
                accAddr(j) = modulo(accAddr(j)+skipInc(i,j),2^accWidth);
                accAddr90(j) = modulo(accAddr90(j)+skipInc(i,j),2^accWidth);                     

                //quantize
                //plus one required since scilab address begins at 1
                qntAddr(j) = floor(accAddr(j)*(2^-lsbWidth))+1;
                qntAddr90(j) = floor(accAddr90(j)*(2^-lsbWidth))+1;

                yLow(j) = tb(qntAddr(j));
                yLow90(j) = tb(qntAddr90(j));
                //temp
                yLow_tap(i,j) = yLow(j);
                yLow90_tap(i,j) = yLow90(j);
            end
        end            
    end
    //end osc
    //__**end main loop**__

    //noise reduction - Taylor Series noise reduction (future)

    scf(1)
    plot(tb,'black')
    xlabel("$0\le x\le 4096$","fontsize",4,"color","black");
    ylabel("Wavetable Wave Shape","fontsize",4,"color","black");
    title("Wavetable Wave Shape","color","black","fontsize",4);

endfunction

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//function: odmk_osc5Stereo - end
////////////////////////////////////////////////////////////////


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//function: odmk_odmk_multiOsc - begin
////////////////////////////////////////////////////////////////

//wavetable oscillator function

function [odmkMultiOsc, odmkMultiOsc90, odmkMultiSqrPulse] = odmk_multiOsc(numSamples, Fs, tableDepth, shape, foArray, phArray)
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // odmk_multiOsc: multi-channel complex wavetable oscillator
    //
    //::Inputs:: 
    //numSamples: number os samples in output arrays
    //Fs: sampling frequency
    //tableDepth: depth of Wavetable
    //shape: Wavetable output waveform shape
    //
    //foArray: frequency Array => [static freq -> length(foArray) = 1] or [variable freq -> length(foArray) = numSamples]
    //phArray: phase Array => [static phase -> length(phArray) = 1] or [variable phase -> length(phArray) = numSamples]  
    //
    //::Real time controls::
    //frequency ; phase
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

    tb = zeros(tableDepth,1);
    tb = tablegen(shape,tableDepth);
    
    nChan = length(foArray(1,:));
    
    //check inputs (assume all channels are the same length)
    //*** fix for independent access to freq and phase modulation
    if (length(foArray(:,1)) == 1) & (length(phArray(:,1)) == 1) then
        foArr = foArray;
        phArr = phArray;
        skipInc = zeros(1,nChan);
        phOffset = zeros(1,nChan);
    elseif (length(foArray) == numSamples) | (length(phArray) == numSamples) then
        foArr = foArray;
        phArr = phArray;
        skipInc = zeros(numSamples,nChan);
        phOffset = zeros(numSamples,nChan);
    else
        error('length of phArray must be 1 or match numSamples')
    end
    mprintf('odmk_multiOsc->nChan %i\n\n',nChan)
    
    //***initialize***
    
    //init osc output vectors
    odmkOsc = zeros(numSamples,nChan);
    odmkOsc90 = zeros(numSamples,nChan);
    odmkSqrPulse = zeros(numSamples,nChan);
    
    accWidth = 48;
    qntWidth = ceil(log2(tableDepth));
    lsbWidth = accWidth-qntWidth;    //expect 36
    //mprintf('lsb_width %i\n\n',lsb_width)
    
    //scale foArray in place]
    foArr = (2^lsbWidth)*foArr;
        
    //phase accumulator increment value (skipping value)
    //skipInc = floor(fo*(2^N/fs)+0.5);
    //skipInc = floor(foscale*(2^qntWidth)/fs+0.5);
    //skipInc = floor((foscale*(2^qntWidth))/fs);
    skipInc = ((2^qntWidth).*foArr)./fs; 
        
    //scale phase offset
    //converts radians into a scaled phase offset to be added to the output of the acc
    //dependent on the table depth
    phOffset = round(modulo(((phArr.*tableDepth)./(2*%pi)),tableDepth));

    //used to add a 90 deg offset for complex sinusoid generation
    offset90 = tableDepth/4; 


    //**main loop**
    //___::((Interpolated WaveTable))::___
    //generates main osc waveform out, 90deg shifted out (sin/cos), square pulse out

    //initialize variables so that the OSC starts at the beginning of the table
    accAddr = phOffset(1,:)*2^lsbWidth;    //init to 'zero plus phase shift'
    accAddr90 = (offset90+phOffset(1,:))*2^lsbWidth;    //init to 'zero plus offset90 plus phase shift'
    //must use modulus here to roll qntaddr
    qntAddr = modulo(1+phOffset(1,:),tableDepth);
    qntAddr90 = modulo(qntAddr+offset90,tableDepth);    

    for h = 1:nChan
        yLow(h) = tb(qntAddr(h));
        //mprintf('odmk_multiOsc->qntAddr %i\n\n',qntAddr90(h))
        //mprintf('odmk_multiOsc->qntAddr90 len %i\n\n',qntAddr90(h))
        yLow90(h) = tb(qntAddr90(h));
    end
    yLow = yLow';
    yLow90 = yLow90';


    for i=1:numSamples;
        //osc main loop

        //switch between fixed or continuously variable frequency 

        if (length(skipInc(:,1))==1) then            
        //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        //generate fixed frequency oscillator output
        /////////////////////////////////////////////////////        
    
            //accAddr plus one
            accAddrP1 = accAddr+(2^lsbWidth);
            accAddr90P1 = accAddr90+(2^lsbWidth);

            //Stereo Processing - calculate new output sample
            for j = 1:nChan
                //avoid addr overflow for linear interpolation addressing - should investigate
                if (qntAddr(j)==tableDepth) then    
                    yHigh(j) = tb(qntAddr(j));
                    yHigh90(j) = tb(qntAddr90(j)+1);
                    odmkMultiOsc(i,j) = tb(tableDepth);
                    //linear interpolation
                    odmkMultiOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*((accAddr90P1(j)-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                    yHigh90_tap(i,j) = yHigh90(j);
                elseif (qntAddr90(j)==tableDepth) then
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j));
                    //linear interpolation
                    odmkMultiOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*((accAddrP1(j)-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkMultiOsc90(i,j) = tb(tableDepth);
                else
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j)+1);
                    odmkMultiOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*((accAddrP1(j)-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkMultiOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*((accAddr90P1(j)-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                end 

                //generate square pulse output
                if odmkMultiOsc(i,j) >= 0 then
                    odmkMultiSqrPulse(i,j) = 1;
                else
                    odmkMultiSqrPulse(i,j) = 0;
                end

                //phase accumulator
                accAddr(j) = modulo(accAddr(j)+skipInc(j),2^accWidth);
                accAddr90(j) = modulo(accAddr90(j)+skipInc(j),2^accWidth);

                //quantize
                //plus one required since scilab address begins at 1
                qntAddr(j) = floor(accAddr(j)*(2^-lsbWidth))+1;
                qntAddr90(j) = floor(accAddr90(j)*(2^-lsbWidth))+1;

                yLow(j) = tb(qntAddr(j));
                yLow90(j) = tb(qntAddr90(j));
            end

        else
        //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        //generate variable frequency oscillator output
        /////////////////////////////////////////////////////              

            //accAddr plus one
            accAddrP1 = accAddr+(2^lsbWidth);
            accAddr90P1 = accAddr90+(2^lsbWidth);

            //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            //generate oscillator output
            /////////////////////////////////////////////////////

            //Parallel Processing - calculate new output sample for each nChan channel
            for j = 1:nChan
                //avoid addr overflow for linear interpolation addressing - should investigate
                if (qntAddr(j)==tableDepth) then    
                    yHigh(j) = tb(qntAddr(j));
                    yHigh90(j) = tb(qntAddr90(j)+1);
                    odmkMultiOsc(i,j) = tb(tableDepth);
                    //linear interpolation
                    odmkMultiOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*(((accAddr90(i,j)+(2^lsbWidth))-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                    yHigh90_tap(i,j) = yHigh90(j);
                elseif (qntAddr90(j)==tableDepth) then
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j));
                    //linear interpolation
                    odmkMultiOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*(((accAddr(i,j)+(2^lsbWidth))-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkMultiOsc90(i,j) = tb(tableDepth);
                else
                    yHigh(j) = tb(qntAddr(j)+1);
                    yHigh90(j) = tb(qntAddr90(j)+1);
                    odmkMultiOsc(i,j) = yLow(j) + (yHigh(j)-yLow(j))*((accAddrP1(j)-(qntAddr(j)*2^lsbWidth))*2^-lsbWidth);
                    odmkMultiOsc90(i,j) = yLow90(j) + (yHigh90(j)-yLow90(j))*((accAddr90P1(j)-(qntAddr90(j)*2^lsbWidth))*2^-lsbWidth);
                end 

                //generate square pulse output
                if odmkMultiOsc(i,j) >= 0 then
                    odmkMultiSqrPulse(i,j) = 1;
                else
                    odmkMultiSqrPulse(i,j) = 0;
                end

                //phase accumulator
                accAddr(j) = modulo(accAddr(j)+skipInc(i,j),2^accWidth);
                accAddr90(j) = modulo(accAddr90(j)+skipInc(i,j),2^accWidth);

                //quantize
                //plus one required since scilab address begins at 1
                qntAddr(j) = floor(accAddr(j)*(2^-lsbWidth))+1;
                qntAddr90(j) = floor(accAddr90(j)*(2^-lsbWidth))+1;

                yLow(j) = tb(qntAddr(j));
                yLow90(j) = tb(qntAddr90(j));
                //temp
                yLow_tap(i,j) = yLow(j);
                yLow90_tap(i,j) = yLow90(j);
            end
        end            
    end
    //end osc
    //__**end main loop**__
    
    //noise reduction - Taylor Series noise reduction (future)

    scf(1)
    plot(tb,'black')
    xlabel("$0\le x\le 4096$","fontsize",4,"color","black");
    ylabel("Wavetable Wave Shape","fontsize",4,"color","black");
    title("Wavetable Wave Shape","color","black","fontsize",4);


endfunction

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//function: odmk_multiOsc - end
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : DDS functions
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
//main: osc2
////////////////////////////////////////////////////////////////


//oscillator outputs sin and cos waveforms - uses 1/4 wave symmetry
//and linear interpolation


//DDS implementation outputs sin, cos. 
//fo = output frequency
//fs = sampling frequency
//td = table depth

//fo = fs/(2^N)*(delta(acc))
//delta(acc) = floor(fo*((2^N)/fs)+0.5)

//sampling frequency
fs = 44100;

//Number of Iterations => defines length of output wav
//numSamples = 131072;    //4096 8192 32768 65536 131072
//numSamples = 32768;    //works w/ single osc
numSamples = 8192;

//define length of FFT for analysis
Nfft = 4096;    //4096; //define length of FFT for analysis


//table depth = number of samples used to store waveform
tableDepth = 4096;   

//waveform shape - 
//1=sin
//2=tri
//3=saw-up
//4=saw-down
//5=chebychev
//6=pulse1
//7=pulse2
//8=pulse3
//9=pulse4
shape = 1;    
 
 
 
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// example osc 1: mono oscillator to mono-stereo output
///////////////////////////////////////////////////////////////

//--<<mono>>--//

//output frequency: fo = fs/(2^N)*(phacc))
//foArray = 528;    //440 528 Hz
foArray = 5000;
//foArray = linspace(930,130,numSamples)'

//define phase value (0-2pi radians)
//phase = [0,%pi/4,%pi/2,3*%pi/4,%pi];
phArray = 0;    //%pi/2;
//phArray = linspace(0,2*%pi,numSamples)'


//init wave file outputs
odmkOscWav = zeros(numSamples,2);
odmkSqrWav = zeros(numSamples,2);


[odmkOsc,odmkOsc90,odmkSqrPulse] = odmk_osc5(numSamples, fs, tableDepth, shape, foArray, phArray) 

//format for writing .wav -> [n,2]
odmkOscWav = [odmkOsc,odmkOsc];
odmkSqrWav = [odmkSqrPulse,odmkSqrPulse];

//analysis
//check spectrum
if numSamples > Nfft then
    //odmkOscFFT = fft(odmkOsc(1:Nfft));
    odmkOscFFT = fft(odmkOsc(1:Nfft));
    odmkOscMag = sqrt(abs(odmkOscFFT(1:Nfft).*odmkOscFFT(1:Nfft)));
else
    odmkOscFFT = fft(odmkOsc);
    odmkOscMag = sqrt(abs(odmkOscFFT.*odmkOscFFT));
end

////////////////////////////////////////////////////////////////
//end: example osc 1
////////////////////////////////////////////////////////////////


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// example osc 2: stereo oscillator (fixed freq & phase)
/////////////////    ;'//////////////////////////////////////////////


//--<<stereo>>--//

//output frequency: fo = fs/(2^N)*(phacc))
//foArray = 528;    //440 528 Hz

//foFixStereo = [528, 528];
foFixStereo = [5000, 10000];

//continually variable Frequency array:
//foArrayStereo = zeros(numSamples,2);
//foArrayStereo(:,1) = linspace(930,130,numSamples)'
//foArrayStereo(:,2) = linspace(130,930,numSamples)'


//define phase value (0-2pi radians)
//phase = [0,%pi/4,%pi/2,3*%pi/4,%pi];

//phFixStereo = [0, 0];
phFixStereo = [%pi/4, 3*%pi/4];

//continually variable Phase array:
//phArrayStereo = zeros(numSamples,2);
//foArrayStereo(:,1) = linspace(930,130,numSamples)'
//foArrayStereo(:,2) = linspace(130,930,numSamples)'


//init wave file outputs
odmkOscStereo = zeros(numSamples,2);
odmkOsc90Stereo = zeros(numSamples,2);
odmkSqrPulseStereo = zeros(numSamples,2);

[odmkOscStereo,odmkOsc90Stereo,odmkSqrPulseStereo] = odmk_osc5Stereo(numSamples, fs, tableDepth, shape, foFixStereo, phFixStereo) 


//create mono signal
odmkOscStereo_L = odmkOscStereo(:,1);
//odmkOsc90Stereo_L = odmkOsc90Stereo(:,1);
//odmkSqrPulseStereo_L = odmkSqrPulseStereo(:,1);


////format for writing .wav -> [n,2]
//odmkOscWav = [odmkOsc,odmkOsc];
//odmkSqrWav = [odmkSqrPulse,odmkSqrPulse];

//analysis
//check spectrum
if numSamples > Nfft then
    //check Left side -> real only FFT
    odmkOscStereoFFT_L = fft(odmkOscStereo_L(1:Nfft));
    odmkOscStereoMag_L = sqrt(abs(odmkOscStereoFFT_L(1:Nfft).*odmkOscStereoFFT_L(1:Nfft)));
    odmkOscStereoFFT_R = fft(odmkOscStereo(1:Nfft,2));
    odmkOscStereoMag_R = sqrt(abs(odmkOscStereoFFT_R(1:Nfft).*odmkOscStereoFFT_R(1:Nfft)));
else
    odmkOscStereoFFT_L = fft(odmkOscStereo_L);
    odmkOscStereoMag_L = sqrt(abs(odmkOscStereoFFT_L.*odmkOscStereoMag_L));
    odmkOscStereoFFT_R = fft(odmkOscStereo(:,2));
    odmkOscStereoMag_R = sqrt(abs(odmkOscStereoFFT_R.*odmkOscStereoMag_R));    
end


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// example osc 3: multi channel oscillator (fixed freq & phase)
/////////////////    ;'//////////////////////////////////////////////


//--<<multi channel>>--//

//output frequency: fo = fs/(2^N)*(phacc))
//foArray = 528;    //440 528 Hz

//foFixStereo = [528, 528];
foMFix = [500, 2500, 4000, 5500, 7000, 9500];

//continually variable Frequency array:
//foArrayStereo = zeros(numSamples,2);
//foArrayStereo(:,1) = linspace(930,130,numSamples)'
//foArrayStereo(:,2) = linspace(130,930,numSamples)'


//define phase value (0-2pi radians)
//phase = [0,%pi/4,%pi/2,3*%pi/4,%pi];

phMFix = [0, 0, 0, 0, 0, 0];
//phMFix = [0, %pi/3, 2*%pi/3, %pi, 4*%pi/3, 5*%pi/3];
//phFixStereo = [%pi/4, 3*%pi/4];

//continually variable Phase array:
//phArrayStereo = zeros(numSamples,2);
//foArrayStereo(:,1) = linspace(930,130,numSamples)'
//foArrayStereo(:,2) = linspace(130,930,numSamples)'

Mchan = length(foMFix(1,:))

//init wave file outputs
odmkMOsc = zeros(numSamples,Mchan);
odmkMOsc90 = zeros(numSamples,Mchan);
odmkMSqrPulse = zeros(numSamples,Mchan);

[odmkMOsc,odmkMOsc90,odmkMSqrPulse] = odmk_multiOsc(numSamples, fs, tableDepth, shape, foMFix, phMFix) 


//create mono signal
odmkMOsc_L = odmkMOsc(:,1);
//odmkOsc90Stereo_L = odmkOsc90Stereo(:,1);
//odmkSqrPulseStereo_L = odmkSqrPulseStereo(:,1);


////format for writing .wav -> [n,2]
//odmkOscWav = [odmkOsc,odmkOsc];
//odmkSqrWav = [odmkSqrPulse,odmkSqrPulse];


//parallel wav spectrum analysis
if numSamples > Nfft then
    for k=1:Mchan
        odmkMOscFFT(:,k) = fft(odmkMOsc((1:Nfft),k));
        odmkMOscMag(:,k) = sqrt(abs(odmkMOscFFT((1:Nfft),k).*odmkMOscFFT((1:Nfft),k)));
    end
else
    for j=1:MChan
        odmkMOscFFT(:,j) = fft(odmkMOsc(:,j))
        odmkMOscMag(:,j) = sqrt(abs(odmkMOscFFT(:,j).*odmkMOscFFT(:,j)))
    end
end

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

//create frequency axis vectors
freqAxis=fs*(0:(Nfft-1))/Nfft;

freqAxisx=fs*(1:(Nfft))/Nfft;
freqAxis_log = log(freqAxisx(1:Nfft/2));
 
freqAxis_zp=fs*(0:(2*Nfft-1))/(2*Nfft);

freqAxis_zpx=fs*(1:(2*Nfft))/Nfft;
freqAxis_zplog = log(freqAxis_zpx);

//scf(1)

truncAmt = 46;

scf(2)
plot(odmkOsc(1:truncAmt),"blue")
plot(odmkOsc90(1:truncAmt),"cyan")
plot(odmkSqrPulse(1:truncAmt),"green")
title("$odmk\_osc4: odmkOsc\ (blue)\ /\ odmkOsc90\ (magenta)\ /\ odmkSqrPulse\ (cyan)$","color","black","fontsize",3)
xlabel 'time'
ylabel 'amplitude'


scf(3)
plot(odmkOscStereo_L(1:truncAmt),"red")
plot(odmkOscStereo(1:truncAmt,2),"magenta")
title("$odmk\_osc4Stereo: odmkOsc90Stereo\_L\ (red)\ /\ odmkOsc90Stereo\_R\ (magenta)$","color","black","fontsize",3)
xlabel 'time'
ylabel 'amplitude'


scf(4)
plot(odmkOscStereo_L(1:truncAmt),"red")
plot(odmkOsc90Stereo(1:truncAmt,1),"magenta")
plot(odmkSqrPulseStereo(1:truncAmt,1),"black")
title("$odmk\_osc4Stereo: odmkOsc90Stereo\_L\ (red)\ /\ odmkOsc90Stereo\_L (magenta)\ /\ odmkSqrPulseStereo\_L\ (black)$","color","black","fontsize",3)
xlabel 'time'
ylabel 'amplitude'




//scf(5)
//plot2d(freqAxis_log,odmkOscMag(1:Nfft/2),logflag='ln')
//plot2d(freqAxis_log,odmkOscStereoMag_L(1:Nfft/2),logflag='ln')
//title 'odmk_osc5 vs odmk_osc5Stereo: Left Output FFT Magnitude'
//xlabel 'frequency'
//ylabel 'amplitude'



scf(6)
//plot2d(freqAxis,20*log10(odmkOscMag),style=[3]);
plot2d(freqAxis,20*log10(odmkOscStereoMag_L),style=[3]);
plot2d(freqAxis,20*log10(odmkOscStereoMag_R),style=[4]); 
title("$odmk\_oscStereo\ FFT\ Magnitude,\ (blue=L\ ;\ green=R)$","color","black","fontsize",3);
xlabel("$Linear\ Frequency : 0\le x\le fs\ (Sample Rate)$","fontsize",3,"color","black");
ylabel("$Log-Scaled\ Magnitude$","fontsize",3,"color","black"); 
a=gca(); // Handle on axes entity
a.grid=[2,2];



scf(7)
for h=1:Mchan
    plot2d(freqAxis(1:truncAmt),odmkMOsc(1:truncAmt,h),style=[h]);
end 
title("$odmk\_multiOsc\ n\ channel\ FFT\ Magnitude$","color","black","fontsize",3);
xlabel 'time'
ylabel 'amplitude'
a=gca(); // Handle on axes entity
a.grid=[2,2];


scf(8)
for h=1:Mchan
    plot2d(freqAxis,odmkMOscMag(:,h),style=[h]);
end 
title("$odmk\_multiOsc\ n\ channel\ FFT\ Magnitude$","color","black","fontsize",3);
xlabel("$Linear Frequency : 0\le x\le fs (Sample Rate)$","fontsize",3,"color","black");
ylabel("$Linear Magnitude$","fontsize",3,"color","black"); 
a=gca(); // Handle on axes entity
a.grid=[2,2];


scf(9)
for h=1:Mchan
    plot2d(freqAxis,20*log10(odmkMOscMag(:,h)),style=[h]);
end 
title("$odmk\_multiOsc\ n\ channel\ FFT\ Magnitude$","color","black","fontsize",3);
xlabel("$Linear Frequency : 0\le x\le fs (Sample Rate)$","fontsize",3,"color","black");
ylabel("$Log-Scaled\ Magnitude$","fontsize",3,"color","black"); 
a=gca(); // Handle on axes entity
a.grid=[2,2];



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : plotting
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : write into wav files:
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////

//bits = 24;
//
////create a new audio file
////nm='unicorn_clone1';    32bit IEEE float
//nm_out1='Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/osc/odmk_osc';
//xwavwrite(odmkOscWav,fs,bits,nm_out1);
//mprintf(' <xwavwrite out> wrote .wav file %s\n\n',nm_out1)
//
////create a new audio file
////nm='unicorn_clone1';    32bit IEEE float
//nm_out2='Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/osc/odmk_sqrpulse';
//xwavwrite(odmkSqrWav,fs,bits,nm_out2);
//mprintf(' <xwavwrite out> wrote .wav file %s\n\n',nm_out2)

///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : write into wav files:
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////