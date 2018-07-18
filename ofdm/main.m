clear
clc
warning('off','all')
% fadingModel = 0 is no fading
% fadingModel = 1 is uniform profile
% fadingModel = 11 is uniform profile with constant gain ( for testing)
% fadingModel = 2 is Exponential profile
% fadingModel = 22 is Exponential profile with constant gain ( for testing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%change the fading profile 
fadingModel=1;              %1: Uniform, 2: Exponential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%change delay spread
maxDelaySpreadInSamples=5;   %5,15,25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNRperCarrierdB = 10;
startSNRperSubcarrierdB=0.0;
stepSNRperSubcarrierdB=4.0;
stopSNRperSubcarrierdB=30.0;

snrArray=[startSNRperSubcarrierdB:stepSNRperSubcarrierdB:stopSNRperSubcarrierdB];

fftSize=128;
numGuardLeft=11;
numGuardRight=10;
cyclicPrefixRatio=0.125;
numUsedSubcariers=fftSize-numGuardLeft-numGuardRight-1;
numOFDMsymbols=2;

% for the time being there are no pilots
numDataSubcariers = numUsedSubcariers;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%change modulation
modulationIndex=4;   %1 BPSK, 2 QPSK, 4 16QAM, 6 64QAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxSimulationFrames=100000;
lowestNeededFER=1.0E-4;
sufficientChannelFrameError=200;
sufficientCodedFrameError=200; 

%********************************
% Covolutional Coding parameters
%********************************
codeRate=0.5;

K=7;
Generator1=171;
Generator2=133;
Generator3=165;
t = poly2trellis(K, [Generator1 Generator2]); % Define trellis.
numTailBits=K-1;
%********************************
% uncoded Eb/No 
uncodedEb_NodBarray=snrArray - 10.0*log10(modulationIndex) ;
       
% Coded Eb/No in dB
codedEb_NodBarray=uncodedEb_NodBarray - 10.0*log10(codeRate) ;


%***************************************
% Prepare what will be written on the figure
%=============================
if fadingModel==0
    TXTchannelModel='AWGN'
elseif fadingModel==1
    TXTchannelModel='Fading Uniform'
elseif fadingModel==2
    TXTchannelModel='Fading Exponential'
else
    TXTchannelModel='Fading ?'   
end
%=============================

%=============================
if modulationIndex==1
    TXTmodType='BPSK'
elseif modulationIndex==2
    TXTmodType='QPSK'
elseif modulationIndex==4
    TXTmodType='16QAM'
elseif modulationIndex==6
    TXTmodType='64QAM'
else
    TXTmodType='mod?';
end
%=============================
% any permutation number now
permutation=1;
if permutation==1
    TXTpermutation=[num2str(fftSize) ' PUSC']
elseif permutation==2
    TXTpermutation=[num2str(fftSize) ' FUSC']
end

% number of source bits must be integer
numSourceBits=numDataSubcariers*modulationIndex*codeRate*numOFDMsymbols - numTailBits ;

%=======================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over enough number of slots for simaultion
numberChannelBitErrors=zeros(1,length(snrArray));
numberCodedBitErrors=zeros(1,length(snrArray));

numberOfChannelFrameErrors=zeros(1,length(snrArray));
numberOfCodedFrameErrors=zeros(1,length(snrArray));
numSimulatedBits=zeros(1,length(snrArray));
numSimulatedSourceBits=zeros(1,length(snrArray));
numSimulatedFrames=zeros(1,length(snrArray));
snrDone=zeros(1,length(snrArray));

for frameNumber = 1 : maxSimulationFrames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(frameNumber,10000) ==0 , frameNumber , end ;
    if sum(snrDone)==length(snrDone) , break, end ;  % all SNR are done
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Let us start with the data
    % random Bits of desired user
    
    sourceBits=round(rand(1, numSourceBits)) ;
    msg = [sourceBits , zeros(1, numTailBits)]; % Random data
    
    % Convolutional Encoder
    codedBits = convenc(msg,t); % Encode the data
    numChannelBits=length(codedBits);
    
    % Random Interleaver for now till we make the standard
    %=====================================================
    interleaverState=randint(1,1,999);
    interleavedCodedBits = interleaverDeinterleaver(codedBits, 0, interleaverState);
    
    % Modulation
     txSymbols= mapping( interleavedCodedBits, modulationIndex );
     
     % make the symbols in matrix format. Each column=1 OFDM symbol
     ofdmInput=reshape(txSymbols,numUsedSubcariers , numOFDMsymbols); 
     
     % Perform OFDM
     txSamples= ofdmTx(ofdmInput, numGuardLeft, numGuardRight, fftSize, cyclicPrefixRatio );

     % Apply fading
    [fadedSamples, gainImpulseResponse]=ApplyFading(txSamples, fadingModel,...
                                           maxDelaySpreadInSamples); 

    % generate unit variance complex Gaussian noise
    noiseSamples=1/sqrt(2)*(randn(size(fadedSamples))+j*randn(size(fadedSamples))) ;
    %===============================================

     
    
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is snr per subraiier in dB
    for snrdBcounter = 1 : length(snrArray)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================================================================
%=================================================================
     % Stop running snr values if not getting enough errors 
     % After "10/lowestNeededFER" frames, FER is below "lowestNeededFER/10"
     % Then stop this SNR value
        minFramesToRun=1000;
        
        FERtest=1;
        if frameNumber > 10 
            FERtest=numberOfCodedFrameErrors(snrdBcounter) / ...
                                numSimulatedFrames(snrdBcounter);
        end

        % condition1 means enough simualtion and can't get needed FER
        condition1=frameNumber > max(minFramesToRun, 100/lowestNeededFER) & ...
           FERtest < lowestNeededFER ;
        
       % condition2 means that this SNR will never get the needed FER
        condition2=frameNumber > max(minFramesToRun, 10/lowestNeededFER) & ...
           FERtest < lowestNeededFER/10 ;
       
        if condition1==1 | condition2==1
     
            if snrDone(snrdBcounter) == 0
                snrDone(snrdBcounter)=1;
                'EB/No completed is:'
                snrPerSubcarrierdB=snrArray(snrdBcounter)
                channelBitErrorRate= numberChannelBitErrors(snrdBcounter) / ...
                                         numSimulatedBits(snrdBcounter)

                codedBitErrorRate= numberCodedBitErrors(snrdBcounter) / ...
                                         numSimulatedSourceBits(snrdBcounter)
                FER=FERtest
            end

            continue  % Skip this SNR value
        end
%=================================================================
%=================================================================
     % Stop running Eb/No values after getting enough errors 
        if numberOfChannelFrameErrors(snrdBcounter) > sufficientChannelFrameError & ...
             numberOfCodedFrameErrors(snrdBcounter) > sufficientCodedFrameError 

            if snrDone(snrdBcounter) == 0
                snrDone(snrdBcounter)=1;
                'EB/No completed is:'
                snrPerSubcarrierdB=snrArray(snrdBcounter)
                channelBitErrorRate= numberChannelBitErrors(snrdBcounter) / ...
                                         numSimulatedBits(snrdBcounter)

                %BER_BPSK=0.5*erfc(sqrt(snrArrayRatio(snrdBcounter)))  % for testing                   

                codedBitErrorRate= numberCodedBitErrors(snrdBcounter) / ...
                                         numSimulatedSourceBits(snrdBcounter)
                                     
                FER=FERtest
            end
            continue  % Skip this snr value
        end
%=================================================================
%=================================================================
        snrPerSubcarrierdB=snrArray(snrdBcounter);  % Eb/No of channel bits
        
        numSimulatedBits(snrdBcounter)=numSimulatedBits(snrdBcounter)+ ...
            numChannelBits;

        numSimulatedSourceBits(snrdBcounter)=numSimulatedSourceBits(snrdBcounter)+ ...
            numSourceBits;
        
        numSimulatedFrames(snrdBcounter)=numSimulatedFrames(snrdBcounter)+ 1;

        % SNR per carrier in dB. Each carrier carries one modulation symbol
        snrPerCarrierRatio=10.0^(snrPerSubcarrierdB/10);

        
    % Noise is added here according to the SNR per subcarrier
    noiseStandardDeviation=sqrt(1.0/fftSize/snrPerCarrierRatio);
    noisySamples=fadedSamples+noiseSamples*noiseStandardDeviation;
    
    % Apply noisy faded samples to OFDM receiver
    ofdmRxOutput= ofdmRx(noisySamples, numGuardLeft, numGuardRight, fftSize, cyclicPrefixRatio );

    % channel estiamtion is done here per OFDM symbol. Now it is ideal
    channelGainsFrequency = channelEstimation( gainImpulseResponse, ...
                                    numGuardLeft, numGuardRight, fftSize );

    %Equalization is done here
    equalizedSymbols  =channelEqualizer(ofdmRxOutput, channelGainsFrequency);
        
    % demapping the received symbols, get soft bits
    % For QAM, it assumes constellations +/-1, +/-3, .. etc
    % we neead the magnitude square of channel gains for QAM
    gain2=abs(channelGainsFrequency).^2;
    rxSoftBits = demapping( equalizedSymbols, modulationIndex, gain2 );
       
        
        % hard bits for channel BER
        detectedChannelBits=(1-sign(rxSoftBits))/2;
        
        % Let us deinterleave
        %====================
        deinterleavedSoftBits = interleaverDeinterleaver(rxSoftBits, 1, interleaverState);
        tblen = 5*K;  % Traceback length
        
        % Viterbi decoder
        %=================
        decodedBits = vitdec(deinterleavedSoftBits,t,tblen,'term','unquant'); % Decode.
        
        % do not take tail bits in BER calculation
        numCodedBitErrorsInFrame=sum(bitxor( sourceBits, decodedBits(1:numSourceBits)));
        numberCodedBitErrors(snrdBcounter)=numberCodedBitErrors(snrdBcounter) + ...
                                       numCodedBitErrorsInFrame;
        %=================
        

        channelErrorPattern=bitxor( interleavedCodedBits, detectedChannelBits);
        numChannelBitErrorsInFrame=sum(channelErrorPattern);
        numberChannelBitErrors(snrdBcounter)=numberChannelBitErrors(snrdBcounter) + ...
                                        numChannelBitErrorsInFrame;
        if numChannelBitErrorsInFrame > 0
            numberOfChannelFrameErrors(snrdBcounter)= ...
                                numberOfChannelFrameErrors(snrdBcounter)+1;
        end

        if numCodedBitErrorsInFrame > 0
            numberOfCodedFrameErrors(snrdBcounter)= ...
                                numberOfCodedFrameErrors(snrdBcounter)+1;
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end    % for snrdBcounter = 1 : length(snrArray)
end    % for frameNumber = 1 : maxSimulationFrames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  toc

 %profile viewer

channel_BER= numberChannelBitErrors ./ numSimulatedBits;
coded_BER= numberCodedBitErrors ./ numSimulatedSourceBits;
coded_FER=numberOfCodedFrameErrors ./ numSimulatedFrames;
%=========================================
errorArray(:,1)=snrArray ;
errorArray(:,2)=channel_BER ;
errorArray(:,3)=coded_BER ;
errorArray(:,4)=coded_FER ;
%========================================

figure;
semilogy( ...
    uncodedEb_NodBarray(1,:) , channel_BER(1,:),'k-o' ,  ... 
    codedEb_NodBarray(1,:) , coded_BER(1,:),   'k-d' , ... 
    codedEb_NodBarray(1,:) , coded_FER(1,:),   'k-s'  ... 
)

ha=gca ;
set(ha,'fontsize',16)
grid
ylim([1e-3 1])
xlabel('E_b/N_o in dB')
ylabel('Bit and Frame Error Rate')
titleText= ['Error Rate for ', TXTmodType, ' in ', TXTchannelModel] ;
title(titleText)
legend('Channel BER' , 'Coded BER' , 'Coded FER')

ha1=text(1, 1e-1, [TXTchannelModel,',' blanks(1) , ...
                         TXTmodType,',' ' FFT size=' num2str(fftSize)]);
set(ha1,'FontSize',14)
set(ha1,'background','white')

ha2=text(1, 1e-2, TXTpermutation);
set(ha2,'FontSize',14)
set(ha2,'background','white')

fileName=['results\', TXTchannelModel, '_', TXTmodType, '_ofdmSymbols', ...
          num2str(numOFDMsymbols), '_', ...
          num2str(maxDelaySpreadInSamples), 'samplesDelay', '_fft' num2str(fftSize), ...
          '_Prefix' num2str(cyclicPrefixRatio)];  
fileNameExt=[fileName, '.dat'];
fileNameExtfig = [fileName, '.fig'];
%saveas(gcf, fileName, 'fig')
%saveas(gcf, fileName, 'jpg')
savefig(gcf, fileNameExtfig)
dlmwrite(fileNameExt, errorArray, 'delimiter','\t')


%=======================================================

