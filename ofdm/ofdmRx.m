function [ output ] = ofdmRx(input, numGuardLeft, numGuardRight, fftSize, guardRatio )

% "input" must be one row for all OFDM Symbols in time domain

% "output" is a matrix with "used subcarriers"
%  one OFDM symbol in each column

% Note that gaurdRatio * fftSize must be an integer (of course)

%=============================================
%%%%%%%%%% IMPORTANT ASSUMTION %%%%%%%%%%%%%%%%
%=============================================
% The first sample is the first sample of the first gaurd interval
% We will start deleting the guard from the first sample
% Synchronization must take care of this
%=============================================


% This is the size of the data at each side of the OFDM symbol
halfDataSize=(fftSize - 1 - numGuardLeft - numGuardRight)/2;
guardInterval=guardRatio * fftSize ;

inputLength=length(input);

% The number of OFDM symbols must give an integer
numOFDMsymbols=floor(inputLength/(fftSize+guardInterval));

% remove additional samples due to fading
UsefulInputLength=numOFDMsymbols * (fftSize+guardInterval);
input=input(1:UsefulInputLength);

% Let's reshape the input such that each column includes one OFDM symbol
% which includes its gaurd interval
x=reshape(input, (fftSize+guardInterval) , numOFDMsymbols) ;

% Now delete the guard interval, i.e., the first "guardInterval" rows
y=x(guardInterval+1:fftSize+guardInterval , :) ;

% Now apply FFT to each column
z=fft(y,fftSize);

% Now apply fftshift to return zero to middle and guard subcarriers
z=fftshift(z,1);

% We now delete the guards from top and bottom and also the zero at middle
outputTop=z(numGuardLeft+1:numGuardLeft+halfDataSize , :);
outputBot=z(numGuardLeft+halfDataSize+2 : ...
            numGuardLeft+halfDataSize+halfDataSize+1 , :);
        
output=vertcat(outputTop, outputBot);
