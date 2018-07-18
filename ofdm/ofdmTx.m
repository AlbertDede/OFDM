function [ output ] = ofdmTx(input, numGuardLeft, numGuardRight, fftSize, guardRatio )

% "input" must be in columns. Each column is one OFDM Symbol in frequency
% domain, i.e., before IFFT

% gaurdRatio * fftSize must be an integer (of course)

% "output" is an array of samples ready to go to the channel


% This is the size of the data at each side of the OFDM
% It must be equal to half the size of the input array
halfDataSize=(fftSize - 1 - numGuardLeft - numGuardRight)/2;
guardInterval=guardRatio * fftSize ;

inputSize=size(input);

%if input is 1 row array, convert it to column
if inputSize(1) == 1
    input=transpose(input);
    inputSize=size(input);
end
numOFDMsymbols=inputSize(2);

% x below is matrix. Each column is one frequency domain OFDM symbol after
% adding the left and right guard subcarriers and the middle zero subcarrier 

x=zeros(fftSize, inputSize(2));

x(numGuardLeft+1:numGuardLeft+halfDataSize , 1:numOFDMsymbols) = ...
    input(1:halfDataSize, 1:numOFDMsymbols);

x(numGuardLeft+halfDataSize+2:numGuardLeft+ 2*halfDataSize+1, 1:numOFDMsymbols) = ...
    input(halfDataSize+1:2*halfDataSize, 1:numOFDMsymbols);

%Prepare to apply to IFFT by doing fftshift on each column
% Now the DC should go to the top element
y=fftshift(x, 1);

% apply IFFT
% Note: fft is done on columns
z=ifft(y, fftSize);

% Put the guard band on the top

zWithGuard=vertcat( z(fftSize-guardInterval+1:fftSize, :), z);

output=reshape(zWithGuard, 1, (fftSize+guardInterval)* numOFDMsymbols) ;

