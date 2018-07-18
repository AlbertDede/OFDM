function [ outputSymbols ] = channelEqualizer(inputSymbols, channelGainsFrequency)

% "inputSymbols" is a matrix. Each column is one OFDM Symbol in frequency domain
% The number of columns in "inputSymbols" is the number of OFDM symbols
%
% "outputSymbols" is a one row array of all modulation symbols

% "channelGainsFrequency" is a row with channel gains in frequency domain

inputSize=size(inputSymbols);
numOFDMsymbols=inputSize(2);
numModulationSymbols=inputSize(1) * inputSize(2);

equalizerMatrix=repmat(channelGainsFrequency', 1, numOFDMsymbols);

equalizedSymbols=inputSymbols .* equalizerMatrix ;

outputSymbols=reshape(equalizedSymbols, 1, numModulationSymbols);