function [outputBits ] = interleaverDeinterleaver(inputBits, ...
                                     inter_0_deinter_1, state)

%===============================================
if inter_0_deinter_1 == 0 
    outputBits = randintrlv(inputBits, state);
end

if inter_0_deinter_1 == 1 
    outputBits = randdeintrlv(inputBits, state);
end
%===============================================
