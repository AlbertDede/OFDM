function channelGains = channelEstimation( actualGain, numGuardLeft, numGuardRight, fftSize )


halfDataSize=(fftSize - 1 - numGuardLeft - numGuardRight)/2;
% ideal channel gains for each used subcarrier
x=fft(actualGain, fftSize);
y=fftshift(x);
%
xLeft=y(numGuardLeft+1:numGuardLeft+halfDataSize);
xRight=y(numGuardLeft+halfDataSize+2 : ...
            numGuardLeft+halfDataSize+halfDataSize+1);

channelGains=[xLeft , xRight];