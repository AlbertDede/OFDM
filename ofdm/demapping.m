function [ softOutput ] = demapping( inputSymbols, modulationIndex, gain2 )

softOutput=zeros(1, length(inputSymbols) * modulationIndex);

% At this time there the same channel gain for all OFDM symbols
numOFDMSymbols=length(inputSymbols)/length(gain2);
gain2=repmat(gain2, 1, numOFDMSymbols);

 % modulationIndex = 1  BPSK
 % modulationIndex = 2  QPSK
 % modulationIndex = 4  16-QAM
 % modulationIndex = 6  64-QAM
 
 if modulationIndex == 1
     scaling=1;
 elseif modulationIndex == 2 
     scaling=sqrt(2);
 elseif modulationIndex == 4 
     scaling=sqrt(10);
 elseif modulationIndex == 6 
     scaling=sqrt(42);
 end
 inputSymbols=inputSymbols * scaling ;
 realPart=real(inputSymbols);
 imagPart=imag(inputSymbols);

 if modulationIndex == 1 % BPSK
     
     softOutput = realPart;
     
 elseif modulationIndex == 2 % QPSK
     
    softOutput(1:2:end)=realPart;
    softOutput(2:2:end)=imagPart;
    
 elseif modulationIndex == 4 % 16 QAM
     
   softOutput(1:4:end)= realPart ;
   softOutput(3:4:end)= imagPart ;

   slope=1; % absolute value of slope
   softOutput(2:4:end)= 2*slope*gain2 - slope*abs(realPart);
   softOutput(4:4:end)= 2*slope*gain2 - slope*abs(imagPart);
       
 elseif modulationIndex == 6
       softOutput(1:6:end)= realPart ;
       softOutput(4:6:end)= imagPart ;
       
       absRealPart=abs(realPart);
       absImagPart=abs(imagPart);
 
       slope1=1;
       softOutput(2:6:end)= 4*slope1*gain2 - slope1*absRealPart;
       softOutput(5:6:end)= 4*slope1*gain2 - slope1*absImagPart;
       
       slope2=1;
       %Work on the real part
       % generate +1 for any element > 4*gain2 and 0 otherwise
       biggerThan4=  (absRealPart > 4*gain2)  ;
       
       % generate +1 for any element >4*gain2 and -1 otherwise
       biggerThan4Bipolar= 2*(biggerThan4-0.5);
       
       mirrotAround_4 = 8*gain2 .* biggerThan4 - absRealPart .* biggerThan4Bipolar ;
       
       softOutput(3:6:end)=slope2 * mirrotAround_4 - 2*slope2*gain2 ;
       
       %Work on the imaginary part
       % generate +1 for any element >4*gain2 and 0 otherwise
       biggerThan4=  (absImagPart > 4*gain2)  ;
       
       % generate +1 for any element >4*gain2 and -1 otherwise
       biggerThan4Bipolar= 2*(biggerThan4-0.5);
       
       mirrotAround_4 = 8*gain2 .* biggerThan4 - absImagPart .* biggerThan4Bipolar ;
       
       softOutput(6:6:end)=slope2 * mirrotAround_4 - 2*slope2*gain2 ;
       
       
       
      
%           if ( abs(real(rxSymbol)) .lt. 4.0*gain2) then
%                  softOutput(i*M_index+2)= abs(real(rxSymbol)) - 2.0*gain2 
%           else
%                  softOutput(i*M_index+2)=  6.0*gain2 - abs(real(rxSymbol)) 
%           end 
          
%           if ( abs(aimag(rxSymbol)) .lt. 4.0*gain2) then
%                  softOutput(i*M_index+5)= abs(aimag(rxSymbol)) - 2.0*gain2 
%           elseif
%                  softOutput(i*M_index+5)=  6.0*gain2 - abs(aimag(rxSymbol)) 
%           end 

          
 end
