function [ outputSymbols ] = mapping( inputBits, modulationIndex )


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
 
numOutputSymbols=length(inputBits)/modulationIndex ;

 if modulationIndex == 1
     outputSymbols = 1-2*inputBits;
 elseif modulationIndex == 2
    BPSK=1-2*inputBits;
    outputSymbols=(BPSK(1:2:end) + j*BPSK(2:2:end)) ;    
 elseif modulationIndex == 4
    outputSymbols=(1-2*inputBits(1:4:end)) .* (1+2*inputBits(2:4:end)) ...
        + j*      (1-2*inputBits(3:4:end)) .* (1+2*inputBits(4:4:end)) ;
     
 elseif modulationIndex == 6
    outputSymbols=(1-2*inputBits(1:6:end)) ...  
	  			      .* (4 +(2*inputBits(2:6:end)-1) .* (1+2*inputBits(3:6:end)))...
                      + j* ...
                  (1-2*inputBits(4:6:end))   ...
	 			      .* (4 +(2*inputBits(5:6:end)-1) .* (1+2*inputBits(6:6:end)));
 end
 
 outputSymbols=outputSymbols/scaling ;
