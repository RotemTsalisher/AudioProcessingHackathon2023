function [res] = rotem_gcc(signal1, signal2, type)
   
   %length, zero-padding for circ-convo=linear-convo, and FFT computations
   %are relevant for both cases

   M = min(length(signal1),length(signal2)); 
   NFFT = 2*M-1; % so that linear convolution = circular convolution 
   x1 = signal1 - mean(signal1); %center signals around avg=0; (avoid DC components for transform);
   x2 = signal2 - mean(signal2); %center signals around avg=0;
   X1 = fft(x1,NFFT); 
   X2 = fft(x2,NFFT);
   S_x1x2 = X1.*conj(X2); %compute multiplication in freq

   switch type
       case 'cc'
            %Computing the standard biased cross-correlation via FFT-IFFT  
            phi_x1x2 = ifft(S_x1x2); 
            % re-arranging the IFFT 
            res = [phi_x1x2(NFFT-M+2:NFFT); phi_x1x2(1:M)]'; 
       case 'phat'
            phi_x1x2 = ifft(S_x1x2 ./ max(abs(S_x1x2),eps)); 
            % re-arranging the IFFT 
            res = [phi_x1x2(NFFT-M+2:NFFT); phi_x1x2(1:M)]';
       otherwise
            error('Unsupported correlation type. Use ''cc'' or ''phat''.');
   end
end