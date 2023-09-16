function [TDOA, peak] = rotem_TDOA_LMS(signal1, signal2, M, mu) 
    % We are supposed to make use of TK operator to improve SNR and remove
    % LF interferences, signals were passed through the operator before
    % being passed to the function.

    % Normalizing max signal amplitudes to +1 
    x1=signal1/max(signal1); 
    x2=signal2/max(signal2); 

    %allocate
    x1_ = zeros(M,1); x2_ = zeros(M,1);
    u = zeros(2*M, 1); %[u(2,M)]' was re arranged as a row vector for computations
    u(M/2) = 1; %initialize the algorithm around M/2 
    e = zeros(1, length(x1));TDOA = zeros(1,length(x1));peak = zeros(1,length(x1));

    %compute the vector u[h2; -h1], TDOA and find peaks in propogation path
    %h1

    for n=1:length(x1)
        x1_ = [x1(n);x1_(1:length(x1_)-1)]; x2_ = [x2(n);x2_(1:length(x2_)-1)]; %reverse sequence x[n] so we get x = [x(n) x(n-1) ...]; (apply eq. 6.18 in the book)
        x = [x1_;x2_]; %stack vector x(n)

        %apply eq. 6.26
        e(n) = u'*x;
        u = u-mu*e(n)*x;
        u(M/2) = 1; %forcing g2 to an impulse response at M/2
        u = u/norm(u); %forcing ||u|| to 1 
 
        %find peaks
        [peak(n),ind] = min(u(M+1:end)); %find peak in main propogation path h1;
        peak(n) = -peak(n); %main propogation path is negetive, so negate the peak.
        TDOA(n) = ind-M/2; %algorithm was initiated around M/2, compansate for the idx shift 
    end
end 