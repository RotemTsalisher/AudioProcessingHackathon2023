% Contemporary Communication Systems Project %
%  ---Contemporary Communication Systems---  %
%  -------------M.F. Mesiya----------------  %
%  --------------Chapter 8-----------------  %
%  ------Solution by Rotem Tsalisher-------  %
%  ----Lecturer: Dr. Bnjamin Gur Salomon---  %

%% q8.20

global i;
i = 1; %figure number running idx;

% (a) and (b) for fs = 2800;
%fs = 2800; q8_20_func(fs); % remove to ivew solution
% (a) and (b) for fs = 1000;
fs = 1000; q8_20_func(fs); % remove to view solution


function q8_20_func(fs) %define a function to run the given task as a func of fs
    global i; %import running figure idx;

    % (a): Generate x[n] with length of N = 4096, and plot it's spectrum;
    
    N = 4096; % length (in samples)
    Ts = 1/fs; 
    f = [500 600 700]; %freq components of sig
    A = [1 0.5 2]; % amplitudes of sig
    T = (N-1)*Ts; % signal duration
    t = 0:Ts:T; %t axis
    n = 0:N-1; %n axis

    xt =@(t) A(1)*sin(2*pi*f(1).*t) + A(2)*sin(2*pi*f(2).*t) + A(3)*sin(2*pi*f(3).*t); %generate signal x(t);
    xn = xt(n*Ts); %by sampeling theorem (x[n] = xt(t=n*Ts));

    X = (fft(xn,N)); %FFT of x[n] (N-length DFT of x);
    magX = abs(X); %compute the magnitude of the FFT
    f_ = (-N/2:N/2-1)*fs/N; %compute the frequency axis ([hz]);


    % % % % % plots % % % % %
    figure(i); i=i+1;
    
    %debugging purposes only
    %subplot(311); plot(t(1:128),xt(t(1:128))); title('x(t) = sin(2{\pi}500t) + 0.5sin(2{\pi}600t) + 2sin(2{\pi}700t)'); xlabel('t[sec]'); ylabel('x(t)'); grid on;
    
    subplot(211); stem(n(1:128),fftshift(xn(1:128))); title('x[n] = sin(2{\pi}500nTs) + 0.5sin(2{\pi}600nTs) + 2sin(2{\pi}700nTs)'); xlabel('n [discrete axis]'); ylabel('x[n]'); grid on;
    subplot(212); plot(f_,fftshift(magX)); title('Magnitude of X(f) - DTFT of x[n]'); xlabel('f[hz]'); ylabel('|X(f)|'); grid on;
    % % % % % % % % % % % % %

    % (b): 7-pole eliptic filter, as a reconstruction filter
    ny_freq = fs/2;
    fpass = 0.9*ny_freq; fstop = 0.99*ny_freq; Rp = 1; Rs = 40; % parameters: frequencies in [hz] and Passband ripple,Stopban attenuation in [dB],
    % a common practice is to set the edge of the passband to be around 0.9 ny_freq, and the stop band freq to be around 0.95~0.99 ny_freq

    Wp = fpass/ny_freq; Ws = fstop/ny_freq; % normlize parameters for ellipord function 
    % (kind of weird to set the parameters as a multiplication
    % of ny_freq, and then divide by it, but I chose to do it anyway, to keep
    % the logical flow clean and to state the logical value of each stage in
    % the proccess)

    [~, Wn] = ellipord(Wp,Ws,Rp,Rs); %calculate cutoff angular freq Wn of filter;
    N_ord = 7; %as requested
    [b,a] = ellip(N_ord,Rp,Rs,Wn,'low'); %calculate coeff of filter; b = num, a = den;


    % % % % % plots % % % % %
    % Remove '%' to view: figure(1): Magnitude and Phase of filter, figure(2):
    % Poles and Zeros of filter on Z-plane
    % figure(i); i=i+1; freqz(b,a,1024,fs);
    % figure(i); i=i+1; zplane(b,a);
    % % % % % % % % % % % % %

    yn = filter(b,a,xn); %% reconstructing yn using the filter we created
    Y = (fft(yn,N)); %FFT of y[n]
    magY = abs(Y); %compute the magnitude of the FFT
    % % % % % plots % % % % %
    figure(i); i=i+1;
    subplot(211); plot(f_, fftshift(magY)); title('Magnitude of Y(f) - DTFT of y[n]'); xlabel('f[hz]'); ylabel('|Y(f)|'); grid on; 
    subplot(212); plot(f_,fftshift(magX)); title('Magnitude of X(f) - DTFT of x[n]'); xlabel('f[hz]'); ylabel('|X(f)|'); grid on;
    % % % % % % % % % % % % %
end
