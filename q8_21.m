% Contemporary Communication Systems Project %
%  ---Contemporary Communication Systems---  %
%  -------------M.F. Mesiya----------------  %
%  --------------Chapter 8-----------------  %
%  ------Solution by Rotem Tsalisher-------  %
%  ----Lecturer: Dr. Bnjamin Gur Salomon---  %
%% q8.21

% (a) generate a Laplacian-distributed message;
% code is given in the example in the textbook: 
global i; %global idx;
i = 1;
n_ = [8 4 6 10 12];
mu = 0; % mean of laplacian = 0;
sigma = 1; % varience of laplacian = 1; ALSO REPRESENTS THE SIGNAL POWER IN THIS SCENARIO
N = 10000;
u = rand(N,1); %random vector with N = 10000 samples;
idx1 = find(u<0.5); idx2 = find(u>0.5);
x(idx1) = log(2*u(idx1));
x(idx2) = -log(2*(1-u(idx2)));
x = x*sigma + mu; %in this case x = x

magmax = max(abs(x)); % max value of absolute values of signal
xmin = -magmax; xmax = magmax; %calculate xmax and xmin for the uniform midrise quantizer;
SQNR = []; % vector to store future computations

for ele = n_
    SQNR = [SQNR my_SQNR(ele,x,xmin,xmax,mu,sigma)]; %apply function to all elements in n_ (representing different number of bits) (f)
end
plot(n_,SQNR,'b:o'); title('SQNR as a function of number of bits'); xlabel('number of bits'); ylabel('SQNR(db)'); ylim([0 62.5]); grid on;
for i = 1:length(SQNR)
    yline(SQNR(i), ':', [num2str(SQNR(i)) '[dB]']);
end
%figure(5); plot(x);
function SQNR = my_SQNR(n,x, xmin, xmax, mu, sigma)
    global i; %import global idx
    % (b) quantize x[n] using an n bit midrise uniform quantizer
    % code is given in the example in the textbook:
    L = 2^n; %resolution;
    delta = (xmax-xmin)/L; %quantization step
    y = floor((x-xmin)/delta)*delta + delta/2 + xmin; %equation of the uniform quantizer;

    % (c) Quantization error; histogram for the error func
    e = y-x; %quantization error is give by e[n] = y[n] - x[n];
    figure(i); i=i+1; hist(e,20); grid on; title('Histogram of error func e[n] = y[n] - x[n]'); xlabel('Error'); ylabel('Amount of occurences'); %20 bit histogram of the error

    % (d) Calculate MSE
    MSE = mean(e.^2); %MSE equation for laplacian with mue = 0;
    disp(['MSE for ' num2str(n) 'bit quantizer is: ' num2str(MSE)]);
    if n == 8
        save('uniform_mse.mat','MSE');
    end

    % (e) Calculate SQNR
    SQNR = 10*log10(sigma^2/MSE); %SQNR equation; sigma.^2 is the signal power; MSE is the quantization noise power
end
