% Contemporary Communication Systems Project %
%  ---Contemporary Communication Systems---  %
%  -------------M.F. Mesiya----------------  %
%  --------------Chapter 8-----------------  %
%  ------Solution by Rotem Tsalisher-------  %
%  ----Lecturer: Dr. Bnjamin Gur Salomon---  %
%% q8.22

%(a) Generating a laplacian as in q8.20, and applying mu=255 law
%companding

u = rand(10000,1); %generating the signal, code is given in q8.20
idx1 = find(u<0.5); idx2 = find(u>0.5);
x(idx1) = log(2*u(idx1));
x(idx2) = -log(2*(1-u(idx2)));

mu = 255; magmax = max(abs(x)); xmin = -magmax; xmax = magmax; %computing parameters for the mu law companding
y = xmax*log10(1+abs(x)*(mu/xmax))/log10(1+mu); %equation of the mu law companding, as defined in chapter 8 in the book, applied to the sequence x[n];

% (b) using 8 bit uniform quantization to acquire yq[n]:
% with no other demand, assuming a midrise quantizer

n=8; % bit depth
L = 2^n; % resolution;4
delta = (xmax-xmin)/(L-1); % quantization step

yq = floor((y)/delta)*delta + delta/2; %apply uniform quantiztion

% (c) apply inverse mu = 255 law mapping to recover uncompressed sig:
xq = (xmax/mu)*(10.^((log10(1+mu)/xmax)*yq)-1).*sign(x); % inverse mue law equation as given in the book (and in the question);
mse_nonuni = mean((x - xq).^2); %calculate mse
disp(['MSE for non-uniform quantization case:' num2str(mse_nonuni)]);

% (d) comparing non uniform quantization mse to uniform quantization mse
% from previous question; NOTE: (!!) PLEASE MAKE SURE TO RUN FILE q8_21 BEFORE RUNNING THIS ONE;
load('uniform_mse.mat'); 
disp(['MSE for uniform quantization case:' num2str(MSE)]);

if mse_nonuni<MSE %conclution (you can run q8.21 and q8.22 a couple of times to make sure all cases work properly);
    disp('As we can see, we got better performance in terms of MSE from the non uniform quantization');
else
    disp('As we can see, we got better performance in terms of MSE from the uniform quantization');
end