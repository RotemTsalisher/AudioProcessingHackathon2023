% Contemporary Communication Systems Project %
%  ---Contemporary Communication Systems---  %
%  -------------M.F. Mesiya----------------  %
%  --------------Chapter 8-----------------  %
%  ------Solution by Rotem Tsalisher-------  %
%  ----Lecturer: Dr. Bnjamin Gur Salomon---  %
%% q8.23

i = 1; %figure idx;
% the signal x(t) = sin(20*pi*t) + 0.25*sin(10*pi*t) is digitized using
% Delta Modulation

% (a) the signal is sampled at fs = 1.28 kHz, step size = 1/15;
f1 = 10; f2 = 5; A1 = 1; A2 = 0.25; fs = 1280; T = 1; Ts = 1/fs; %calculating parameters for sampled sig
t = 0:Ts:T; N = length(t); %time axis; N= num of samples;
xt = A1*sin(2*pi*f1.*t) + A2*sin(2*pi*f2.*t); %generate sig
%xt = sin(2*pi*t); fs = 128; %FOR Q 8.25 !!

step = 1/15; %given value for step size
xq = zeros(1,N); y = zeros(1,N); %allocate for better performance;
xq(1) = 0; y(1) = 1; %initial conditions;
for k = 2:N
    w1 = xt(k) - xq(k-1); %diff between original sig sample and previous xq sample
    y(k) = sign(w1); %if x(t) is larger - take a step upwards. else: downwards
    xq(k) = xq(k-1) + y(k)*step;
end %apply DM computations as presented in chapter 8 in the book
figure(i); i=i+1; stem(y(1:30)); grid on; ylim([-1.1, 1.1]); title('y[n] - Delta-Modulation representation of x[n]');

% (b) recornstruct xhat from xq by passing through a LP filter
fc = 50; %chose a value that is a little over 2*wmax normalized by a factor of 1/pi;
N_ord = 6; Wn = 2*fc/fs; % compute parameters for filter;
[b,a] = butter(N_ord,Wn); % no demands on type of filter, Rs, Rp and order (chose N= 8 randomly);
% chose butterworth because the function doesn't require those args;

xhat = filter(b,a,xq); %reconstruction filter;
figure(i); subplot(211); i=i+1; plot(t,xt,'b--',t, xhat, 'black');grid on; legend('Original Sig x(t)', 'Reconstructed Sig xq'); xlabel('Time [sec]'); ylabel('x(t), xq');title('x(t) digitized using DM, {\Delta} = 1/15');

% (c) same procedure, with smaller step size step = 1/20;

step_b = 1/20;
xq_b = zeros(1,N); y_b = zeros(1,N); xq_b(1) = 0; y_b(1) = 1;

for k = 2:N
    w1_b = xt(k) - xq_b(k-1);
    y_b(k) = sign(w1_b);
    xq_b(k) = xq_b(k-1) + y_b(k)*step_b;
end 

xhat_b = filter(b, a, xq_b); %used the same filter calculated in (b);
subplot(212); plot(t,xt,'b--',t, xhat_b, 'red');grid on; legend('Original Sig x(t)', 'Reconstructed Sig xq_c'); xlabel('Time [sec]'); ylabel('x(t), xq');title('x(t) digitized using DM,  {\Delta} = 1/20');

% (c) calculate MSE for both digitized signals
mse_a = mean((xt - xhat).^2); % mse for digitized signal with step size = 1/15;
mse_b = mean((xt - xhat_b).^2); % mse for digitized signal with step size = 1/20;
comment = 'We can see that a smaller step size gave us better aprox with less noise. A smaller step size also means it will be more difficult to follow sudden big rises \ drops in the signal';
disp(['mse for step size = 1/15: ' num2str(mse_a)]);
disp(['mse for step size = 1/20: ' num2str(mse_b)]);
disp(comment);

% (d) comparing DM mse to uniform quantization mse
% from previous question; NOTE: (!!) PLEASE MAKE SURE TO RUN FILE q8_21 BEFORE RUNNING THIS ONE;
load('uniform_mse.mat'); 
disp('===================='); disp('Comparing DM mse to uniform quantization MSE');
disp(['MSE for uniform quantization case: ' num2str(MSE)]); disp(['MSE for DM with step size = 1/15: ' num2str(mse_a)]); disp(['MSE for DM with step size = 1/20: ' num2str(mse_b)]);
disp('We got a considerably larger MSE using DM');disp('====================');

i = 1; %reset;