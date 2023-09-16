% Contemporary Communication Systems Project %
%  ---Contemporary Communication Systems---  %
%  -------------M.F. Mesiya----------------  %
%  --------------Chapter 8-----------------  %
%  ------Solution by Rotem Tsalisher-------  %
%  ----Lecturer: Dr. Bnjamin Gur Salomon---  %
%% q8.24

i = 1; %figure idx;

% (a) ADM of x(t) = 10sin(2*pi*t), with I.C xq(1)=0, y(1) =0, sampeling frequency = 128 ;
fs = 128; T = 1; Ts = 1/fs; t = 0:Ts:1; xt = 10*sin(2*pi.*t); %calculate parameters and define sig

% ADM parameters
N = length(t); y = zeros(1,N); xq = zeros(1,N); step = zeros(1,N); step_min = 0.125;

%apply ADM (same as DM but with the adaptive step size);
for k = 2:N
    w1 = xt(k) - xq(k-1); %diff between original sig sample and previous xq sample
    y(k) = sign(w1); %if x(t) is larger - take a step upwards. else: downwards
    step(k) = ((abs(step(k-1))) / y(k)) * (y(k) + 0.5 * y(k-1)); %make regulare computation of step size, then determine if it's too small;
    if(step(k-1)<step_min)
        step(k) = step_min;
    end

    xq(k) = xq(k-1) + y(k)*step(k);
end %apply DM computations as presented in chapter 8 in the book


% plot the ADM representation (y[n]);
figure(i); i=i+1; stem([0:length(y)-1],y); grid on; ylim([-1.05 1.05]); ylabel('y[n]'); title('y[n] - ADM representation with adaptive {\Delta}');

% (b) reconstruct xq by accumulating y[n], recunstruct xhat by passing xq
% in a LP filter

fc = 10; Wn = 2*fc/fs; N_ord = 6; [b,a] = butter(N_ord, Wn); %calculate parameters of filter and filter coeff
xhat = filter(b,a,xq); %pass xq through LP filter to get xhat;
figure(i); clf(figure(i)); i = i+1;plot(t,xt,'b--');hold on; plot(t,xhat,'black'); grid on; title('Original signal x(t) and xhat - reconstructed from xq');  ylabel('x(t), xhat'); legend('x(t)', 'xhat'); xlabel('Time [sec]');

% (c) repeat (a), (b) with step = 1/20 (same as applied in q8.23):

step_b = 1/20;
xq_b = zeros(1,N); y_b = zeros(1,N); xq_b(1) = 0; y_b(1) = 0;

for k = 2:N
    w1_b = xt(k) - xq_b(k-1);
    y_b(k) = sign(w1_b);
    xq_b(k) = xq_b(k-1) + y_b(k)*step_b;
end 

xhat_b = filter(b, a, xq_b); %used the same filter calculated in (b);

% plot the DM representation (y[n]);
figure(i); i=i+1; stem([0:length(y_b)-1],y_b); grid on; ylim([-1.05 1.05]); ylabel('y[n]'); title('y[n] - DM representation with {\Delta} = 1/20');
figure(i); i=i+1; plot(t,xt,'b--',t, xhat_b, 'red');grid on; legend('Original Sig x(t)', 'Reconstructed Sig xq_c'); xlabel('Time [sec]'); ylabel('x(t), xq');title('x(t) digitized using DM,  {\Delta} = 1/20');

% (d) calculate mse and compare it to mse of DM
mse_ADM = mean((xt-xhat).^2); disp(['mse for ADM system: ' num2str(mse_ADM)]);
mse_DM = mean((xt-xhat_b).^2); disp(['mse for DM system: ' num2str(mse_DM)]);

i=1; %reset;