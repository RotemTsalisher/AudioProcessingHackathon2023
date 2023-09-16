% Contemporary Communication Systems Project %
%  ---Contemporary Communication Systems---  %
%  -------------M.F. Mesiya----------------  %
%  --------------Chapter 8-----------------  %
%  ------Solution by Rotem Tsalisher-------  %
%  ----Lecturer: Dr. Bnjamin Gur Salomon---  %
%% q8.25
i =1; %figure idx;
% (a) given over-sampeling rate fs = 128; generate signal x(t) =
% sin(2*pi*t); 
fs = 128; T = 1; Ts = 1/fs; t = 0:Ts:1; N = length(t); xt = sin(2*pi.*t); %calculate parameters and generate sig;

%sigma delta modulation representation, code is given in the problem:
wo = 0; y = zeros(1,N+1); y(1) = 1; yn = zeros(1,N+1); 
for k=2:1:N+1
    w1 = xt(k-1) - y(k-1) + wo;
    y(k) = sign(w1);
    wo = w1;
end
yn = y(2:N+1);

%plot the sigma-delta modulation representation yn[n];
subplot(211); stem(yn); grid on; ylabel('yn[n]'); title('Sigma-Delta modulation representation signal');

% (b) recover xhat by passing yn in a LP reconstructing filter
fc = 1.5; Wn = fc/(fs/2); N_ord = 31; b = fir1(N_ord,Wn); %filter parameters;
xhat = filter(b,1,yn); %reconstruct sig;
N_hat = length(xhat); t_hat = [0:N_hat-1]*Ts;subplot(212); plot(t,xt, 'r--'); grid on; hold on; plot(t_hat,xhat, 'black');
legend('Original sig x(t)', 'Recoverd sig xhat(t)'); xlabel('Time[sec]'); ylabel('x(t), xhat(t)'); title('Original signal x(t), Recovered signal xhat(t)'); hold off;
% comment: It seems like the oversampling system is a good fit, but there are a few
% places (like the drop around t ~ 0.22-0.26) that the oversampling method
% got us to wrongly aprox the original signal.

% (c) calculate mse and compare to DM performances
mse = mean((xt-xhat).^2); disp(['mse for oversampeling system: ' num2str(mse)]);

% comparison: we can see that in previous system, the mse was much bigger,
% and the current system produces a much lower error in terms of MSE.
