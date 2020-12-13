clc
close all
clear all
%%
qam = 16; 
SNR = 30; % SNR at output of channel
equalizer_length = 35;
Delta =15; %approx. delay by channel + half length of equalizer
mu = 0.4;
epsilon = 1e-6; 

training = 500;  % length of training with QPSK data
decision_directed = 5000; % length of decision-directed with QAM data
N=training + decision_directed; % total number of symbols
train_data = zeros(1,training+Delta); %training data
s=zeros(1,N);

%Channel
channel=[0.5 1.2 1.5 -1];

%Noise variance during decision-directed
%QAM symbols don't have unit-variance
sigma_s = sqrt(10);
sigma_v_dd = sqrt((sigma_s)^2*norm(channel)^2/10^(SNR/10)); % during decision-directed
sigma_v_tt = sqrt(norm(channel)^2/10^(SNR/10)); % during training

v=zeros(1,N);
v(1:training)=(sigma_v_tt^2/2)*(randn(1,training)+1i*randn(1,training));
v(training:N)=(sigma_v_dd^2/2)*(randn(1,N-training+1)+1i*randn(1,N-training+1));


%Data
s(1:training)=(sign(randn(1,training))+1i*sign(randn(1,training)))/(sqrt(2)); 
train_data(Delta+1:training+Delta)=s(1:training); % QPSK training data
   % generates QPSK data; it has unit variance
   % The training vector starts with Delta zeros.
   

for i=1:decision_directed
 xint=randi([0,(qam-1)],1);
 xcoor=qammod(xint,qam);% gives real and imaginary coordinates
 s(training+i)=xcoor; 
end 
 
y=filter(channel,1,s); % filters data through channel
r=y+v;%awgn(y,SNR,'measured');

% Equalizer initialization
w  = zeros(equalizer_length,1); % equalizer coefficients (column) 
u  = zeros(1,equalizer_length); 
e = zeros(1,N); % error vector
num_errors=0;

%  Adaptive Equalization

% Training mode
for i = 1:training+Delta
   u  = [r(i) u(1:equalizer_length-1)];
   hat_s(i) = u*w;
   d(i) = train_data(i); % training
   e(i) = d(i)- hat_s(i);
   w = w + mu*e(i)*u'/(norm(u)^2+epsilon);
end

% Decision-directed mode
for i = training+Delta+1:N
   u  = [r(i) u(1:equalizer_length-1)];
   hat_s(i) = u*w;
   check_s(i)=slicer16(hat_s(i));  
   d(i) = check_s(i); % decision_directed
   e(i) = d(i)- hat_s(i);
   w = w + mu*e(i)*u'/(norm(u)^2+epsilon);
   if (check_s(i) ~= s(i-Delta))
     num_errors = num_errors + 1;
   end
end  
figure

subplot(221); 
plot(real(train_data(Delta+1:training)),imag(train_data(Delta+1:training)),'.');
grid; 
title('Training sequence'); 

subplot(222); 
plot(real(s(training+Delta+1:N)),imag(s(training+Delta+1:N)),'.');
grid; 
title('Transmitted sequence'); 
 
subplot(223); 
plot(real(r(training+Delta+1:N)),imag(r(training+Delta+1:N)),'.');grid; 
title('Received sequence'); 
 
subplot(224); 
plot(real(hat_s(training+Delta+1:N)),imag(hat_s(training+Delta+1:N)),'.');grid; 
title('Equalizer output');