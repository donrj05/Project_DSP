sine = dsp.SineWave('Frequency',375,'SampleRate',8000,'SamplesPerFrame',1000)
s = sine();
v = 0.8*randn(sine.SamplesPerFrame,1); % Random noise part.
ar = [1,1/2];          % Autoregression coefficients.
ARfilt = dsp.IIRFilter('Numerator',1,'Denominator',ar)
v1 = ARfilt(v);
x = s + v1;
ma = [1, -0.8, 0.4, -0.2];
MAfilt = dsp.FIRFilter('Numerator',ma)
v2 = MAfilt(v);
L = 7;
lms = dsp.LMSFilter(L,'Method','LMS')
nlms = dsp.LMSFilter(L,'Method','Normalized LMS')
[mumaxlms,mumaxmselms] = maxstep(lms,x)
[mumaxnlms,mumaxmsenlms] = maxstep(nlms,x)
lms.StepSize  = mumaxmselms/30 
nlms.StepSize = mumaxmsenlms/20 
[~,elms,wlms] = lms(v2,x);
[~,enlms,wnlms] = nlms(v2,x);
reset(MAfilt);
bw = firwiener(L-1,v2,x); % Optimal FIR Wiener filter
MAfilt = dsp.FIRFilter('Numerator',bw)
yw = MAfilt(v2); % Estimate of x using Wiener filter
ew = x - yw; % Estimate of actual sinusoid
n = (1:1000)';
subplot(2,2,1);
plot(n(900:end),[ew(900:end), elms(900:end),enlms(900:end)])
legend('Wiener filter denoised sinusoid',...
    'LMS denoised sinusoid','NLMS denoised sinusoid')
xlabel('Time index (n)')
ylabel('Amplitude')
hold on
plot(n(900:end),x(900:end),'k:')
xlabel('Time index (n)')
ylabel('Amplitude')
hold off
[bw.' wlms wnlms]
[ylms,elms,wlms] = lms(v2,x);
[ynlms,enlms,wnlms] = nlms(v2,x);
reset(ARfilt)
reset(sine);
release(sine);
n = (1:5000)';
sine.SamplesPerFrame = 5000
s = sine();
nr = 25;
v = 0.8*randn(sine.SamplesPerFrame,nr);
ARfilt = dsp.IIRFilter('Numerator',1,'Denominator',ar)
v1 = ARfilt(v);
x = repmat(s,1,nr) + v1;
reset(MAfilt);
MAfilt = dsp.FIRFilter('Numerator',ma)
v2 = MAfilt(v);
reset(lms);
reset(nlms);
M = 10; % Decimation factor
mselms = msesim(lms,v2,x,M);
msenlms = msesim(nlms,v2,x,M);
subplot(2,2,2)
plot(1:M:n(end),mselms,'b',1:M:n(end),msenlms,'g')
legend('LMS learning curve','NLMS learning curve')
xlabel('Time index (n)')
ylabel('MSE')
reset(lms);
[mmselms,emselms,meanwlms,pmselms] = msepred(lms,v2,x,M);
x = 1:M:n(end);
y1 = mmselms*ones(500,1);
y2 = emselms*ones(500,1);
y3 = pmselms;
y4 = mselms;
plot(x,y1,'m',x,y2,'b',x,y3,'k',x,y4,'g')
legend('MMSE','EMSE','Predicted LMS learning curve',...
    'LMS learning curve')
xlabel('Time index (n)')
ylabel('MSE')