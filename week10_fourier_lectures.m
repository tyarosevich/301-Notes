%% Fourier lectures

clc;clear all;close all;
% THis is the discrete or fast fourier which pertains to data, as opposed
% to continuous functions.

% Some functions are well approximated by Tyalor Series expansions, like
% x^4 - x^2 (if you examine the polynomial of deltax that emerges)

% Some are ill suited however.
% For example sine waves (because they are hard to represent with
% polynomials. Every time you want to capture another period, you need
% another order in the taylor expansion, interesting.
% 
% So instead of expanding a function in terms of polynomials, it is
% possible to represent our function as a sum of sine waves. This is the
% basic idea of a fourier series.
% 
% There is an entire theory of expanding functions in terms of sines, or
% sin(k*omega*x) and this is called fourier analysis.
% So we are intersted in extracting frequencies from data. 
% Suppose we have some data points and we want to extract some frequency
% components at each point. We would use a discrete fourier transform to do
% this. So at each point f(x_n) we would have fhat(x_n);

% The DFT is what converts all the f points into the frequency extractions,
% or fourier coefficients. 
% The DFT is:

%  fhat_k = sum n=0 to N-1  f_n * e^ (-2 * pi * i * K * n / N)

% recal Euler's formula: e^ (i * theta) = cos(theta) + i * sin(theta)
% Think of this as e to the i*theta is a coordinate in the complex plane,
% ON the complex unit circle. So your real part is cos(theta) in the real direction and the
% imaginary part is sin(theta) or the imaginary direction
% Returning to fourier, 2 pi is a full revolution of that circle in the
% imaginary plane. The K*n/N tells us what part of the full revolution in
% question. So 2piKn/N is some frequency w, and the fhats are the
% coefficients of the different frequencies, or the magnitudes of each sin
% wave. 

% This is reversable. 
% f_n = 1/N sum k=0 to N-1 fhat * e ^ (2*pi*i*k*n/N)

% So in the example we have f_0_hat * some 0 frequency function, this is
% the sort of avergae slope of the whole thing.
% Then you have the second coefficient making some piece of a sine wave,
% etc etc, and each one is higher freq than the last.


% So how to compute? consider w_n = e ^ -2pi * i/N
% N determines what portion of the circle we are revolving around.
% So this all takes the form Ax = b where b are our frequency
% coeffiecients.
% Just look it up. It's the discrete fourier transform matrix. Also called
% the Vandermonde matrix.
% Note that the fhat coefficients tell us how much of the signal is given
% by each frequency band.
% This is an O(n^2) operation.
% Next lecture we'll learn the fast fourier transform, which is O(n*log(n))
% operation.
% FFT and DFT have become synonymous.
%% Simple signal with 2 freq example
dt = .001;
t = 0:dt:1;
x = sin(2*pi*50*t) + sin(2*pi*120*t); %w1 = 50hz, w2 = 120hz
figure(1)
plot(t,x,'LineWidth',1.2)

figure
N = length(t);
Y = fft(x, N); %fft is a fast discrete fourier

% what we're really intersted in is the magnitude of each frequency, which
% tells us how much of the freq is in each fhat. This is called Power
% Frequency D something
PSD = Y.*conj(Y)/N;
freq = 1/(dt*N)*(0:N); %this creates the x axis of frequencies in Hz
L = 1:floor(N/2); % only plot the first half of freqs)
plot(freq(L), PSD(L))
title('Power Spectrum')
xlabel('Frequency(Hz)')

% Note all the power is at 50 hz and 120 hz, which is as expected. 

%% Second lecture on DFT

% Example: IO of audio:

% Audio is sampled at 44.1 khz, or 44,100 samples per second
% So 10 seconds of audio is about 441,000 samples, so that's N
% So n^2 = 1,944,800,000,000,000

% Or even 1 second n^2 is 2 billion about, and n * log(n) is 1/4 million.
% 8000x speedup
% FFT was invented by Cooley and Tuky around 1965 and Cooley was at IBM and
% Tuky was at Princeton. HOWEVER it was actually discovered 160 years
% earlier by Gauss in 1805. 

%% FFT 

% SO it turns out that DFT can be implemented much more efficiently if my N
% is a power of 2.

% So we're gonna say, for example, F_1024 is a matrix taking my data x to
% my fourier coefficients xhat. 

% So at the core of this is that if xhat = F * x then we can break this
% down into :
% xhat = [I -D; I -D] [ F_512 0;0 F_512] [xeven; xodd]
% The upshot is that the 512 matrices are much easier to compute, and we
% recapture the correct values with the first matrix. D is a diagonal
% matrix of w values basically. 
% AHhhh and then the F512s get broken up in to 256, etc etc. THen no doubt
% the algorithm is recursive. 

% Now even if N does not = a power of 2, you just pad the the data set with
% zeros until it's a power of 2. 
%% Now filtering noise
clc;clear all;close all;
dt = .001;
t = 0:dt:1;
x = sin(2*pi*50*t) + sin(2*pi*120*t); %w1 = 50hz, w2 = 120hz
figure(1)
subplot(2,1,1)
plot(t,x,'LineWidth',1.2)
% add some noise
y = x + 2.5 * randn(size(t));
hold on
plot(t,y,'r','LineWidth', 1.2)
axis([0 .25 -5 5])
legend('Clean', 'Noisy')


figure
N = length(t);
Y = fft(y, N); %fft is a fast discrete fourier

% what we're really intersted in is the magnitude of each frequency, which
% tells us how much of the freq is in each fhat. This is called Power
% Frequency D something
PSD = Y.*conj(Y)/N;
freq = 1/(dt*N)*(0:N); %this creates the x axis of frequencies in Hz
L = 1:floor(N/2); % only plot the first half of freqs)
plot(freq(L), PSD(L))
title('Power Spectrum')
xlabel('Frequency(Hz)')

% Note all the power is at 50 hz and 120 hz, which is as expected.
%% Now we're going to filter out noise
indices = PSD > 50; %this will create vector of booleans for indices of PSD when greater than 50
PSD = PSD.*indices; %this of courese zeros out powers below 50, or 'removes the noise floor'
hold on
plot(freq(L), PSD(L), 'r');
legend('original','filtered')

Y = Y.*indices; %zero out the small fourier coefficients in Y
yfilt = ifft(Y); % inverse FFT to get filtered time domain signal.
figure(1)
subplot(2,1,2)
plot(t,x,'b', 'LineWidth', 1.2)
hold on
plot(t,yfilt,'r')
legend('Clean', 'Filtered')
