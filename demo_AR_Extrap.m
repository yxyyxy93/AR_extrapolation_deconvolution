%  from Jeroen Vandendriessche:
% modified by yxy

clc;
clear;
close all;
fclose all;

% gaussian pulse
tc = gauspuls('cutoff', 5e6, 0.6, [], -60);
t = (-tc : 1e-8 : tc);
yi = gauspuls(t, 5e6, 0.6);
plot(yi);
title("The intial pulse");

% get the sampling frequency from the pulse
fs = length(t) / (2 * tc); 

% get the spectrum of the pulse
n = 2^nextpow2(length(yi));
Y = fft(yi, n);
f = fs*(0:(n/2))/n;
P = abs(Y(1:n/2+1)/n);
% this function is used to select the bandwidth.
% -3db, -6db, or .....
[f1, f2] = fx_calculate_bw(f, P, -3);

% reflection amplitudes
refAmps = [zeros(99,1);1;...
    zeros(199,1);...
    0.5;zeros(199,1);...
    -1;zeros(199,1);...
    0.5;zeros(49,1);...
    0.5;zeros(199,1);...
    0.5;zeros(24,1);...
    0.5;zeros(199,1);...
    0.5;zeros(12,1);...
    0.5;zeros(100,1)];

figure, plot(refAmps);
title("The reference spikes");

% test signal without noise
testSignal = conv(refAmps, yi);
figure, plot(testSignal);
title("The reference signals");

% add various levels of noise
noisedb = flip([0, 10, 20, 30, 40, 50]);
Signal(:,1) = testSignal;

for i = 1:length(noisedb)
    Signal(:, i+1) = awgn(testSignal, noisedb(i), 10 * log10(rms(testSignal)^2));
%     figure, plot(Signal);
end

% choose a Signal here
deconv = class_deconv(Signal(:, 3), yi', fs);

% deconvolved_s = deconv.deconvolution ;
% figure, plot(deconvolved_s);
% wiener_deconv_s = deconv.wiener_deconlution(1e-2);
% figure, plot(wiener_deconv_s);

% you could choose Q factor, model name('bg' or 'yk'), and order of the AR model here.
wienerdeconv_ARextrap_s = deconv.wienerdeconv_ARextrap(1e-8, [f1 f2], 10, 'yk');
figure, plot(real(wienerdeconv_ARextrap_s));
title("The AR_expolated signal");

