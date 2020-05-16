# AR_extrapolation_deconvolution
This is a demo of deconvoluton, wiener deconvolution, and Autoregressive spectrum extrapolation of 1D signal;

1. run 'demo_AR_Extrap'.

choose bw in 'demo_AR_Extrap'
% -3db, -6db, or .....
[f1, f2] = fx_calculate_bw(f, P, -3);

define signal and ref signal of the object in 'demo_AR_Extrap'
% choose a Signal here
deconv = class_deconv(Signal(:, 3), yi', fs);

choose order, method, and Q factor for AR extrapolation and wiener deconvolution in 'demo_AR_Extrap'
% you could choose Q factor, model name('bg' or 'yk'), and order of the AR model here.
wienerdeconv_ARextrap_s = deconv.wienerdeconv_ARextrap(1e-8, [f1 f2], 10, 'yk');