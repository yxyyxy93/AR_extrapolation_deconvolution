classdef class_deconv
    % used for deconvolution
    % including the Autoregressive
    
    properties
        recorded_signal
        ref_signal
        fs
    end
    
    methods
        function obj = class_deconv(recorded_signal, ref_signal, fs)
            % a class to implement the deconvolution, including the
            % AutoRregressive extrapolation
            obj.recorded_signal = recorded_signal;
            obj.ref_signal = ref_signal;
            obj.fs = fs;
        end
        
        function deconv = deconvolution(obj)
            % normal deconvolution
            % just divide the spectrum, ignoring the noises
            n = 2^nextpow2(length(obj.recorded_signal));
            deconv_spectrum = fft(obj.recorded_signal, n) ./ (fft(obj.ref_signal, n));
            deconv = ifft(deconv_spectrum);
            
        end
        
        function [deconv, deconv_spectrum] = wiener_deconlution(obj, Q_factor)
            % wiener deconvolution
            % Q_factor: the Q factor in wiener deconvolution function
            % output the deconv_spectrum for AR extrapolation
            if nargin == 1
                Q_factor = 1e-2; % should be 10^(-2) here
            end
            
            kernel = obj.ref_signal;
            ori_signal = obj.recorded_signal;
            shift = round(length(kernel) / 2);
            n = 2^nextpow2(length(obj.recorded_signal));
            H = fft(kernel, n);
            Y_omega  = fft(ori_signal, n);
            % This constant is sometimes called the ‘‘noise desensitizing factor’’
            Q = sqrt(Q_factor * max(H .* conj(H)).^2);
            deconv_spectrum = Y_omega .* conj(H) ./ (H .* conj(H) + Q.^2);
            deconv = ifft(deconv_spectrum);
            % The shift maybe not needed
            %             deconv = circshift(devolved, shift);
        end
        
        function deconv = wienerdeconv_ARextrap(obj, Q_factor, f_window, k, method)
            % Autoregressive extrapolation
            % Q_factor: Q_factor in wiener deconvolution
            % f_window: the spectrum window to fit the model;
            % k: order the the AR model
            
            % calculate the wiener deconvolved specturm
            [~, deconv_spectrum] = wiener_deconlution(obj, Q_factor);
            ns = length(deconv_spectrum);
            f_window_index = round(ns * f_window / obj.fs); % transfer the frequency window to the index.
            m = f_window_index(1);
            n = f_window_index(2);
            spectrum_window = deconv_spectrum(m: n);
            
            % debug
            figure, plot(abs(deconv_spectrum))
            hold on;
            scatter(f_window_index, abs(deconv_spectrum(f_window_index)));
            title("The deconvolved spectrum");
            
            % AR model
            % choose methods
            if method=='bg'
                [a, e, rc] = arburg(spectrum_window, k);
            elseif method=='yk'
                [a, e, rc] = aryule(spectrum_window,k);
            else
                error("please choose method as 'bg' or 'yk'. ");
            end
            % transverse a and rc
            a = a';
            a_conj = conj(a);
            rc_conj = conj(rc);

            % extrapolate the miss values outside the window
            AR_spetrum = deconv_spectrum;
            % p = 1, 2, ... m - 1
            for p = m-1:-1:1
                AR_spetrum(p) = - sum( deconv_spectrum(p + 1: p + k) .* rc_conj(1: end) );
            end
            % q = n + 1, n + 2, ... Ns
            for q = n + 1:1:ns
                AR_spetrum(q) = - sum( deconv_spectrum(q - 1: -1: q - k) .* rc(1: end) ) ;
            end
            
            % debug
            figure, plot(abs(AR_spetrum)) ;
            deconv = ifft(AR_spetrum);
            title("The AR_extrapolated spectrum");
            
        end
    end
end

