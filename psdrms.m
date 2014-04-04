function [signals_psd, signals_integrated_rms, freq, freq_sel] = psdrms(signals, Ts, freq_low, freq_high, window, noverlap, nfft, method)

if nargin < 8
    method = 'rms';
end

% Sampling frequency
Fs = 1/Ts;

% Remove DC component
npts = size(signals,1);
signals = signals - repmat(mean(signals), npts, 1);

% Compute PSD and RMS integrated
for j=1:size(signals, 2)
    [signals_psd(:,j), freq] = pwelch(signals(:,j), window, noverlap, nfft, Fs, 'onesided');
    
    index = find((freq >= freq_low-eps) & (freq <= freq_high+eps));
    
    if isempty(index) || length(index) < 2
        index = 1:size(signals_psd,1);
    end

    freq_sel = freq(index);
    df = freq_sel(2) - freq_sel(1);

    aux = signals_psd(index,j);
    
    if strcmpi(method, 'rms')
        signals_integrated_rms(:,j) = sqrt(cumtrapz(aux)*df);
    elseif strcmpi(method, 'rms_reversed')
        aux = aux(end:-1:1);
        signals_integrated_rms(:,j) = sqrt(cumtrapz(aux)*df);
        signals_integrated_rms(:,j) = signals_integrated_rms(end:-1:1,j);
    elseif strcmpi(method, 'power')
        signals_integrated_rms(:,j) = cumtrapz(aux)*df;
    elseif strcmpi(method, 'psd')
        signals_integrated_rms(:,j) = aux;
    elseif strcmpi(method, 'sqrtpsd')
        signals_integrated_rms(:,j) = sqrt(aux);
    end
end