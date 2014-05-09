function [spectrum, freq] = psdrms(signals, Fs, freq_low, freq_high, window, noverlap, nfft, method)

if nargin < 8
    method = 'rms';
end

% Remove DC component
npts = size(signals,1);
signals = signals - repmat(mean(signals), npts, 1);

% Compute PSD and RMS integrated
for j=1:size(signals, 2)
    [signals_psd(:,j), freq_aux] = pwelch(signals(:,j), window, noverlap, nfft, Fs, 'onesided');
    
    index = find((freq_aux >= freq_low-eps) & (freq_aux <= freq_high+eps));
    
    if isempty(index) || length(index) < 2
        index = 1:size(signals_psd,1);
    end

    freq = freq_aux(index);
    df = freq(2) - freq(1);

    aux = signals_psd(index,j);
    
    if strcmpi(method, 'rms')
        spectrum(:,j) = sqrt(cumtrapz(aux)*df);
    elseif strcmpi(method, 'rms_reversed')
        aux = aux(end:-1:1);
        spectrum(:,j) = sqrt(cumtrapz(aux)*df);
        spectrum(:,j) = spectrum(end:-1:1,j);
    elseif strcmpi(method, 'power')
        spectrum(:,j) = cumtrapz(aux)*df;
    elseif strcmpi(method, 'psd')
        spectrum(:,j) = aux;
    elseif strcmpi(method, 'sqrtpsd')
        spectrum(:,j) = sqrt(aux);
    end
end