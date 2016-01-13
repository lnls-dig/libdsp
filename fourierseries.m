function [amp, f, ph] = fourierseries(data, Fs, window)
%FOURIERSERIES   Fourier series of a time series.
%
%   [fft_amplitude, f, fft_phase] = FOURIERSERIES(data, Fs, window)

%   Copyright (C) 2014 CNPEM
%   Licensed under GNU Lesser General Public License v3.0 (LGPL)

if any(size(data) == 1)
    data = data(:);
end

npts = size(data,1);

if nargin < 2 || isempty(Fs)
    Fs = 1;
end
if nargin < 3 || isempty(window)
    window = rectwin(npts);
end
window = window(:);

data = data.*repmat(window, 1, size(data,2));
amp = abs(fft(data))/npts;
ph = angle(fft(data));

half_npts = ceil((npts+1)/2);
amp = amp(1:half_npts, :);
ph = ph(1:half_npts, :);

if rem(npts,2) > 0
    amp = [amp(1,:); 2*amp(2:end,:)];
else
    amp = [amp(1,:); 2*amp(2:end-1,:); amp(end,:)];
end

df = Fs/npts;

f = 0:df:(half_npts-1)*df;
