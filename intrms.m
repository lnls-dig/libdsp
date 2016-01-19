function [Xrms, f_sel] = intrms(x, Fs, flow, fhigh, window, noverlap)
%INTRMS   Integrated RMS-value over a frequency band.
%
%   [Xrms, f] = intrms(x, Fs, flow, fhigh, window, noverlap)
%
%   Inputs:
%       x:          time series (raw data values)
%       Fs:         sampling time (Hz) (default = 1)
%       flow:       low frequency limit for integration (default = 0)
%       fhigh:      high frequency limit for integration (default = Inf)
%       window:     windowing function for internal DFT calculation (default = rectangular window with 10% of number of points of x)
%       noverlap:   number of samples which overlap between two consecutive windows on PSD calculation
%
%   Outputs:
%       Xrms:       integrated RMS-value in frequency domain
%       f:          frequency vector

%   Copyright (C) 2014 CNPEM
%   Licensed under GNU Lesser General Public License v3.0 (LGPL)
%
%   Revisions:
%       2014-04-04 Daniel Tavares (LNLS/DIG) - Initial realease (as psdrms)
%       2016-01-19 Daniel Tavares (LNLS/DIG) - Renamed intrms and simplified the interface

if nargin < 2 || isempty(Fs)
    Fs = 1;
end
if nargin < 4 || isempty(flow)  || isempty(fhigh)
    flow = 0;
    fhigh = Inf;
end
if nargin < 5 || isempty(window)
    window = rectwin(floor(size(x,1)/10));
end
if nargin < 6 || isempty(noverlap)
    noverlap = floor(0.5*length(window));
end    

% Remove DC component
npts = size(x,1);
x = x - repmat(mean(x), npts, 1);

Pxx = zeros(ceil((size(window,1)+1)/2), size(x,2));
for j=1:size(x, 2)
    % Compute PSD (Power Spectrum Density in units of Power/Hz)
    [Pxx(:,j), f] = pwelch(x(:,j), window, noverlap, length(window), Fs, 'onesided');
    
    % Select only frequencies and PSD values at selected interval (flow <= f <= fhigh)
    index = find((f >= flow-eps) & (f <= fhigh+eps));
    if isempty(index) || length(index) < 2
        index = 1:size(Pxx,1);
    end
    f_sel = f(index);
    Pxx_sel = Pxx(index,j);
    
    % Frequency resolution
    df = f_sel(2) - f_sel(1);    
    
    % Integrate PSD and square root the result to get the integrated RMS
    Xrms(:,j) = sqrt(cumtrapz(Pxx_sel)*df);
end