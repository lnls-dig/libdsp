function sm = mps(s)
% sm = mps(s)
% Create minimum-phase spectrum sm from complex spectrum s
%
% Original code written by J.O. Smith (CCRMA-Stanford)
sm = exp(fft(fold(ifft(log(clipdb(s,-100))))));

function [clipped] = clipdb(s,cutoff)
% [clipped] = clipdb(s,cutoff)
% Clip magnitude of s at its maximum + cutoff in dB.
% Example: clip(s,-100) makes sure the minimum magnitude
% of s is not more than 100 dB below its maximum magnitude.
% If s is zero, nothing is done.

clipped = s;
as = abs(s);
mas = max(as(:));
if mas==0
    return
end
if cutoff >= 0
    return
end
thresh = mas*10^(cutoff/20); % dB to linear
clipped = s;
clipped(as < thresh) = thresh;

function rw = fold(r)
% Fold left wing of vector in "FFT buffer format" onto right wing

[m,n] = size(r);
if m*n ~= m+n-1
    error('fold.m: input must be a vector');
end
flipped = 0;
if (m > n)
    n = m;
    r = r.';
    flipped = 1;
end
if n < 3
    rw = r;
    return;
elseif mod(n,2)==1
    nt = (n+1)/2;
    rw = [r(1), r(2:nt) + conj(r(n:-1:nt+1)), zeros(1,n-nt)];
else
    nt = n/2;
    rf = [r(2:nt),0];
    rf = rf + conj(r(n:-1:nt+1));
    rw = [r(1), rf, zeros(1,n-nt-1)];
end

if flipped
    rw = rw.';
end