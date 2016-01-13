function [y,t] = fourierseries2time(amp, ph, Fs, npts)

if nargin < 4
    npts = 2*length(f);
end

Y = amp.*exp(1j*ph);

Y = Y(:);

ffty = Y*npts/2;
ffty(1) = 2*ffty(1);

if npts == 2*(length(Y)-1)
    ffty(end) = 2*ffty(end);
end
    
ffty = [ffty; zeros(npts-length(Y),1)];
y = ifft(ffty, 'symmetric');

t = (0:npts-1)/Fs;