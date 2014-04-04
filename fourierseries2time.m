function [y,t] = fourierseries2time(f, coeff, tfinal)

coeff = coeff(:);

ffty = 2*coeff*length(coeff);
ffty(1) = coeff(1);
ffty = [ffty; zeros(length(ffty),1)];
y = ifft(ffty, 'symmetric');

npts = length(y);

Fs = 2*f(end)/npts*(npts+1);

t = (0:npts-1)/Fs;

if nargin > 2
    tindex = find(t <= tfinal);
    t = t(tindex);
    y = y(tindex);    
end
