function Yq = frinterp(f,Y,fq)

lf = log(f);
lYmag = log(abs(Y));
Yph = unwrap(angle(Y));

lfq = log(fq);

lYmagi = interp1(lf,lYmag,lfq);
Yphi = interp1(lf,Yph,lfq);

Yq = exp(lYmagi).*exp(1j*Yphi);