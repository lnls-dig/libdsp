function res = ddc(signals, nco, lpf)
    i = signals.*real(nco);
    q = signals.*imag(nco);

    ifilt = filter(lpf, i);
    qfilt = filter(lpf, q);

    res = ifilt + 1j*qfilt;
end
