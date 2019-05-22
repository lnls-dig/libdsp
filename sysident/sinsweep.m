function [u, t, stable_intervals, periods] = sinsweep(param)

nspecs = length(param.f_spec);

if isscalar(param.a_spec)
    param.a_spec = repmat(param.a_spec,1,nspecs);
end
if isscalar(param.t_spec)
    param.t_spec = repmat(param.t_spec,1,nspecs);
end

npts_wait = ceil(param.t_wait*param.Fs);

% Build transient waveforms (fadein, fadeout, frequency/amplitude 
% transitions) using integrated Hann window w(t): 
% w(t) = (t-sin(2*pi*t)/(2*pi), t = [0..1]

% Transition
npts_transition = ceil(param.t_transition*param.Fs);
t = linspace(0, 1, npts_transition)';
transition_waveform = t-sin(2*pi*t)/2/pi;

% Fadein
if isfield(param, 't_fadein')
    npts_fadein = ceil(param.t_fadein*param.Fs);
    t = linspace(0, 1, npts_fadein)';
    fadein_waveform = t-sin(2*pi*t)/2/pi;
else
    npts_fadein = 0;
end

% Fadeout
if isfield(param, 't_fadeout')
    npts_fadeout = ceil(param.t_fadeout*param.Fs);
    t = linspace(1, 0, npts_fadeout)';
    fadeout_waveform = t-sin(2*pi*t)/2/pi;
else
    npts_fadeout = 0;
end

% Calculate actual excitation frequencies that will result in an integer
% number of cycles for an integer number of samples
f_actual = zeros(nspecs,1);
npts_stable = zeros(nspecs,1);
periods = zeros(nspecs,1);
for i=1:length(param.f_spec)
    if param.f_spec(i) > 0
        periods(i) = ceil(param.Fs/param.f_spec(i));
        f_actual(i) = param.Fs/periods(i);
        nperiods = ceil(param.t_spec(i)*param.Fs/periods(i));
        npts_stable(i) = nperiods*periods(i);
    else
        error('All excitation frequencies must be greater than 0.');
    end
end

npts_total = nspecs*(npts_transition+npts_wait) + sum(npts_stable) + npts_fadein + npts_fadeout;

% Initialize frequency and amplitude arrays and auxiliary variables
f = zeros(npts_total, 1);
a = zeros(npts_total, 1);
if isfield(param, 't_fadein')
    f_prev = f_actual(1);
    a_prev = param.a_spec(1);
    marker = npts_fadein;
else
    f_prev = 0;
    a_prev = 0;
    marker = 0;
end

% Fadein
if isfield(param, 't_fadeout')
    f(1:npts_fadein) = fadein_waveform*f_actual(1);
    a(1:npts_fadein) = fadein_waveform*param.a_spec(1);
end

% Fadeout
if isfield(param, 't_fadeout')
    f(end-npts_fadeout+1:end) = fadeout_waveform*f_actual(end);
    a(end-npts_fadeout+1:end) = fadeout_waveform*param.a_spec(end);
end

% Frequency and amplitude transitions, wait times and stable periods
stable_intervals = zeros(2, nspecs);
for i=1:length(param.f_spec)
    f(marker + (1:npts_transition + npts_wait + npts_stable(i))) = [transition_waveform*(f_actual(i)-f_prev)+f_prev; repmat(f_actual(i), npts_stable(i) + npts_wait, 1)];
    a(marker + (1:npts_transition + npts_wait + npts_stable(i))) = [transition_waveform*(param.a_spec(i)-a_prev)+a_prev; repmat(param.a_spec(i), npts_stable(i) + npts_wait, 1)];
    stable_intervals(:,i) = marker + npts_transition  + npts_wait + [1 npts_stable(i)];
    marker = marker + npts_transition + npts_wait + npts_stable(i);
    f_prev = f_actual(i);
    a_prev = param.a_spec(i);
end

% Build time (t) and excitation (u) arrays
dt = 1/param.Fs;
t = (0:length(f)-1)'/param.Fs;
u = a.*sin(cumsum(2*pi*f)*dt);