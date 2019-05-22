function [u, t, stable_interval, nperiods, f_actual] = sinsum(param)

nspecs = length(param.f_spec);

if isscalar(param.a_spec)
    param.a_spec = repmat(param.a_spec,1,nspecs);
end
if ~isfield(param, 'ntries')
    param.ntries = 10;
end

npts_wait = ceil(param.t_wait*param.Fs);
npts_stable = ceil(param.t_spec*param.Fs);

% Build transient waveforms (fadein, fadeout) using integrated Hann window w(t):
% w(t) = (t-sin(2*pi*t)/(2*pi), t = [0..1]

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
nperiods = zeros(nspecs,1);
for i=1:length(param.f_spec)
    if param.f_spec(i) > 0
        nperiods(i) = round(param.f_spec(i)/param.Fs*npts_stable);
        if nperiods(i) < 1
            nperiods(i) = 1;
        end
        f_actual(i) = nperiods(i)*param.Fs/npts_stable;
    else
        error('All excitation frequencies must be greater than 0.');
    end
end

npts_total = npts_wait + npts_stable + npts_fadein + npts_fadeout;

% Initialize frequency and amplitude arrays and auxiliary variables
a_fade = ones(npts_total, 1);

% Fadein
if isfield(param, 't_fadeout')
    a_fade(1:npts_fadein) = fadein_waveform;
end

% Fadeout
if isfield(param, 't_fadeout')
    a_fade(end-npts_fadeout+1:end) = fadeout_waveform;
end

stable_interval = npts_fadein + npts_wait + [1; npts_stable];

% Build time (t) and excitation (u) arrays
t = (0:length(a_fade)-1)'/param.Fs;
% Pre-compute ometa*t values for all sinusoidals - random phases will be
% added later
wt = t*2*pi*f_actual';
a = repmat(param.a_spec(:)', npts_total, 1);
crest_min = Inf;
for i=1:param.ntries
    ph = repmat(rand(1,length(f_actual)), npts_total, 1);    
    u_temp = sum(a.*sin(wt + ph), 2);
    
    % Search for sum of sinusoidals with minimum crest factor
    crest = peak2rms(u_temp(stable_interval(1):stable_interval(end)));
    if crest < crest_min
        u = u_temp;
        crest_min = crest;
    end
end

% Apply fadein and fadeout amplitude modulation
u = a_fade.*u;