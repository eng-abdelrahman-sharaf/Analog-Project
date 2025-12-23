clc; clear; close all;

%% audio loading
try
    [m, Fs] = audioread('eric');
catch
    [m, Fs] = audioread('eric.wav');
end
if size(m,2)>1, m=m(:,1); end
N = length(m);
t = (0:N-1)/Fs;

%% original signal (time domain)
figure; plot(t,m);
title('Original Audio Signal (Time Domain)'); 
xlabel('Time (s)'); ylabel('Amplitude');
ylim([-.25 .25])

%% original spectrum
M_f = fftshift(fft(m));
f = (-N/2:N/2-1)*(Fs/N);
figure; plot(f/1e3,abs(M_f));
title('Original Spectrum'); xlabel('KHz');
xticks(-16:4:16)

%% ideal lpf 4 khz
BW = 4000;
H = (abs(f)<=BW).';
M_filt = M_f .* H;

figure; plot(f/1e3,abs(M_filt));
title('Filtered Spectrum (BW = 4 kHz)'); xlabel('KHz');
xticks(-16:4:16)

m_filt = real(ifft(ifftshift(M_filt)));
figure; plot(t,m_filt);
title('Filtered Signal (Time Domain)');
sound(m_filt,Fs); pause(3);
ylim([-.25 .25])
%% experiment 1: dsb
Fc = 100e3;
Fs_new = 5*Fc;

m_up = resample(m_filt,Fs_new,Fs);
t_up = (0:length(m_up)-1).'/Fs_new;
A = 2*max(abs(m_up));   % dc bias for modulation index = 0.5

% dsb-sc & dsb-tc
dsb_sc = m_up .* cos(2*pi*Fc*t_up);
dsb_tc = (A + m_up) .* cos(2*pi*Fc*t_up);

N2 = length(dsb_sc);
f2 = (-N2/2:N2/2-1)*(Fs_new/N2);

figure; plot(f2/1e3,abs(fftshift(fft(dsb_sc))));
title('DSB-SC Spectrum'); xlabel('KHz');

figure; plot(f2/1e3,abs(fftshift(fft(dsb_tc))));
title('DSB-TC Spectrum'); xlabel('KHz');

%% envelope detection
env_tc = abs(hilbert(dsb_tc));
env_sc = abs(hilbert(dsb_sc));

% remove dc gain
env_tc = env_tc - mean(env_tc);
env_sc = env_sc - mean(env_sc);

env_tc_ds = resample(env_tc,Fs,Fs_new);
env_sc_ds = resample(env_sc,Fs,Fs_new);

% plot envelope waveforms
t_env = (0:length(env_tc_ds)-1)/Fs;

figure; 
plot(t_env, env_tc_ds);
title('Envelope Detector Output (DSB-TC)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;
ylim([-.25 .25])

figure; 
plot(t_env, env_sc_ds);
title('Envelope Detector Output (DSB-SC)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

sound(env_tc_ds,Fs); pause(3);
sound(env_sc_ds,Fs); pause(3);

%% coherent detection with noise
SNRs = [0 10 30]; % Signal to Noise Ratio

for snr = SNRs
    noisy = add_awgn(dsb_sc,snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod,Fs,Fs_new);

    % time-domain plot
    t_demod = (0:length(demod_ds)-1)/Fs;
    figure; 
    plot(t_demod, demod_ds);
    title(['DSB-SC Coherent Detection (Time), SNR = ',num2str(snr),' dB']);
    ylim([-.25 .25])

    % frequency-domain plot (required)
    DEM_F = fftshift(fft(demod_ds));
    f_dem = (-length(demod_ds)/2:length(demod_ds)/2-1)*(Fs/length(demod_ds));
    figure; plot(f_dem/1e3,abs(DEM_F));
    title(['DSB-SC Coherent Detection (Freq), SNR = ',num2str(snr),' dB']);
    xlabel('KHz');
    xticks(-16:4:16)

    sound(demod_ds,Fs); pause(3);
end


%% coherent detection with noise
SNRs = [0 10 30]; % Signal to Noise Ratio

for snr = SNRs
    noisy = add_awgn(dsb_sc,snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    

    % time-domain plot
    t_demod = (0:length(demod)-1)/Fs;
    figure; 
    plot(t_demod, demod);
    title(['DSB-SC Coherent Detection (Time), SNR = ',num2str(snr),' dB']);
    ylim([-.25 .25])

    % frequency-domain plot (required)
    DEM_F = fftshift(fft(demod));
    f_dem = (-length(demod)/2:length(demod)/2-1)*(Fs/length(demod));
    figure; plot(f_dem/1e3,abs(DEM_F));
    title(['DSB-SC Coherent Detection (Freq), SNR = ',num2str(snr),' dB']);
    xlabel('KHz');
    xticks(-16:4:16)

    demod_ds = resample(demod,Fs,Fs_new);
    sound(demod_ds,Fs); pause(3);
end
%% frequency error
freq_err = dsb_sc .* cos(2*pi*100100*t_up);  % Carrier frequency offset (Beat frequency)

% time-domain plot
t_demod = (0:length(demod_ds)-1)/Fs;
figure; 
plot(t_demod, demod_ds);
title(['DSB-SC Coherent Detection (Time), freq error = +0.1 KHz']);
ylim([-.25 .25])

% frequency-domain plot (required)
DEM_F = fftshift(fft(demod_ds));
f_dem = (-length(demod_ds)/2:length(demod_ds)/2-1)*(Fs/length(demod_ds));
figure; plot(f_dem/1e3,abs(DEM_F));
title(['DSB-SC Coherent Detection (Freq), freq error = +0.1 KHz']);
xlabel('KHz');
xticks(-16:4:16)

sound(resample(freq_err,Fs,Fs_new),Fs); pause(3);
%% phase error
phase_err = dsb_sc .* cos(2*pi*Fc*t_up + deg2rad(20));

% time-domain plot
t_demod = (0:length(demod_ds)-1)/Fs;
figure; 
plot(t_demod, demod_ds);
title(['DSB-SC Coherent Detection (Time), phase error = 20 deg']);
ylim([-.25 .25])

% frequency-domain plot (required)
DEM_F = fftshift(fft(demod_ds));
f_dem = (-length(demod_ds)/2:length(demod_ds)/2-1)*(Fs/length(demod_ds));
figure; plot(f_dem/1e3,abs(DEM_F));
title(['DSB-SC Coherent Detection (Freq), phase error = 20 deg']);
xlabel('KHz');
xticks(-16:4:16)

sound(resample(phase_err,Fs,Fs_new),Fs); pause(3);

%% experiment 2: ssb
DSB_F = fftshift(fft(dsb_sc));

H_lsb = ((f2 >= Fc-BW & f2 <= Fc) | (f2 <= -Fc+BW & f2 >= -Fc)).';
SSB_F = DSB_F .* H_lsb;

figure;
plot(f2/1e3, abs(SSB_F));
title('Ideal SSB-LSB Spectrum');
xlabel('Frequency (KHz)');

ssb = real(ifft(ifftshift(SSB_F)));

%% ssb coherent detection
demod_ssb = ssb .* cos(2*pi*Fc*t_up);
demod_ssb_ds = resample(demod_ssb, Fs, Fs_new);

t_ssb = (0:length(demod_ssb_ds)-1)/Fs;

figure;
plot(t_ssb, demod_ssb_ds);
title('SSB Coherent Detection (Time Domain)');
xlabel('Time (s)');

SSB_DEM_F = fftshift(fft(demod_ssb_ds));
f_ssb_dem = (-length(demod_ssb_ds)/2:length(demod_ssb_ds)/2-1)*(Fs/length(demod_ssb_ds));

figure;
plot(f_ssb_dem/1e3, abs(SSB_DEM_F));
title('SSB Coherent Detection (Frequency Domain)');
xlabel('Frequency (KHz)');

sound(demod_ssb_ds, Fs); pause(3);

%% practical butterworth ssb filter
Wp = [(Fc-BW)/(Fs_new/2) Fc/(Fs_new/2)]; % normalized Passband edge frequencies
[b,a] = butter(4, Wp, 'bandpass');
ssb_practical = filter(b, a, dsb_sc);

demod_p = ssb_practical .* cos(2*pi*Fc*t_up);
demod_p_ds = resample(demod_p, Fs, Fs_new);

t_dem_p = (0:length(demod_p_ds)-1)/Fs;
figure;
plot(t_dem_p, demod_p_ds);
title('Practical SSB Demodulated Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');

DEM_P_F = fftshift(fft(demod_p_ds));
f_dem_p = (-length(demod_p_ds)/2:length(demod_p_ds)/2-1)*(Fs/length(demod_p_ds));
figure;
plot(f_dem_p/1e3, abs(DEM_P_F));
title('Practical SSB Demodulated Spectrum');
xlabel('Frequency (KHz)'); ylabel('|X(f)|');

sound(demod_p_ds, Fs); pause(3);

%% ssb with noise
for snr = SNRs
    noisy = add_awgn(ssb, snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod, Fs, Fs_new);

    t_dem = (0:length(demod_ds)-1)/Fs;

    figure;
    plot(t_dem, demod_ds);
    title(['SSB Coherent Detection, SNR = ', num2str(snr), ' dB']);
    xlabel('Time (s)');

    DEM_F = fftshift(fft(demod_ds));
    f_dem = (-length(demod_ds)/2:length(demod_ds)/2-1)*(Fs/length(demod_ds));

    figure;
    plot(f_dem/1e3, abs(DEM_F));
    title(['SSB Spectrum After Demodulation, SNR = ', num2str(snr), ' dB']);
    xlabel('Frequency (KHz)');

    sound(demod_ds, Fs); pause(3);
end

%% ssb-tc + envelope detection
ssb_tc = ssb + A*cos(2*pi*Fc*t_up);
env_ssb = abs(hilbert(ssb_tc));
env_ssb_ds = resample(env_ssb, Fs, Fs_new);

t_env = (0:length(env_ssb_ds)-1)/Fs;

figure;
plot(t_env, env_ssb_ds);
title('SSB-TC Envelope Detector Output');
xlabel('Time (s)');

sound(env_ssb_ds, Fs); pause(3);

%% experiment 3: nbfm
kf = 0.05;  % Narrowband condition: beta << 1
Ac = 1;

m_integration = cumsum(m_up)/Fs_new;

phi = 2*pi*kf * m_integration;

NBFM = Ac * cos(2*pi * Fc * t_up + phi);

figure; plot(f2/1e3, abs(fftshift(fft(NBFM))));
title('NBFM Spectrum');
xlabel('Frequency (KHz)');

%% correct nbfm demodulation

NBFM_diff = diff(NBFM)*Fs_new;         % Differentiator
NBFM_env = abs(hilbert(NBFM_diff));    % envelope detector

NBFM_env = NBFM_env - mean(NBFM_env);  % remove dc

NBFM_ds = resample(NBFM_env, Fs, Fs_new);
t_ds = (0:length(NBFM_ds)-1)/Fs;

figure;
plot(t_ds, NBFM_ds);
title('NBFM Demodulated Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;
ylim([-.16 .16]);

sound(NBFM_ds, Fs);

%% local function
function y = add_awgn(x, SNRdB)
    P = mean(x.^2);
    N = P / (10^(SNRdB/10));
    y = x + sqrt(N) * randn(size(x));
end
