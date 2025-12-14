clc; clear; close all;

%% =========================================================
% ANALOG COMMUNICATION PROJECT
% Experiments 1 & 2
% Audio file name: eric   (extension-less WAV)
%% =========================================================

%% =========================================================
% SAFE AUDIO LOADING (extension-less support)
%% =========================================================
try
    [m, Fs] = audioread('eric');
catch
    [m, Fs] = audioread('eric.wav');
end

if size(m,2) > 1
    m = m(:,1);   % Convert to mono if stereo
end

t = (0:length(m)-1)/Fs;

%% =========================================================
% COMMON PART – Spectrum & 4 kHz Band-Limiting
%% =========================================================

N = length(m);
M_f = fftshift(fft(m));
f = (-N/2:N/2-1)*(Fs/N);

figure;
plot(f, abs(M_f));
title('Original Message Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

BW = 4000;
H = (abs(f) <= BW).';
M_filt = M_f .* H;

figure;
plot(f, abs(M_filt));
title('Filtered Spectrum (BW = 4 kHz)');
xlabel('Frequency (Hz)');

m_filt = real(ifft(ifftshift(M_filt)));
sound(m_filt, Fs);
pause(3);

%% =========================================================
% EXPERIMENT 1 – DOUBLE SIDEBAND MODULATION
%% =========================================================

Fc = 100e3;
Fs_new = 5 * Fc;

m_up = resample(m_filt, Fs_new, Fs);
t_up = (0:length(m_up)-1).' / Fs_new;

A = 2 * max(abs(m_up));           % DC bias

% DSB-SC
dsb_sc = m_up .* cos(2*pi*Fc*t_up);

% DSB-TC
dsb_tc = (A + m_up) .* cos(2*pi*Fc*t_up);

% Frequency domain
N2 = length(dsb_sc);
f2 = (-N2/2:N2/2-1)*(Fs_new/N2);

DSB_SC_F = fftshift(fft(dsb_sc));
DSB_TC_F = fftshift(fft(dsb_tc));

figure;
plot(f2, abs(DSB_SC_F));
title('DSB-SC Spectrum');
xlabel('Frequency (Hz)');

figure;
plot(f2, abs(DSB_TC_F));
title('DSB-TC Spectrum');
xlabel('Frequency (Hz)');

%% Envelope Detection
env_sc = abs(hilbert(dsb_sc));
env_tc = abs(hilbert(dsb_tc));

env_sc_ds = resample(env_sc, Fs, Fs_new);
env_tc_ds = resample(env_tc, Fs, Fs_new);

sound(env_tc_ds, Fs);
pause(3);
sound(env_sc_ds, Fs);
pause(3);

%% Coherent Detection (DSB-SC) with Noise
SNRs = [0 10 30];

for snr = SNRs
    noisy = add_awgn(dsb_sc, snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod, Fs, Fs_new);
    sound(demod_ds, Fs);
    pause(3);
end

%% Frequency Error
carrier_err = cos(2*pi*100100*t_up);
freq_err = dsb_sc .* carrier_err;
freq_err_ds = resample(freq_err, Fs, Fs_new);
sound(freq_err_ds, Fs);
pause(3);

%% Phase Error
carrier_phase = cos(2*pi*Fc*t_up + deg2rad(20));
phase_err = dsb_sc .* carrier_phase;
phase_err_ds = resample(phase_err, Fs, Fs_new);
sound(phase_err_ds, Fs);
pause(3);

%% =========================================================
% EXPERIMENT 2 – SINGLE SIDEBAND MODULATION
%% =========================================================

dsb_sc = m_up .* cos(2*pi*Fc*t_up);
DSB_F = fftshift(fft(dsb_sc));

figure;
plot(f2, abs(DSB_F));
title('DSB-SC Spectrum (Experiment 2)');
xlabel('Frequency (Hz)');

%% Ideal SSB (LSB only)
H_ssb = (f2 >= Fc-BW & f2 <= Fc).';
SSB_F = DSB_F .* H_ssb;

figure;
plot(f2, abs(SSB_F));
title('SSB-LSB Spectrum (Ideal Filter)');
xlabel('Frequency (Hz)');

ssb = real(ifft(ifftshift(SSB_F)));

%% Coherent Detection
demod_ssb = ssb .* cos(2*pi*Fc*t_up);
demod_ssb_ds = resample(demod_ssb, Fs, Fs_new);
sound(demod_ssb_ds, Fs);
pause(3);

%% Practical Butterworth Filter
[b,a] = butter(4, BW/(Fs_new/2));
ssb_practical = filter(b, a, dsb_sc);

demod_practical = ssb_practical .* cos(2*pi*Fc*t_up);
demod_practical_ds = resample(demod_practical, Fs, Fs_new);
sound(demod_practical_ds, Fs);
pause(3);

%% SSB with Noise
for snr = SNRs
    noisy = add_awgn(ssb, snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod, Fs, Fs_new);
    sound(demod_ds, Fs);
    pause(3);
end

%% SSB-TC with Envelope Detector
ssb_tc = (A + m_up) .* cos(2*pi*Fc*t_up);
env_ssb = abs(hilbert(ssb_tc));
env_ssb_ds = resample(env_ssb, Fs, Fs_new);
sound(env_ssb_ds, Fs);

%% =========================================================
% END OF FILE
%% =========================================================

%% =========================================================
% LOCAL FUNCTION: AWGN (NO TOOLBOX REQUIRED)
%% =========================================================
function y = add_awgn(x, SNR_dB)
    signal_power = mean(x.^2);
    SNR_linear = 10^(SNR_dB/10);
    noise_power = signal_power / SNR_linear;
    noise = sqrt(noise_power) * randn(size(x));
    y = x + noise;
end

