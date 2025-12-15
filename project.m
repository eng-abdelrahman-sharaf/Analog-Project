
clc; clear; close all;

%% ================= AUDIO LOADING =================
try
    [m, Fs] = audioread('eric');
catch
    [m, Fs] = audioread('eric.wav');
end
if size(m,2)>1, m=m(:,1); end
N = length(m);
t = (0:N-1)/Fs;

%% ================= ORIGINAL SPECTRUM =================
M_f = fftshift(fft(m));
f = (-N/2:N/2-1)*(Fs/N);
figure; plot(f,abs(M_f));
title('Original Spectrum'); xlabel('Hz'); ylabel('|M(f)|');

%% ================= IDEAL LPF 4 kHz =================
BW = 4000;
H = (abs(f)<=BW).';
M_filt = M_f .* H;

figure; plot(f,abs(M_filt));
title('Filtered Spectrum (BW = 4 kHz)'); xlabel('Hz');

m_filt = real(ifft(ifftshift(M_filt)));
figure; plot(t,m_filt);
title('Filtered Signal (Time Domain)');
sound(m_filt,Fs); pause(3);

%% ================= EXPERIMENT 1: DSB =================
Fc = 100e3;
Fs_new = 5*Fc;

m_up = resample(m_filt,Fs_new,Fs);
t_up = (0:length(m_up)-1).'/Fs_new;
A = 2*max(abs(m_up));   % DC bias for modulation index = 0.5

% DSB-SC & DSB-TC
dsb_sc = m_up .* cos(2*pi*Fc*t_up);
dsb_tc = (A + m_up) .* cos(2*pi*Fc*t_up);

N2 = length(dsb_sc);
f2 = (-N2/2:N2/2-1)*(Fs_new/N2);

figure; plot(f2,abs(fftshift(fft(dsb_sc))));
title('DSB-SC Spectrum'); xlabel('Hz');

figure; plot(f2,abs(fftshift(fft(dsb_tc))));
title('DSB-TC Spectrum'); xlabel('Hz');

%% ================= Envelope Detection =================
env_tc = abs(hilbert(dsb_tc));
env_sc = abs(hilbert(dsb_sc));

figure; plot(env_tc(1:2000));
title('Envelope Detector Output (DSB-TC)');

figure; plot(env_sc(1:2000));
title('Envelope Detector Output (DSB-SC)');

env_tc_ds = resample(env_tc,Fs,Fs_new);
env_sc_ds = resample(env_sc,Fs,Fs_new);

sound(env_tc_ds,Fs); pause(3);
sound(env_sc_ds,Fs); pause(3);

%% ================= Coherent Detection with Noise =================
SNRs = [0 10 30];

for snr = SNRs
    noisy = add_awgn(dsb_sc,snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod,Fs,Fs_new);

    % Time-domain plot
    figure; plot(demod_ds(1:2000));
    title(['DSB-SC Coherent Detection (Time), SNR = ',num2str(snr),' dB']);

    % Frequency-domain plot (REQUIRED)
    DEM_F = fftshift(fft(demod_ds));
    f_dem = (-length(demod_ds)/2:length(demod_ds)/2-1)*(Fs/length(demod_ds));
    figure; plot(f_dem,abs(DEM_F));
    title(['DSB-SC Coherent Detection (Freq), SNR = ',num2str(snr),' dB']);
    xlabel('Hz');

    sound(demod_ds,Fs); pause(3);
end

%% ================= Frequency & Phase Error =================
freq_err = dsb_sc .* cos(2*pi*100100*t_up);  % Carrier frequency offset (Beat frequency)
phase_err = dsb_sc .* cos(2*pi*Fc*t_up + deg2rad(20));

sound(resample(freq_err,Fs,Fs_new),Fs); pause(3);
sound(resample(phase_err,Fs,Fs_new),Fs); pause(3);

%% ================= EXPERIMENT 2: SSB =================
DSB_F = fftshift(fft(dsb_sc));

% Ideal LSB filter
H_lsb = ((f2>=Fc-BW & f2<=Fc) | (f2<=-Fc+BW & f2>=-Fc)).';
SSB_F = DSB_F .* H_lsb;

figure; plot(f2,abs(SSB_F));
title('Ideal SSB-LSB Spectrum'); xlabel('Hz');

ssb = real(ifft(ifftshift(SSB_F)));

%% ================= SSB Coherent Detection (Ideal) =================
demod_ssb = ssb .* cos(2*pi*Fc*t_up);
demod_ssb_ds = resample(demod_ssb,Fs,Fs_new);

figure; plot(demod_ssb_ds(1:2000));
title('SSB Coherent Detection (Time)');

SSB_DEM_F = fftshift(fft(demod_ssb_ds));
f_ssb_dem = (-length(demod_ssb_ds)/2:length(demod_ssb_ds)/2-1)*(Fs/length(demod_ssb_ds));
figure; plot(f_ssb_dem,abs(SSB_DEM_F));
title('SSB Coherent Detection (Frequency)');

sound(demod_ssb_ds,Fs); pause(3);

%% ================= Practical Butterworth Filter =================
[b,a] = butter(4, BW/(Fs_new/2));
ssb_practical = filter(b,a,dsb_sc);

demod_p = ssb_practical .* cos(2*pi*Fc*t_up);
demod_p_ds = resample(demod_p,Fs,Fs_new);
sound(demod_p_ds,Fs); pause(3);

%% ================= SSB with Noise =================
for snr = SNRs
    noisy = add_awgn(ssb,snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod,Fs,Fs_new);

    figure; plot(demod_ds(1:2000));
    title(['SSB Coherent Detection, SNR = ',num2str(snr),' dB']);

    DEM_F = fftshift(fft(demod_ds));
    f_dem = (-length(demod_ds)/2:length(demod_ds)/2-1)*(Fs/length(demod_ds));
    figure; plot(f_dem,abs(DEM_F));
    title(['SSB Spectrum After Demodulation, SNR = ',num2str(snr),' dB']);

    sound(demod_ds,Fs); pause(3);
end

%% ================= SSB-TC + Envelope Detection =================
ssb_tc = ssb + A*cos(2*pi*Fc*t_up);
env_ssb = abs(hilbert(ssb_tc));
env_ssb_ds = resample(env_ssb,Fs,Fs_new);

figure; plot(env_ssb_ds(1:2000));
title('SSB-TC Envelope Detector Output');

sound(env_ssb_ds,Fs); pause(3);

%% ================= EXPERIMENT 3: NBFM =================
kp = 0.05;  
% NBFM condition: beta = kp * max(|m(t)|) << 1

int_m = cumsum(m_up)/Fs_new;
fm = cos(2*pi*Fc*t_up + kp*int_m);

figure; plot(f2,abs(fftshift(fft(fm))));
title('NBFM Spectrum'); xlabel('Hz');

%% ================= FM Demodulation =================
fm_diff = diff(fm);                 % Differentiator
fm_env = abs(hilbert(fm_diff));     % Envelope Detector
fm_ds = resample(fm_env,Fs,Fs_new);

figure; plot(fm_ds(1:2000));
title('NBFM Demodulated Signal');

sound(fm_ds,Fs);

%% ================= LOCAL FUNCTION =================
function y = add_awgn(x,SNRdB)
    P = mean(x.^2);
    N = P / (10^(SNRdB/10));
    y = x + sqrt(N)*randn(size(x));
end
