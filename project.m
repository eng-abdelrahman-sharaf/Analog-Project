clc; clear; close all;

%% ================= AUDIO LOADING =================
try
    [m, Fs] = audioread('eric');
catch
    [m, Fs] = audioread('eric.wav');
end

if size(m,2)>1
    m = m(:,1);
end

N = length(m);
t = (0:N-1)/Fs;

%% ================= ORIGINAL SIGNAL =================
figure; plot(t,m);
title('Original Audio Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');

M_f = fftshift(fft(m));
f = (-N/2:N/2-1)*(Fs/N);

figure; plot(f,abs(M_f));
title('Original Spectrum'); xlabel('Hz'); ylabel('|M(f)|');

%% ================= IDEAL LPF =================
BW = 4000;
H = (abs(f)<=BW).';
M_filt = M_f .* H;

figure; plot(f,abs(M_filt));
title('Filtered Spectrum (BW = 4 kHz)');

m_filt = real(ifft(ifftshift(M_filt)));
figure; plot(t,m_filt);
title('Filtered Signal (Time Domain)');

play_sound(m_filt, Fs, 'Filtered baseband signal');
pause(3);

%% ================= DSB MODULATION =================
Fc = 100e3;
Fs_new = 5*Fc;

m_up = resample(m_filt, Fs_new, Fs);
t_up = (0:length(m_up)-1).' / Fs_new;

A = 2*max(abs(m_up));   % modulation index = 0.5

dsb_sc = m_up .* cos(2*pi*Fc*t_up);
dsb_tc = (A + m_up) .* cos(2*pi*Fc*t_up);

N2 = length(dsb_sc);
f2 = (-N2/2:N2/2-1)*(Fs_new/N2);

figure; plot(f2,abs(fftshift(fft(dsb_sc))));
title('DSB-SC Spectrum');

figure; plot(f2,abs(fftshift(fft(dsb_tc))));
title('DSB-TC Spectrum');

%% ================= ENVELOPE DETECTION =================
env_tc = abs(hilbert(dsb_tc));
env_sc = abs(hilbert(dsb_sc));

env_tc_ds = resample(env_tc, Fs, Fs_new);
env_sc_ds = resample(env_sc, Fs, Fs_new);

t_env = (0:length(env_tc_ds)-1)/Fs;

figure; plot(t_env, env_tc_ds);
title('Envelope Detector Output (DSB-TC)');

figure; plot(t_env, env_sc_ds);
title('Envelope Detector Output (DSB-SC)');

play_sound(env_tc_ds, Fs, 'Envelope detected DSB-TC');
pause(3);

play_sound(env_sc_ds, Fs, 'Envelope detected DSB-SC');
pause(3);

%% ================= COHERENT DETECTION WITH NOISE =================
SNRs = [0 10 30];

for snr = SNRs
    noisy = add_awgn(dsb_sc, snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod, Fs, Fs_new);

    t_demod = (0:length(demod_ds)-1)/Fs;

    figure; plot(t_demod, demod_ds);
    title(['DSB-SC Coherent Detection (Time), SNR = ', num2str(snr), ' dB']);

    DEM_F = fftshift(fft(demod_ds));
    f_dem = (-length(demod_ds)/2:length(demod_ds)/2-1)*(Fs/length(demod_ds));

    figure; plot(f_dem, abs(DEM_F));
    title(['DSB-SC Coherent Detection (Freq), SNR = ', num2str(snr), ' dB']);

    play_sound(demod_ds, Fs, ...
        sprintf('DSB-SC coherent demodulation (SNR = %d dB)', snr));
    pause(3);
end

%% ================= FREQUENCY & PHASE ERROR =================
freq_err = dsb_sc .* cos(2*pi*100100*t_up);
phase_err = dsb_sc .* cos(2*pi*Fc*t_up + deg2rad(20));

play_sound(resample(freq_err, Fs, Fs_new), Fs, 'DSB-SC with frequency offset');
pause(3);

play_sound(resample(phase_err, Fs, Fs_new), Fs, 'DSB-SC with phase error (20 deg)');
pause(3);

%% ================= SSB GENERATION =================
DSB_F = fftshift(fft(dsb_sc));

H_lsb = ((f2 >= Fc-BW & f2 <= Fc) | ...
         (f2 <= -Fc+BW & f2 >= -Fc)).';

SSB_F = DSB_F .* H_lsb;

figure; plot(f2, abs(SSB_F));
title('Ideal SSB-LSB Spectrum');

ssb = real(ifft(ifftshift(SSB_F)));

%% ================= SSB COHERENT DEMOD =================
demod_ssb = ssb .* cos(2*pi*Fc*t_up);
demod_ssb_ds = resample(demod_ssb, Fs, Fs_new);

t_ssb = (0:length(demod_ssb_ds)-1)/Fs;

figure; plot(t_ssb, demod_ssb_ds);
title('SSB Coherent Detection (Time Domain)');

play_sound(demod_ssb_ds, Fs, 'SSB coherent demodulation');
pause(3);

%% ================= PRACTICAL SSB FILTER =================
Wp = [(Fc-BW)/(Fs_new/2) Fc/(Fs_new/2)];
[b,a] = butter(4, Wp, 'bandpass');

ssb_practical = filter(b, a, dsb_sc);
demod_p = ssb_practical .* cos(2*pi*Fc*t_up);
demod_p_ds = resample(demod_p, Fs, Fs_new);

figure; plot((0:length(demod_p_ds)-1)/Fs, demod_p_ds);
title('Practical SSB Demodulated Signal');

play_sound(demod_p_ds, Fs, 'Practical SSB demodulated signal');
pause(3);

%% ================= SSB WITH NOISE =================
for snr = SNRs
    noisy = add_awgn(ssb, snr);
    demod = noisy .* cos(2*pi*Fc*t_up);
    demod_ds = resample(demod, Fs, Fs_new);

    figure; plot((0:length(demod_ds)-1)/Fs, demod_ds);
    title(['SSB Coherent Detection, SNR = ', num2str(snr), ' dB']);

    play_sound(demod_ds, Fs, ...
        sprintf('SSB coherent demodulation (SNR = %d dB)', snr));
    pause(3);
end

%% ================= SSB-TC ENVELOPE =================
ssb_tc = ssb + A*cos(2*pi*Fc*t_up);
env_ssb = abs(hilbert(ssb_tc));
env_ssb_ds = resample(env_ssb, Fs, Fs_new);

figure; plot((0:length(env_ssb_ds)-1)/Fs, env_ssb_ds);
title('SSB-TC Envelope Detector Output');

play_sound(env_ssb_ds, Fs, 'SSB-TC envelope detector output');
pause(3);

%% ================= NBFM =================
kf = 0.05;
Ac = 1;

m_int = cumsum(m_up)/Fs_new;
phi = 2*pi*kf * m_int;

NBFM = Ac * cos(2*pi*Fc*t_up + phi);

figure; plot(f2/1e3, abs(fftshift(fft(NBFM))));
title('NBFM Spectrum');

%% ================= NBFM DEMOD =================
NBFM_diff = diff(NBFM)*Fs_new;
NBFM_env = abs(hilbert(NBFM_diff));
NBFM_env = NBFM_env - mean(NBFM_env);

NBFM_ds = resample(NBFM_env, Fs, Fs_new);

figure; plot((0:length(NBFM_ds)-1)/Fs, NBFM_ds);
title('NBFM Demodulated Signal');

play_sound(NBFM_ds, Fs, 'NBFM demodulated signal');

%% ================= LOCAL FUNCTIONS =================
function play_sound(x, Fs, label)
    fprintf('\n▶ Playing: %s\n', label);
    sound(x, Fs);
end

function y = add_awgn(x, SNRdB)
    P = mean(x.^2);
    N = P / (10^(SNRdB/10));
    y = x + sqrt(N) * randn(size(x));
end
