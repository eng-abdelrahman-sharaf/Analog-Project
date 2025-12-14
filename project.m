
clear; clc; close all;

% Resolve audio file path (handle missing extension)
projectDir = fileparts(mfilename('fullpath'));
baseName = 'eric';
candidatePaths = { fullfile(projectDir, baseName), ...
				   fullfile(projectDir, [baseName '.wav']) };

audioPath = '';
for i = 1:numel(candidatePaths)
	if exist(candidatePaths{i}, 'file')
		audioPath = candidatePaths{i};
		break;
	end
end

if isempty(audioPath)
	error('Audio file "eric" not found in project directory.');
end

% Read audio
[x, fs] = audioread(audioPath);

% Ensure column vector
if size(x,2) > 1
	% Convert to mono if stereo: mean of channels
	x = mean(x, 2);
end

% Basic stats
dur = numel(x)/fs;
fprintf('Loaded %s | Fs = %d Hz | Duration = %.2f s | Channels: 1\n', audioPath, fs, dur);

% Normalize to -1..1 without clipping
mx = max(abs(x));
if mx > 0
	xNorm = x ./ mx;
else
	xNorm = x;
end

% Optional simple denoise (light spectral gating)
applyDenoise = true;
if applyDenoise
	% Parameters
	frameLen = round(0.03*fs);   % 30 ms
	hop = round(0.01*fs);        % 10 ms
	win = hamming(frameLen, 'periodic');
	nfft = 2^nextpow2(frameLen);

	% Estimate noise floor from first 0.25s
	noiseSamples = min(numel(xNorm), round(0.25*fs));
	noiseSeg = xNorm(1:noiseSamples);
	% Average magnitude spectrum over noise frames
	idx = 1;
	noiseMag = zeros(nfft/2+1,1);
	count = 0;
	while idx+frameLen-1 <= numel(noiseSeg)
		frame = noiseSeg(idx:idx+frameLen-1) .* win;
		X = fft(frame, nfft);
		mag = abs(X(1:nfft/2+1));
		noiseMag = noiseMag + mag;
		count = count + 1;
		idx = idx + hop;
	end
	if count > 0
		noiseMag = noiseMag / count;
	end

	% Process full signal with spectral gate
	y = zeros(size(xNorm));
	idx = 1;
	out = zeros(numel(xNorm)+nfft,1);
	outIdx = 1;
	while idx+frameLen-1 <= numel(xNorm)
		frame = xNorm(idx:idx+frameLen-1) .* win;
		X = fft(frame, nfft);
		mag = abs(X(1:nfft/2+1));
		phase = angle(X(1:nfft/2+1));
		% Gate: attenuate bins below threshold
		thresh = 1.5 * noiseMag;
		gain = mag ./ max(mag, thresh);
		magG = mag .* gain;
		Xh = magG .* exp(1j*phase);
		% Reconstruct full spectrum (Hermitian symmetry)
		Xfull = [Xh; conj(Xh(end-1:-1:2))];
		frameRec = real(ifft(Xfull));
		% Overlap-add
		out(outIdx:outIdx+nfft-1) = out(outIdx:outIdx+nfft-1) + frameRec;
		outIdx = outIdx + hop;
		idx = idx + hop;
	end
	y = out(1:numel(xNorm));
	% Renormalize
	my = max(abs(y));
	if my > 0
		y = y ./ my;
	end
else
	y = xNorm;
end

% Playback (comment out if not desired)
fprintf('Playing original...\n');
sound(xNorm, fs);
pause(min(dur, 3)); % preview 3 seconds
fprintf('Playing processed...\n');
sound(y, fs);

% Visualizations
t = (0:numel(xNorm)-1)/fs;

figure('Name','Eric: Waveform');
subplot(2,1,1);
plot(t, xNorm);
xlabel('Time (s)'); ylabel('Amplitude'); title('Original (normalized)'); grid on;
subplot(2,1,2);
plot(t, y);
xlabel('Time (s)'); ylabel('Amplitude'); title('Processed (denoised)'); grid on;

% Spectrum (magnitude)
N = 2^nextpow2(numel(xNorm));
Xmag = abs(fft(xNorm, N));
Ymag = abs(fft(y, N));
f = (0:N-1)*(fs/N);

figure('Name','Eric: Spectrum');
plot(f(1:floor(N/2)), 20*log10(Xmag(1:floor(N/2))+eps), 'b'); hold on;
plot(f(1:floor(N/2)), 20*log10(Ymag(1:floor(N/2))+eps), 'r');
legend('Original','Processed'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;
title('Magnitude Spectrum');
xlim([0, fs/2]);

% Save processed file
outPath = fullfile(projectDir, 'eric_processed.wav');
audiowrite(outPath, y, fs);
fprintf('Saved processed audio to %s\n', outPath);

