%% STEP 1: IIR BPF Filter Bank
clc; clear;

% Parameters
fs = 22050;          % Sampling frequency
N = 8;               % Number of channels
f_low = 200;         % Lowest center frequency
f_high = 5500;       % Highest center frequency

% Calculate center frequencies (log-spaced)
center_freqs = logspace(log10(f_low), log10(f_high), N);

% Bandwidth factor
bw_factor = 0.15;

% Storage for filter coefficients
b_all = cell(N,1);
a_all = cell(N,1);

% Create invisible figure for frequency response
fig = figure('Visible', 'off'); hold on;

% Design IIR Bandpass Filters
for i = 1:N
    fc = center_freqs(i);
    bw = fc * bw_factor;
    f1 = max(fc - bw/2, 10);
    f2 = min(fc + bw/2, fs/2 - 100);
    Wn = [f1 f2] / (fs/2);
    [b, a] = butter(2, Wn, 'bandpass');
    b_all{i} = b;
    a_all{i} = a;

    [H, f] = freqz(b, a, 1024, fs);
    plot(f, 20*log10(abs(H)));
end

title('Frequency Response of All Bandpass Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
legend(arrayfun(@(x) ['BPF ', num2str(x)], 1:N, 'UniformOutput', false));
saveas(fig, 'filter_bank_response.png');
close(fig);

%% STEP 2: Envelope Extraction
[x, fs] = audioread('greasyWW22k.wav');  % Input signal
x = x(:,1);  % Mono

envelopes = cell(N, 1);
filtered_signals = cell(N, 1);

for i = 1:N
    x_bpf = filter(b_all{i}, a_all{i}, x);
    x_mag = abs(x_bpf);
    [b_lpf, a_lpf] = butter(3, 400 / (fs/2));
    env = filter(b_lpf, a_lpf, x_mag);
    envelopes{i} = env;
    filtered_signals{i} = x_bpf;

    % Save envelope plot (invisible)
    t = (0:length(env)-1)/fs;
    fig = figure('Visible', 'off');
    plot(t, env);
    title(['Envelope for Channel ', num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 0.1]);
    saveas(fig, ['envelope_channel_' num2str(i) '.png']);
    close(fig);
end

%% STEP 3: DC Notch Filtering
dc_removed = cell(N, 1);
a = 0.98;  % Notch sharpness

for i = 1:N
    b_notch = [1, -(1 + a), 1];
    a_notch = [1, -a, 0];
    env_dc_removed = filter(b_notch, a_notch, envelopes{i});
    dc_removed{i} = env_dc_removed;

    % Save comparison plot (invisible)
    t = (0:length(env_dc_removed)-1)/fs;
    fig = figure('Visible', 'off');
    plot(t, envelopes{i}, 'b--'); hold on;
    plot(t, env_dc_removed, 'r');
    title(['Channel ', num2str(i), ' - Envelope (w/ and w/o DC)']);
    legend('Original Envelope', 'DC Removed');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 0.1]);
    saveas(fig, ['envelope_dc_removed_' num2str(i) '.png']);
    close(fig);
end

%% STEP 4: Envelope Modulation and Summation
modulated_channels = cell(N, 1);
y_sum = zeros(size(dc_removed{1}));
t = (0:length(y_sum)-1)/fs;

for i = 1:N
    fc = center_freqs(i);
    carrier = cos(2 * pi * fc * t)';
    modulated = dc_removed{i} .* carrier;
    modulated_channels{i} = modulated;
    y_sum = y_sum + modulated;

    % Save modulated plot (invisible)
    fig = figure('Visible', 'off');
    plot(t, modulated);
    title(['Modulated Signal - Channel ', num2str(i)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 0.1]);
    saveas(fig, ['modulated_channel_' num2str(i) '.png']);
    close(fig);
end

% Normalize output
y_sum = y_sum / max(abs(y_sum));


%% STEP 5: Playback, Save, and Spectrograms

%% PLAYBACK: Compare Original and Reconstructed Signals

disp('Playing original input signal...');
soundsc(x, fs);        % Play original input
pause(length(x)/fs + 1);  % Wait until done + 1s gap

disp('Playing reconstructed output signal...');
soundsc(y_sum, fs);    % Play reconstructed signal
pause(length(y_sum)/fs + 1);  % Wait until done


% ---- Save as .wav file ----
audiowrite('output_cochlear.wav', y_sum, fs);

% ---- Plot and Save Spectrogram of Original ----
figure('Visible', 'off');
spectrogram(x, 512, 256, 1024, fs, 'yaxis');
title('Original Input Signal - Spectrogram');
saveas(gcf, 'spectrogram_original.png');
close(gcf);

% ---- Plot and Save Spectrogram of Reconstructed ----
figure('Visible', 'off');
spectrogram(y_sum, 512, 256, 1024, fs, 'yaxis');
title('Synthesized Output Signal - Spectrogram');
saveas(gcf, 'spectrogram_output.png');
close(gcf);
