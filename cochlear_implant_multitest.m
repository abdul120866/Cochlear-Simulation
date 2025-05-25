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

% === MASTER TEST LOOP ===
test_files = {
    'greasyWW22k.wav', 'female';
    'catsdogs22k.wav', 'male';
    'furElise22k.wav', 'music'
};

for f = 1:size(test_files, 1)
    filename = test_files{f, 1};
    label = test_files{f, 2};

    fprintf('\n===== Testing: %s (%s) =====\n', filename, label);

    %% STEP 2: Envelope Extraction
    [x, fs] = audioread(filename);
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

        % Save envelope plot
        t = (0:length(env)-1)/fs;
        fig = figure('Visible', 'off');
        plot(t, env);
        title(['Envelope - ', label, ' - Ch ', num2str(i)]);
        xlabel('Time (s)'); ylabel('Amplitude'); xlim([0 0.1]);
        saveas(fig, ['envelope_', label, '_ch' num2str(i) '.png']);
        close(fig);
    end

    %% STEP 3: DC Notch Filtering
    dc_removed = cell(N, 1);
    a = 0.98;

    for i = 1:N
        b_notch = [1, -(1 + a), 1];
        a_notch = [1, -a, 0];
        env_dc_removed = filter(b_notch, a_notch, envelopes{i});
        dc_removed{i} = env_dc_removed;

        % Save DC-removed comparison
        t = (0:length(env_dc_removed)-1)/fs;
        fig = figure('Visible', 'off');
        plot(t, envelopes{i}, 'b--'); hold on;
        plot(t, env_dc_removed, 'r');
        title(['DC Removed - ', label, ' - Ch ', num2str(i)]);
        legend('Original', 'DC Removed');
        xlabel('Time (s)'); ylabel('Amplitude'); xlim([0 0.1]);
        saveas(fig, ['dc_removed_', label, '_ch' num2str(i) '.png']);
        close(fig);
    end

    %% STEP 4: Modulation and Summation
    modulated_channels = cell(N, 1);
    y_sum = zeros(size(dc_removed{1}));
    t = (0:length(y_sum)-1)/fs;

    for i = 1:N
        fc = center_freqs(i);
        carrier = cos(2 * pi * fc * t)';
        modulated = dc_removed{i} .* carrier;
        modulated_channels{i} = modulated;
        y_sum = y_sum + modulated;

        % Save modulated signal
        fig = figure('Visible', 'off');
        plot(t, modulated);
        title(['Modulated - ', label, ' - Ch ', num2str(i)]);
        xlabel('Time (s)'); ylabel('Amplitude'); xlim([0 0.1]);
        saveas(fig, ['modulated_', label, '_ch' num2str(i) '.png']);
        close(fig);
    end

    y_sum = y_sum / max(abs(y_sum));  % Normalize

    %% STEP 5: Output and Playback
    audiowrite(['output_', label, '.wav'], y_sum, fs);

    figure('Visible', 'off');
    spectrogram(x, 512, 256, 1024, fs, 'yaxis');
    title(['Original Spectrogram - ', label]);
    saveas(gcf, ['spec_orig_', label, '.png']);
    close(gcf);

    figure('Visible', 'off');
    spectrogram(y_sum, 512, 256, 1024, fs, 'yaxis');
    title(['Output Spectrogram - ', label]);
    saveas(gcf, ['spec_output_', label, '.png']);
    close(gcf);

    fprintf('Now playing original (%s)...\n', label);
    soundsc(x, fs); pause(length(x)/fs + 1);
    fprintf('Now playing reconstructed output (%s)...\n', label);
    soundsc(y_sum, fs); pause(length(y_sum)/fs + 1);
end
