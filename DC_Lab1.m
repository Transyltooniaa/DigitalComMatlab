% Parameters
fs = 1000;           % Sampling Frequency (Hz)
T = 2;               % Duration of the input signal (seconds)
t = 0:1/fs:T-1/fs;   % Time vector

% Generate Gaussian noise
noi=se_level = 0.2;   % Adjust the noise level (Standard Deviation)
noise = noise_level * randn(size(t));

% 1. Generate a sinusoid input signal with frequency 50Hz and add noise to the signal
input_signal = sin(2 * pi * 50 * t) + noise;

% Shorter duration for the reference signal
T_ref = 0.5;         % Duration of the reference signal (seconds)
t_ref = 0:1/fs:T_ref-1/fs; % Time vector for the reference signal

% 2. Generate reference signal same as input but of a shorter duration
reference_signal = sin(2 * pi * 50 * t_ref);

% 3. Flip the reference signal for the matched filter using MATLAB built-in function
matched_filter = fliplr(reference_signal);

% 4. Perform convolution operation with a flipped version of reference signal and input signal
filtered_signal = conv(input_signal, matched_filter, 'same');

% 5. Perform Cross-Correlation using MATLAB built-in function for input signal and reference signal
[correlation, lags] = xcorr(input_signal, reference_signal);

% Normalize the cross-correlation output manually
norm_factor = sqrt(sum(input_signal.^2)) * sqrt(sum(reference_signal.^2));
normalized_correlation = correlation / norm_factor;

% Align Cross-Correlation with the matched filter output signal
lag_offset = (length(reference_signal)-1); % Center of Cross-Correlation
aligned_lags = lags - lag_offset; % Shift lags for centering
aligned_correlation = normalized_correlation(lags >= 1 & lags <= length(input_signal));

% Create Time Vector for Cross-Correlation Plot
t_corr = linspace(t(1), t(end), length(aligned_correlation));

% Normalize the matched filter output and cross-correlation result to their respective max
filtered_signalNormalized = filtered_signal / max(abs(filtered_signal));
correlationNormalized = aligned_correlation / max(abs(aligned_correlation));

% Plot all the graphs in subplots
figure;

% Plot the input signal
subplot(4, 1, 1);
plot(t, input_signal);
title('Input Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the reference signal
subplot(4, 1, 2);
plot(t_ref, reference_signal);
title('Reference Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the filtered signal (Matched filter output)
subplot(4, 1, 3);
plot(t, filtered_signalNormalized);
title('Matched Filter Output (Normalized)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the cross-correlation result
subplot(4, 1, 4);
plot(t_corr, correlationNormalized);
title('Cross-Correlation (Normalized)');
xlabel('Time (s)');
ylabel('Correlation');


figure;
plot(t,filtered_signalNormalized,'b','DisplayName','Matched Filter Output');
hold on;
plot(t_corr,correlationNormalized,'r','DisplayName','Cross-Correlation Ouptut');
legend;
title('Comparision of Normalized Outputs');
xlabel('Time(s)');
ylabel('Normalized Amplitude');

