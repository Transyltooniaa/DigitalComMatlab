% Parameters
N = 10000; % Number of symbols
SNR_dB = 20; % Signal-to-Noise Ratio in dB
EbNo = 10^(SNR_dB/10); % Convert SNR from dB to linear scale

% Generate random binary data
data = randi([0 1], 1, N); 

% BPSK Modulation (0 -> -1, 1 -> 1)
modulated_signal = 2*data - 1;  % BPSK: 1 -> 1, 0 -> -1

% Generate Rayleigh magnitude and uniform phase
magnitude = sqrt(-2 * log(rand(1, N)));  % Rayleigh-distributed magnitude
phase = 2 * pi * rand(1, N);  % Uniformly distributed phase between [0, 2*pi]

% Construct the manual complex channel h = magnitude * exp(j*phase)
h = magnitude .* exp(1i * phase);  % Rayleigh magnitude and uniform phase

% Received signal: fading channel (h) applied to modulated signal
received_signal = h .* modulated_signal; 

% Additive White Gaussian Noise (AWGN)
noise = (randn(1, N) + 1i * randn(1, N)) / sqrt(2*EbNo); 
received_signal_with_noise = received_signal + noise;

% Equalization: remove channel effect by dividing by h
equalized_signal = received_signal_with_noise ./ h;

% BPSK Demodulation after equalization (real part decision)
demodulated_signal = real(equalized_signal) > 0;  % Decision based on real part

% Bit Error Rate (BER)
BER = sum(data ~= demodulated_signal) / N;  % Compare transmitted and received bits
disp(['Bit Error Rate (BER) after Equalization: ' num2str(BER)]);

% Plot channel: magnitude (Rayleigh) and phase (uniform)
figure;
subplot(2,1,1);
histogram(magnitude, 'Normalization', 'pdf');
hold on;
x = linspace(0, max(magnitude), 100);
sigma = 1;  % Rayleigh parameter
rayleigh_pdf = (x ./ sigma^2) .* exp(-(x.^2) / (2 * sigma^2));
plot(x, rayleigh_pdf, 'r', 'LineWidth', 2);
title('Magnitude of h (Rayleigh Distribution)');
xlabel('Magnitude');
ylabel('Probability Density');
legend('Empirical Magnitude', 'Theoretical Rayleigh PDF');
grid on;

subplot(2,1,2);
histogram(phase, 'Normalization', 'pdf');
hold on;
phase_theoretical = ones(1, 100) / (2 * pi);  % Uniform PDF
plot(linspace(0, 2*pi, 100), phase_theoretical, 'r', 'LineWidth', 2);
title('Phase of h (Uniform Distribution)');
xlabel('Phase (radians)');
ylabel('Probability Density');
legend('Empirical Phase', 'Theoretical Uniform PDF');
grid on;

% Plot equalized signal magnitude and phase
figure;
subplot(2,1,1);
plot(abs(equalized_signal), 'b');
title('Magnitude of Equalized Signal');
xlabel('Symbol Index');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(angle(equalized_signal), 'r');
title('Phase of Equalized Signal');
xlabel('Symbol Index');
ylabel('Phase (radians)');
grid on;
